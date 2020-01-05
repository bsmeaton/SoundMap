# -*- coding: utf-8 -*-
"""
Created on Fri Aug  4 11:31:20 2017

@author: BenjaSmea
"""
from __future__ import division 
import numpy as np
import sys
# Sound Attenuation calculator to ISO 9613-2 1996
sourcelist = np.atleast_2d(np.genfromtxt("SourceList.txt", delimiter=",", dtype=str,comments='#'))
walllist = np.atleast_2d(np.genfromtxt("WallList.txt", delimiter=",", dtype=str,comments='#'))
coordlist = np.atleast_2d(np.genfromtxt("inputCoordList.txt", delimiter=",", dtype=float,comments='#'))

#wall programsgvb 
def line(p1, p2):
    A = (p1[1] - p2[1])
    B = (p2[0] - p1[0])
    C = (p1[0]*p2[1] - p2[0]*p1[1])
    return A, B, -C

def intersection(L1, L2):
    D  = L1[0] * L2[1] - L1[1] * L2[0]
    Dx = L1[2] * L2[1] - L1[1] * L2[2]
    Dy = L1[0] * L2[2] - L1[2] * L2[0]
    if D != 0:
        x = Dx / D
        y = Dy / D
        return x,y
    else:
        return False
#SYSTEM VARIABLES
#Reciever height
zcoord = 1.5
hum = 0.7
T = 293.15
T0 = 293.15

# KEYOPS
speccoords = 0
outputdata = 0
Aweightingkey = 1
# calculation with A weighting
atmosatten = np.array((0.1,0.3,1,3.1,7.4,12.7,23.1,59.3))
freqspec = np.array((63,125,250,500,1000,2000,4000,8000))
Aweighting = np.array((-26.22,-16.19,-8.67,-3.25,0.00,1.20,0.96,-1.14))

#Point field set out
xminus = -500
xplus = 500
yminus = -500
yplus = 500
xspacing = 10
yspacing = 10

xc = xminus
yc = yminus
if speccoords == 0:
    coordlist = np.array((xc,yc,float(0)))
    while (xc <= xplus):
        yc = yminus
        while (yc <= yplus):
            coords = np.array((xc,yc,float(0)))
            coordlist = np.vstack((coordlist,coords))
            yc += yspacing
        xc += xspacing
def NoiseCalc(sourcelist):
    noiselist = []
    for a in range(0,len(sourcelist)):
        # Decibil value of source spectrum
        dbsource = 10*np.log10(10**(float(sourcelist[a][5])/10)+10**(float(sourcelist[a][6])/10)+10**(float(sourcelist[a][7])/10)+10**(float(sourcelist[a][8])/10)+10**(float(sourcelist[a][9])/10)+10**(float(sourcelist[a][10])/10)+10**(float(sourcelist[a][11])/10)+10**(float(sourcelist[a][12])/10))
        # Source to reciever distance
        lindist = np.sqrt((float(sourcelist[a][2])-xcoord)**2+(float(sourcelist[a][3])-ycoord)**2+(float(sourcelist[a][4])-zcoord)**2)
        # Source to recieve distance along ground plane
        dp = np.sqrt((float(sourcelist[a][2])-xcoord)**2+(float(sourcelist[a][3])-ycoord)**2)
        # Average height
        h = (float(sourcelist[a][4])+zcoord)/2
        #Attenuation due to geometric divergence
        geoatten = max((20*np.log10(lindist)+11,0))
        #Attenuation due to atmospheric attenuation
        freq = np.array((63,125,250,500,1000,float(2000),4000,float(8000)))
        Frn = ((T/T0)**(-0.5))*(9+280*hum*np.exp(-4.17*(((T/T0)**(-1/3))-1)))
        Fro = 24+(4.04*10**4)*hum*((0.02+hum)/(0.391+hum))
        ac = (869*(freq**2))*((1.84*10**-11)*((T/T0)**0.5)+((T/T0)**(-5/2))*(0.01275*(np.exp(-2239.1/T)/((Fro+(freq**2/Fro))))+0.1068*(np.exp(-3352/T)/(Frn+(freq**2/Frn)))))
        atmos = (atmosatten*(lindist/1000))
        atmos2 = (ac*lindist)/100
        #Attenuation due to ground attenuation
        G=0
        if dp < 30*(float(sourcelist[a][4])+zcoord):
            q = 0
        else:
            q = 1-(30*(float(sourcelist[a][4])+zcoord)/dp)
        ah = 1.5+3*np.exp(-0.12*(h-5)**2)*(1-np.exp(-dp/50))+5.7*np.exp(-0.09*h**2)*(1-np.exp(-2.8*10**-6*dp**2))
        bh =1.5+8.6*np.exp(-0.09*h**2)*(1-np.exp(-dp/50))
        ch =1.5+14*np.exp(-0.46*h**2)*(1-np.exp(-dp/50))
        dh = 1.5+5*np.exp(-0.9*h**2)*(1-np.exp(-dp/50))
        As = np.array((-1.5,-1.5+G*ah,-1.5+G*bh,-1.5+G*ch,-1.5+G*dh,-1.5*(1-G),-1.5*(1-G),-1.5*(1-G)))
        Am = np.array((-3*q,-3*q*(1-G),3*q*(1-G),3*q*(1-G),3*q*(1-G),3*q*(1-G),3*q*(1-G),3*q*(1-G)))
        Agr = As+As+Am
        # Wall Attenuation
        Dz0 = np.array((0,0,0,0,0,0,0,0))
        for c in range(0,len(walllist)):
            #Points describing the line of the wall (Note: Current issue with entirely vertical or horizontal walls make slight slant)
            wallline = np.array([[float(walllist[c][2]),float(walllist[c][3])],[float(walllist[c][4]),float(walllist[c][5])]])
            #Points describing the source and reciever
            sourceline = np.array([[xcoord,ycoord],[float(sourcelist[a][2]),float(sourcelist[a][3])]])
            L1 = line(wallline[0], wallline[1])
            L2 = line(sourceline[0], sourceline[1])
            R = intersection(L1, L2)
            #check that source to reciever line and wall line intersect within bounds
            if R:

                # Test that source to reciever intersection points are within the bounds.
                if min(wallline[0][0],wallline[1][0]) <= R[0] <= max(wallline[0][0],wallline[1][0]) and min(wallline[0][1],wallline[1][1]) <= R[1] <= max(wallline[0][1],wallline[1][1]) and min(sourceline[0][0],sourceline[1][0]) <= R[0] <= max(sourceline[0][0],sourceline[1][0]):
                    # Barrier attenuation constants
                    #print("points = " + str(sourcelist[a][2]) + "," + str(sourcelist[a][3])+ "," + str(sourcelist[a][4]) + " Intersection" + str(R) + " X and Y" + str(L1) + "," + str(L2))
                    c2 = 20
                    c3 = 1
                    # Distance from source to wall top
                    dss = np.sqrt((float(sourcelist[a][2])-R[0])**2+(float(sourcelist[a][3])-R[1])**2+(float(sourcelist[a][4])-float(walllist[c][6]))**2)
                    # Distance from wall top to reciever
                    dsr = np.sqrt((xcoord-R[0])**2+(ycoord-R[1])**2+(zcoord-float(walllist[c][6]))**2)
                    #linear distance from source
                    sowall = np.sqrt((float(sourcelist[a][2])-R[0])**2+(float(sourcelist[a][3])-R[1])**2)
                    #source reciever line height at wall
                    linh = (np.abs(float(sourcelist[a][4])-zcoord)*(1-(sowall/dp)))+min((float(sourcelist[a][4])),zcoord)
                    z = dss+dsr-lindist
                    kmet = np.exp(-(1/2000)*np.sqrt((dss*dsr*lindist)/(2*z)))
                    l = 340/freqspec
                    if linh < float(walllist[c][6]):
                        Dz = 10*np.log10(3+(c2/l)*c3*z*kmet) + Agr
                    else:
                        Dz = np.array((0,0,0,0,0,0,0,0))
                    for d in range(0,len(Dz)):
                        if Dz[d] >20:
                           Dz[d] =20
                    Dz0 = np.vstack((Dz0,Dz))

                    dbwall = 10*np.log10(10**(float(Dz[0])/10)+10**(float(Dz[1])/10)+10**(float(Dz[2])/10)+10**(float(Dz[3])/10)+10**(float(Dz[4])/10)+10**(float(Dz[5])/10)+10**(float(Dz[6])/10)+10**(float(Dz[7])/10))

        #Find the maximum attenuation from all barrier in path
        maxDz = Dz0[np.argmax(np.mean(np.atleast_2d(Dz0), axis=1))]
        #Source spectrum array
        sourcearray = np.array(((float(sourcelist[a][5])),(float(sourcelist[a][6])),(float(sourcelist[a][7])),(float(sourcelist[a][8])),(float(sourcelist[a][9])),(float(sourcelist[a][10])),(float(sourcelist[a][11])),(float(sourcelist[a][12]))))
        # corrected sound power level 
        if outputdata == 1:
            print('Distance = ' + str(lindist))
        noisecorr = sourcearray - geoatten-atmos-Agr-maxDz
        noisecorrA = sourcearray - geoatten-atmos-Agr-maxDz + Aweighting
        if Aweightingkey == 1:
            noisecorrf = noisecorrA
        else:
            noisecorrf = noisecorr
        dbsource = 10*np.log10(10**(float(noisecorrf[0])/10)+10**(float(noisecorrf[1])/10)+10**(float(noisecorrf[2])/10)+10**(float(noisecorrf[3])/10)+10**(float(noisecorrf[4])/10)+10**(float(noisecorrf[5])/10)+10**(float(noisecorrf[6])/10)+10**(float(noisecorrf[7])/10))
        base10 = 10**(dbsource/10)
        noiselist = np.append(noiselist,base10)
    noiselevel = 10*np.log10(np.sum(noiselist))
    return noiselevel

def write(writestr,boundrow,boundhead):
    write = open(writestr, 'w')
    write.write(str(boundhead)+'\n')
    for a in range(0,len(boundrow)):
        for b in range(0,len(boundrow[a])):
            write.write(str(float(boundrow[a][b])))
            if b < len(boundrow[a])-1:
                write.write(',')
        write.write('\n')
    write.close()  
    
if speccoords == 0:
    print('Node Number: ' + str((np.abs(xminus)+np.abs(xplus)/xspacing)*(np.abs(yminus)+np.abs(yplus)/yspacing)))   
for b in range(0,len(coordlist)):  
    xcoord = coordlist[b][0]
    ycoord = coordlist[b][1]   
    nl = NoiseCalc(sourcelist)
    per = int(b/int(len(coordlist))*100)
    sys.stdout.write("\r%d%%" % per)
    sys.stdout.flush()
    coordlist[b][2]=nl

coordstring = 'coordlist.txt'
coordhead = 'CSV file'
print(coordlist)
write(coordstring,coordlist,coordhead)