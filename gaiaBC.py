#!/usr/bin/env python
import sys
from math import *
import random
import numpy as np
from scipy.linalg import expm
from astropy.io import fits

pi=np.pi

def FINDCLOSE(par1,pararr):
    nclose=np.argsort(abs(np.add(-par1,pararr)))
    return nclose[0]

def rot_euler(v, xyz):
    ''' Rotate vector v (or array of vectors) by the euler angles xyz '''
    # https://stackoverflow.com/questions/6802577/python-rotation-of-3d-vector
    for theta, axis in zip(xyz, np.eye(3)):
        v = np.dot(np.array(v), expm(np.cross(np.eye(3), axis*-theta)))
    return v


def BCcal(wavelength,flux,lambdaF,weight,AVstar,Rv,Teff,par2vega,DrainearrLam,DrainearrK): # sed, filter
#  BC= Mbol_sun -2.5 log10 [4!PI 10pc^2 Sigma Teff^4 / Lsun] + 2.5 .... Girardi + 2002
    Mbolsun = 4.77
    bolconstant = -2.5*log10(4*pi*(3.0857**2.0)*5.67051/3.826)
 
    n_wave=len(wavelength)
# I put the zero weight at the edges just to solve the problem of edges in linear interpolation
    weight[0]=0.0
    weight[len(weight)-1]=0.0
    wavelength_weight=np.interp(wavelength,lambdaF,weight)
    par1=0.0
    par2=0.0
    alamarr=[] #fltarr((size(wavelength))[-1])

    kappased=np.interp(wavelength,DrainearrLam,DrainearrK)
    alamarr=AVstar*kappased

  
    for i in range(1, len(wavelength)):
      ltemp=wavelength[i]*1.0e-4 
      par1 += wavelength[i]*flux[i]*(10.0**(-0.4*alamarr[i]))*wavelength_weight[i]*(wavelength[i]-wavelength[i-1])
#     par2 += wavelength[i]*wavelength_weight[i]*(wavelength[i]-wavelength[i-1])  ;!!!! par2 will be calculating from Vega flux
    BCfilter = Mbolsun + bolconstant - 10.*log10(Teff)+2.5*log10(par1/par2vega)
    return BCfilter

def makeGaussian(size, fwhm, center):
    """ Make a square gaussian kernel.
    size is the length-array of a side of the square
    fwhm is full-width-half-maximum, which
    can be thought of as an effective radius.
    """
    x = np.arange(0, size[0], 1, float)
    y = np.arange(0, size[1], 1, float)
    y = y[:,np.newaxis]
    if center is None:
        x0 = y0 = size // 2
    else:
        x0 = center[0]
        y0 = center[1]
    return np.exp(-4*np.log(2) * ((x-x0)**2 + (y-y0)**2) / fwhm**2)



for a in sys.argv[1:]:
    spa = a.strip().split('=')
    if len(spa) == 2:
        if spa[0] == 'sim': simfolder = str(spa[1])
        elif spa[0] == 'age': agestarinput = float(spa[1])
        elif spa[0] == 'AV': AVstar = float(spa[1])
        elif spa[0] == 'fov': fov = float(spa[1])
        elif spa[0] == 'distance': distance = float(spa[1])
        else: print 'unknown variable'
    else: print 'unknown argument'


vegafile='vegaf'
filterfileG='Gaia/GAIA0.G.dat'
filterfileGB='Gaia/GAIA0.Gbp.dat'
filterfileGR='Gaia/GAIA0.Grp.dat'

#simfolder='/Users/Zeinab/Desktop/nbody6-zeinab/zeinab-sims/M103/Q05/Rh05/D20S00bin00M103r05Q05/'
#agestarinput=2.0 #[Myr] between 0-100
kzstarinput=0.015
rhostarinput=0.0 #[Msun/pc2]
alphai=0.0
bettai=0.0
gammai=0.0
res=0.01

fovx=fov
fovy=fov
Rv=3.1


filecloud= 'NoCloud'
Columndensities='user'
OBtreatment='yes'


evolutionary='Z0p015.dat'
foldersed='SEDs/'
Mbolsun = 4.77

filestar=simfolder+simfolder[-22:-1]+'/'+simfolder[-22:-1]+'-'+np.str("%05g" %(round(agestarinput)))+'-sin'

print 'ID'.center(5),'RA'.center(13),'Dec'.center(13),'Parallax'.center(13),'Proper motion'.center(13),'Proper motion'.center(13),'G mag'.center(13),'(BP - RP)'.center(13),'Radial'.center(13),'Mass'.center(13),'Teff'.center(13),'logL/Lo'.center(13)
print 'star'.center(5),'[arcsec]'.center(13),'[arcsec]'.center(13),'[mas]'.center(13),'in RA[mas/yr]'.center(13) ,'in Dec[mas/yr]'.center(13), ' '.center(13), '[mag]'.center(13),'Velocity[km/s]'.center(13),'[Mo]'.center(13),'[K]'.center(13),' '


lambdaFG,weightG=np.loadtxt(filterfileG,unpack=True)
lambdaFGB,weightGB=np.loadtxt(filterfileGB,unpack=True)
lambdaFGR,weightGR=np.loadtxt(filterfileGR,unpack=True)

if (Rv == 3.1): 
    Drainemodel='Draine3p1.txt'
    DraineKappaV=8.551E+03
elif (Rv == 4.0): 
    Drainemodel='Draine4.txt'
    DraineKappaV=8.492E+03
elif (Rv == 5.5): 
    Drainemodel='Draine5p5.txt'
    DraineKappaV=7.313E+03
else: print 'For Dmodel, R_V should be 3.1 or 4.0 or 5.5. If you need other Rv values please choose Fmodel'
DrainearrLamu,drainealbedo,drainecos,draineC,DrainearrKu,drainecos2=np.loadtxt(Drainemodel,usecols=(0,1,2,3,4,5),unpack=True)
DrainearrLamu=DrainearrLamu*1.0E+4
DrainearrKu=DrainearrKu/DraineKappaV
DrainearrLam, DrainearrK = zip(*sorted(zip(DrainearrLamu,DrainearrKu)))

xpix=round(fovx/res)
ypix=round(fovy/res)

wavelengthvega,fluxvega=np.loadtxt(vegafile,unpack=True)
wavelength_weightvegaG=np.interp(wavelengthvega,lambdaFG,weightG)
wavelength_weightvegaGB=np.interp(wavelengthvega,lambdaFGB,weightGB)
wavelength_weightvegaGR=np.interp(wavelengthvega,lambdaFGR,weightGR)
par2vegaG=0.0
par2vegaGB=0.0
par2vegaGR=0.0

for i in range(1, (len(wavelengthvega))):
  par2vegaG += fluxvega[i]*wavelengthvega[i]*wavelength_weightvegaG[i]*(wavelengthvega[i]-wavelengthvega[i-1])
  par2vegaGB += fluxvega[i]*wavelengthvega[i]*wavelength_weightvegaGB[i]*(wavelengthvega[i]-wavelengthvega[i-1])
  par2vegaGR += fluxvega[i]*wavelengthvega[i]*wavelength_weightvegaGR[i]*(wavelengthvega[i]-wavelengthvega[i-1])

#print par2vegaG,par2vegaGB,par2vegaGR

IDstar,xstar,ystar,zstar,vxstar,vystar,vzstar,massstar=np.loadtxt(filestar,usecols=(0,1,2,3,4,5,6,10),unpack=True)
nstar=len(xstar)
logagestar=np.full(nstar,log10(agestarinput*1.0E+6))
kzstar=np.full(nstar,kzstarinput)
rhostar=np.full(nstar,rhostarinput)

positionvector=zip(xstar,ystar,zstar)
velocityvector=zip(vxstar,vystar,vzstar)

zeinab=rot_euler(positionvector,np.multiply([alphai,bettai,gammai],pi/180.))
xstar=zeinab[0:nstar,0]
ystar=zeinab[0:nstar,1]
zstar=zeinab[0:nstar,2]

zeinabv=rot_euler(velocityvector,np.multiply([alphai,bettai,gammai],pi/180.))
vxstar=zeinabv[0:nstar,0]
vystar=zeinabv[0:nstar,1]
vzstar=zeinabv[0:nstar,2]

distancestar=np.add(distance,-zstar)
pc2pixstar=206264.806247/distancestar/res
Teffstar=np.zeros(nstar)
loggstar=np.zeros(nstar)
loglstar=np.zeros(nstar)
sedstar=np.zeros(nstar,dtype=np.uint64)
fluxstar=np.zeros(nstar)
newx=xstar*pc2pixstar #convert x[pc] into pixel position
newy=ystar*pc2pixstar
newz=zstar*pc2pixstar
columncloud=np.zeros(nstar) #column density of the cloud in fron of each star

if (Columndensities == 'sph'):
#    reading the cloud:
    xcloud,ycloud,zcloud,vxcloud,vycloud,vzcloud,masspar,hpar=np.loadtxt(filecloud,unpack=True)
    ncloud=len(xcloud)
    
    positioncvector=zip(xcloud,ycloud,zcloud)
    zeinabc=rot_euler(positioncvector,np.multiply([alphai,bettai,gammai],pi/180.))
    xcloud=zeinabc[0:nstar,0]
    ycloud=zeinabc[0:nstar,1]
    zcloud=zeinabc[0:nstar,2]

    distancecloud=np.add(distance,-zcloud)
    pc2pixcloud=206264.806247/distancecloud/res
    newxcloud=xcloud*pc2pixcloud #convert x[pc] into pixel position
    newycloud=ycloud*pc2pixcloud
    newzcloud=zcloud*pc2pixcloud
    newhcloud=hpar*pc2pixcloud

#  reading the isochrones
ziso,logageiso,miniiso,mactiso,logliso,logteff,loggiso=np.loadtxt(evolutionary,usecols=(0,1,2,3,4,5,6),unpack=True)
teffiso=10.**logteff
niso=len(miniiso)


#  reading list of SEDs:
teffsed,loggsed,metallicity,lh,vtur,sedname=np.loadtxt(foldersed+'kseds.dat',dtype={'names': ('col1', 'col2', 'col3','col4', 'col5', 'col6'), 'formats':(np.float,np.float,np.float,np.float,np.float,'|S60')},usecols=(0,1,2,3,4,5),unpack=True)
nseds=len(teffsed)

if (OBtreatment == 'yes'):
    teffsedOB,loggsedOB,metallicityOB,sednameOB=np.loadtxt(foldersed+'Tseds.dat',dtype={'names': ('col1', 'col2', 'col3','col4'), 'formats':(np.float,np.float,np.float,'|S60')},usecols=(0,1,2,3),unpack=True)
    nsedsOB=len(teffsedOB)

nstars=0
for ii in range(nstar):
    if ((abs(newx[ii]) < (xpix/2)-1) and (abs(newy[ii]) < (ypix/2)-1)):
        nstars += 1

#        cloudden=rhostar[ii]*((1.989/1.673534)/((3.0857**2))) #*1.0e21 #convert column density unit [Msun/pc^2] --> [10^21 * Mhydrogen/cm^2]
#        AVstar=cloudden/2.21 #e21 ;Guver&Ozel2009: The relation between Optical Extinction and Hydrogen column density
        nage=FINDCLOSE(logagestar[ii],logageiso)
        selectedage=logageiso[nage]
        nmetalicity=FINDCLOSE(kzstar[ii],ziso)
        selectedz=ziso[nmetalicity]
        marrtemp=np.full(niso,999.99)
        for kk in range(niso):  
            if ((ziso[kk] == selectedz) and (logageiso[kk] == selectedage)):  marrtemp[kk]=mactiso[kk]
        ns=FINDCLOSE(massstar[ii],marrtemp)
        Teffstar[ii]=teffiso[ns]
        loggstar[ii]=loggiso[ns]
        loglstar[ii]=logliso[ns]

        deltaT=abs(Teffstar[ii]-teffsed)
        deltagarr=np.full(nseds,99.)


        for jj in range(nseds): 
            if (deltaT[jj] == min(deltaT)):  deltagarr[jj]=abs(loggstar[ii]-loggsed[jj])
        sedstar[ii]=FINDCLOSE(min(deltagarr),deltagarr)
        readsed=sedname[sedstar[ii]]
        wavelength,flux=np.loadtxt(foldersed+sedname[sedstar[ii]],comments=['fn:', '#'],unpack=True)

        if ((OBtreatment == 'yes') and (Teffstar[ii] >= 15000.)):
            deltaT=abs(Teffstar[ii]-teffsedOB)
            deltagarr=np.full(nsedsOB,99.)
            for jj in range(nsedsOB):  
                if (deltaT[jj] == min(deltaT)):  deltagarr[jj]=abs(loggstar[ii]-loggsedOB[jj])
            sedstar[ii]=FINDCLOSE(min(deltagarr),deltagarr)

            readsed=sednameOB[sedstar[ii]]
            wavelength,flux=np.loadtxt(foldersed+sednameOB[sedstar[ii]],comments=['fn:', '#'],unpack=True)

        bcG=BCcal(wavelength,flux,lambdaFG,weightG,AVstar,Rv,Teffstar[ii],par2vegaG,DrainearrLam,DrainearrK)
        magG=Mbolsun-2.5*loglstar[ii]-(bcG)+5.0*log10(distancestar[ii]/10.0)

        bcGB=BCcal(wavelength,flux,lambdaFGB,weightGB,AVstar,Rv,Teffstar[ii],par2vegaGB,DrainearrLam,DrainearrK)
        magGB=Mbolsun-2.5*loglstar[ii]-(bcGB)+5.0*log10(distancestar[ii]/10.0)
        
        bcGR=BCcal(wavelength,flux,lambdaFGR,weightGR,AVstar,Rv,Teffstar[ii],par2vegaGR,DrainearrLam,DrainearrK)
        magGR=Mbolsun-2.5*loglstar[ii]-(bcGR)+5.0*log10(distancestar[ii]/10.0)


        print("%5g %13.4f %13.4f  %13.4f %13.4f %13.4f %13.4f %13.4f %13.4f %13.2f %13.2f %13.4f  " %(IDstar[ii],newx[ii]*res,newy[ii]*res,1000.0/distancestar[ii],vxstar[ii]*pc2pixstar[ii]*res*0.001,vystar[ii]*pc2pixstar[ii]*res*0.001, magG, magGB-magGR, vzstar[ii],massstar[ii],Teffstar[ii],loglstar[ii]))
  
