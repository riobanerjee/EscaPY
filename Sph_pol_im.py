# PROGRAM NAME : Sph_im_pol
# AUTHOR : Agnibha Banerjee
# DATE : 5 Aug 2018
#
# This is a program to simulate the movement of photons
# (only scattering an absorption) through a spherically
# symmetric stellar atmosphere, using a 1-D opacity grid
# in the radial direction. The output is the final
# intensities at different locations on the surface
# of the star projected on a 2-D plane.
# (INCLUDES POLARISATION)

from math import *
from random import *
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


# All the parameters are inititalized below. The
# density goes as r**(rhoexp). The final image is
# imagenx x imagenx pixels in size.
# limb=0 ignores limb darkening law while limb=1
# implements it. The last 4 are required for 
# creating the polarisation matrices.
# ngrid is the number of grid cells created.
# The rest of the variable names are
# self explanatory.

seed(76542)
albedo=1.0	
nphoton=100000
taumax=10.
rmin=1.	
rmax=10.
limb=0.
rhoexp=-1.
ngrid=10
imagenx=20

a=0.75
g=1.0
pl=0
pc=0

r=np.linspace(rmin,rmax,ngrid)	 #grid of shells in r		


# newphot emits a new photon isotropically
# from the origin into the hemisphere above it.
# 		Parameters : None
# 		Returns : 
# mu : Cosine of angle between direction of emission and +z axis
# phi : Azimuthal angle on XY plane from X axis
# xt,yt,zt : Cartesian Co-ordinates


def newphot():							
	mu=sqrt(random())
	phi=2*pi*random()
	xt=yt=zt=0.
	return (mu,phi,xt,yt,zt)


# isoscatt isotropically scatters a photon in all directions.
# 		Parameters : None
# 		Returns :
# mu : Cosine of angle between direction of emission and +z axis
# phi : Azimuthal angle on XY plane from X axis


def isoscatt():			
	mu=2*random()-1
	phi=2*pi*random()
	return (mu,phi)


# lmatrix creates the rotation matrix:
#       | 1  		0 			0 		 0 |
#		| 0 	cos(2*ang)	sin(2*ang)   0 |
#		| 0 	-sin(2*ang) cos(2*ang)	 0 |
#		| 0		    0		    0		 1 |
#
#		Parameters :
# ang : The angle used in the rotation matrix
# 		Returns :
# l : The rotation matrix


def lmatrix(ang):
	l=np.zeros([4,4])
	l[0,0]=l[3,3]=1.
	l[1,1]=l[2,2]=cos(2*ang)
	l[1,2]=sin(2*ang)
	l[2,1]=-l[1,2]
	l=a*l
	return l


# rmatrix_electron creates the polarisation 
# matrix for Thomson(electron) scattering.
#		| p1  p2  0  0 |
#		| p2  p1  0  0 |
#	 	| 0   0  p3 -p4|
#		| 0   0  p4  p3|
# where, p1 = cosang**2+1
#		 p2 = cosang**2-1
#		 p3 = 2*cosang
#		 p4 = 0
#		
#		Parameters :
# cosang : Cosine of the angle of scattering
# 		Returns :
# r : The polarisation matrix


def rmatrix_electron(cosang):
	r=np.zeros([4,4])
	r[0,0]=r[1,1]=cosang**2+1
	r[0,1]=r[1,0]=cosang**2-1
	r[2,2]=r[3,3]=2*cosang
	r=a*r
	return r


# rmatrix_dust creates the polarisation matrix
# for dust scattering.
# This function is still under construction and
# not to be used.


def rmatrix_dust(cosang):
	r=np.zeros([4,4])
	r[0,0]=r[1,1]=(1-g**2)/((1+g**2-2*g*cosang)**1.5)
	r[0,1]=r[1,0]=-pl*r[0,0]*(1-cosang**2)/(1+cosang**2)
	r[2,2]=r[3,3]=r[0,0]*2*cosang/(1+cosang**2)
	r[2,3]=pc*r[0,0]*(1-cosang**2)/(1+cosang**2)
	r[3,2]=-r[2,3]
	r=a*r
	return r


# raylscatt scatters a photon using Rayleigh Scattering.
# The final Stokes' parameters(S) is obtained using :
# 	S = L(pi-i2)RL(-i1)S'
#		Parameters :
# muini : Initial mu
# phini : Initital phi
# sinit : Initial Stokes' parameters
#		Returns  :
# mufin : mu after scattering
# phifin : phi after scattering
# sfinal : Stokes' parameters after scattering


def raylscatt(muini,phini,sinit):
	sfinal=np.zeros(4)
	i1=i2=mufin=phifin=0.0
	cosT=0.0
	while (True):

		i1=2*pi*random()
		l1=lmatrix(-i1)
		cosT=1-2*random()
		mufin=sqrt(1-muini**2)*sqrt(1-cosT**2)*cos(i1)+muini*cosT
		p=sqrt(1-muini**2)*sin(i1)/sqrt(1-mufin**2)
		if (p>1):
			i2=pi/2
		elif(p<-1):
			i2=-pi/2
		else:
			i2=asin(p)
		t=sin(i1)*sqrt(1-cosT**2)/sqrt(1-mufin**2)
		if (t>1):
			phifin=phini-pi/2
		elif (t<-1):
			phifin=phini+pi/2
		else:
			phifin=phini-asin(t)
		l2=lmatrix(pi-i2)
		r=rmatrix_electron(cosT)
		sfinal=np.matmul(np.matmul(l2,np.matmul(r,l1)),sinit)
		if(random()<sfinal[0]/2):
			break
	sfinal=sfinal/sfinal[0]
	return(mufin,phifin,sfinal)


# locate searches an ordered (in ascending order)
# list using binary search technique.
# 		Parameters :
# a : The ordered list to be searched
# s : The value to be located
# 		Returns :
# left : a[left] <= s <= a[right] (The left index of the cell
#        in which the value s is)


def locate(a,s):
	left=0
	right=len(a)-1
	while(a[right]-1>a[left]):
		mid=(left+right)/2
		if (a[mid]>s):
			right=mid
		else:
			left=mid
	return left


# darkening is a function to apply limb-darkening law.
# 		Parameters : None
# 		Returns :
# mu : Cosine of angle between direction of emission and +z axis


def darkening():
	ntab=110
	dmu=1.0/ntab
	xmu=np.zeros(ntab)
	prob=np.zeros(ntab)
	for i in range(1,ntab):
		xmu[i]=dmu*i
		prob[i]=0.5*(xmu[i]**3+xmu[i]**2)  #probability distribution table of mu
	rand=random()
	n=locate(prob,rand)
	mu=xmu[n]-(prob[n]-rand)*dmu/(prob[n]-prob[n-1])
	return mu


# newphotsph emits a new photon at the surface of a sphere
# of radius rmin with random direction pointing outwards.
# 		Parameters : None
# 		Returns :
# costi : Cosine of angle between direction of emission and +z axis
# phi : Azimuthal angle on XY plane from X axis
# xt,yt,zt : Cartesian Co-ordinates


def newphotsph():
	mu,phi,xt,yt,zt=newphot()
	cosb=1-2*random()
	sinb=sqrt(1-cosb**2)
	if limb!=0:	
		mu=darkening()
	lp=2*pi*random()
	if(abs(cosb)<0.9999999):
		costi=mu*cosb-sqrt(1-mu**2)*sinb*cos(phi)
		if(costi>1):
			costi=1.0
		elif(costi<-1):
			costi=-1.0
		sinti=sqrt(1-costi**2)
		if(sinti>0.0000001):
			cospn=(mu-costi*cosb)/(sinti*sinb)
			if(cospn>1.01):
				cospn=1.0
			elif(cospn<-1.01):
				cospn=-1.0
			if(phi<pi):
				phinew=acos(cospn)
				if(phinew>pi):
					phinew=pi
			else:
				phinew=2*pi-acos(cospn)
				if(phinew<pi):
					phinew=pi
		else:
			phinew=0.
		phi=phinew
	else:
		if(cosb<0):
			costi=-mu
	phi=phi+lp
	if(phi>2*pi):
		phi=phi-2*pi
	zt=cosb*rmin*1.0001
	xt=sinb*cos(lp)*rmin*1.0001
	yt=sinb*sin(lp)*rmin*1.0001
	return(costi,phi,xt,yt,zt)


# makegrid makes a 1-D grid of opacity and optical depth
# as a function of radius.
# 		Parameters : None
# 		Returns : 
# opac : Grid of opacity
# tauint : Grid of optical depth


def makegrid():

    
    dr=(rmax-rmin)/(ngrid-1)
    opac=[0]*ngrid
    tauint=[0]*ngrid
    
    if (rhoexp==-1):
    	kmax=taumax/(log(rmax)-log(rmin))
    else:
    	kmax=taumax*(1+rhoexp)/((rmax**(1.+rhoexp))-(rmin**(1.+rhoexp)))
    
    rtau=0.0

    
    for i in range(ngrid-1):
    	if (rhoexp==-1):
    		opac[i]=kmax*(log(r[i+1])-log(r[i]))/dr
    	else:
    		opac[i]=kmax*((r[i+1])**(1+rhoexp)-(r[i])**(1+rhoexp))/((1+rhoexp)*dr)
    	rtau=rtau+opac[i]*dr
    	tauint[i+1]=rtau
    return opac,tauint

opac,tauint=makegrid()  #Making opac and tauint globally accessible


# findWall finds out which face of the grid cell (outer or inner) is closest
# to the photon along it's direction of motion.
# 		Parameters :
# r1 : Radius of inner shell
# r2 : Radius of outer shell
# nx,ny,mu : Direction cosines in the x, y and z directions
# xp,yp,zp : Cartesian Co-ordinates of the photon
# 		Returns:
# t : Distance to the nearest face
# ioffset : 1 if outer shell is closest, -1 if inner shell is closest


def findWall(r1,r2,nx,ny,mu,xp,yp,zp):
	root=[-999,-999,-999,-999]
	ind=[0,0,1,1]
	bb=(xp*nx+yp*ny+zp*mu)						
	cc=xp*xp+yp*yp+zp*zp-r1*r1
	dd1=bb*bb-cc
	if (dd1>=0):  		#finds distance to inner face
		dd1=sqrt(dd1)
		root[0]=(-bb+dd1)
		root[1]=(-bb-dd1)
	cc=xp*xp+yp*yp+zp*zp-r2*r2
	dd2=bb*bb-cc
	if (dd2>=0):		#finds distance to outer face
		dd2=sqrt(dd2)
		root[2]=(-bb+dd2)
		root[3]=(-bb-dd2)
	npos=0
	posroot=[0.,0.,0.,0.]
	for i in range(4):
		if(root[i]>0):
			posroot[npos]=root[i]
			npos+=1
	if (npos==0):
		t=0.00001
		print "But Why??",r1,r2,sqrt(bb+r2*r2),nx,ny,mu,sqrt(nx*nx+ny*ny+mu*mu)	
    		ioffset=1
        else:			#only positive roots are accepted
		if (npos==1):
            		t = posroot[0]
        	elif (npos==2):
            		t = min(posroot[0],posroot[1])
        	elif (npos==3):
       	    		t = min(posroot[0],posroot[1],posroot[2])
        	else:
            		t = min(posroot[0],posroot[1],posroot[2],posroot[3])
		for i in range(4):
			if(root[i]==t):
				iface=ind[i]
				ioffset=-1+(2*iface)
		
	
	return(t,ioffset)


# tau_1d updates the position and direction of the photon
# and the intensity moments until the sampled optical depth 
# is reached.
# 		Parameters :
# tau : The sampled optical depth
# xold,yold,zold : Cartesian co-ordinates of photon
# nx,ny,mu : Direction cosines in x, y and z directions
# ii : The grid cell number in which the photon is
# exitflag : 1 if the photon has escaped, 0 if not
# jp,hp,kp : Intensity moments in outwards (+) direction
# jm,hm,km : Intensity moments in inwards (-) direction
# aflag : 1 if the photon has been absorbed (r<rmin), 0 if not
# 		Returns:
# exitflag : 1 if the photon has escaped, 0 if not
# xp,yp,zp : Updated cartesian co-ordinates of photon
# ii : Updated grid cell number in the photon is now
# jp,hp,kp : Intensity moments in outwards (+) direction
# jm,hm,km : Intensity moments in inwards (-) direction
# aflag : 1 if the photon has been absorbed (r<rmin), 0 if not


def tau_1d(tau,xold,yold,zold,nx,ny,mu,ii,exitflag,jp,hp,kp,jm,hm,km,aflag):
	tsum=0.0
	xp=xold
	yp=yold
	zp=zold
	while(True):
		ds,ioffset=findWall(r[ii],r[ii+1],nx,ny,mu,xp,yp,zp)

		if(ds<1e-5):
			ds=1e-5

		opaci=opac[ii]
		dtau=opaci*ds

		if(tsum+dtau)>=tau:
			ds=(tau-tsum)/opaci
			xp=xp+ds*nx
			yp=yp+ds*ny
			zp=zp+ds*mu
			return(exitflag,xp,yp,zp,ii,jp,hp,kp,jm,hm,km,aflag)
		
		ii=ii+ioffset

		xp=xp+ds*nx*1.01           #For numerical precision
   		yp=yp+ds*ny*1.01
   		zp=zp+ds*mu*1.01
   		tsum=tsum+dtau

		if(ii==ngrid-1):
			exitflag=1
			costii=(xp*nx+yp*ny+zp*mu)/rmax
			if (costii>0):
				jp[ii]=jp[ii]+1./costii
				hp[ii]=hp[ii]+1.
				kp[ii]=kp[ii]+costii
			return(exitflag,xp,yp,zp,ii,jp,hp,kp,jm,hm,km,aflag)

		rtot=sqrt(xp*xp+yp*yp+zp*zp)
		iinew=locate(r,rtot)

		if(iinew==ii):
			if(ii>=0):
				costii=(xp*nx+yp*ny+zp*mu)/rtot
				if(costii>0):
					jp[ii]=jp[ii]+1./costii
					hp[ii]=hp[ii]+1.
					kp[ii]=kp[ii]+costii
				else:
					jm[ii+1]=jm[ii+1]-1./costii
					hm[ii+1]=hm[ii+1]-1.
					km[ii+1]=km[ii+1]-costii
		else:
			ii=iinew
	
		
		if (rtot<rmin*1.001):
			#aflux+=1
			costii=(xp*nx+yp*ny+zp*mu)/rtot
			jm[0]-1./costii
			hm[0]=hm[0]-1.
			km[0]=km[0]-costii
			aflag=1
			return (exitflag,xp,yp,zp,ii,jp,hp,kp,jm,hm,km,aflag)


# transfer simulates the whole radiative transfer process for
# all the photons.
# 		Parameters :
# taumax : The maximum optical depth of the star
# 		Returns :
# energy : Number of photons that have escaped from the surface
# jp,hp,kp : Intensity moments in outwards (+) direction
# jm,hm,km : Intensity moments in inwards (-) direction
# image : Pixel intenstity counts for the output image


def transfer(taumax):
	
	i=1
	energy=0
	fractx=2.*rmax/float(imagenx)
	jp=np.zeros(ngrid)
	hp=np.zeros(ngrid)
	kp=np.zeros(ngrid)
	jm=np.zeros(ngrid)
	hm=np.zeros(ngrid)
	km=np.zeros(ngrid)	
	image=np.zeros([imagenx,imagenx])
	poli=np.zeros([imagenx,imagenx])
	polq=np.zeros([imagenx,imagenx])
	polu=np.zeros([imagenx,imagenx])
	polv=np.zeros([imagenx,imagenx])
	while(i<=nphoton):
		if (i%1000==0):
			print i,"\tphotons completed."
		mu,phi,xt,yt,zt=newphotsph()
		exitflag=0
		ii=0
		xold=yold=zold=0.
		jp[ii]=jp[ii]+1./mu
		hp[ii]=hp[ii]+1.
		kp[ii]=kp[ii]+mu
		s=np.array([1,0,0,0])
		while (True) :
		
			nx=sqrt(1-mu**2)*cos(phi)
			ny=sqrt(1-mu**2)*sin(phi)
			aflag=0
			xold=xt
			yold=yt
			zold=zt
			tau=-log(random())    #Sampled from exponential distribution
			exitflag,xt,yt,zt,ii,jp,hp,kp,jm,hm,km,aflag=tau_1d(tau,xold,yold,zold,nx,ny,mu,ii,exitflag,jp,hp,kp,jm,hm,km,aflag)
			
			if (aflag==1):
				i+=1
				break			
			if (exitflag==1) :
				i+=1
				energy+=1
				x=yold*cos(phi)-xold*sin(phi)
				y=zold*sqrt(1-mu**2)-yold*mu*sin(phi)-xold*mu*cos(phi)
				ix=int((x+rmax)/fractx)
				iy=int((y+rmax)/fractx)
				if (ix<imagenx) and (iy<imagenx) and (ix>=0) and (iy>=0):
					image[ix,iy]=image[ix,iy]+1.	#creates the 2-D projection image with pixel intensities
					poli[ix,iy]=poli[ix,iy]+s[0]
					polq[ix,iy]=polq[ix,iy]+s[1]
					polu[ix,iy]=polu[ix,iy]+s[2]
					polv[ix,iy]=polv[ix,iy]+s[3]
				break
	
			if (random()<albedo):
				mu,phi,s=raylscatt(mu,phi,s)	#photon scattered when not absorbed	
			else:
				#aflux+=1
				i+=1

	return energy,jp,hp,kp,jm,hm,km,image,poli,polq,polu,polv


# output plots the final data that has been computed. 
# Can be altered to plot the following: 
# 1) The image of pixel intensity(projection).
# 2) The resultant moments vs optical depth.
# 3) The polarisation angles on the final image.
# 		Parameters : None
# 		Returns : None


def output():
	en,jp,hp,kp,jm,hm,km,image,poli,polq,polu,polv=transfer(taumax)
	#angle=np.zeros([imagenx/5,imagenx/5])        #Averaging 5 X 5 squares
	#for i in range(imagenx/5):
	#	for j in range(imagenx/5):
	#		for k in range(5*i,5*(i+1)):
	#			for l in range (5*j,5*(j+1)):
	#				angle[i,j]=angle[i,j]+0.5*np.arctan(polu[k,l]/polq[k,l])
	#angle=angle/25.
	#print angle.shape
	np.seterr(divide='ignore',invalid='ignore') #ignoring all pixels with no photons
	angle=0.5*np.arctan(polu/polq)
	deg=np.sqrt(polq**2+polu**2+polv**2)/poli
	#deg=deg/image
	energy= float(en)/nphoton	#Normalising all the outputs
	jp=jp/nphoton
	hp=hp/nphoton
	kp=kp/nphoton
	jm=jm/nphoton
	hm=hm/nphoton
	km=km/nphoton
	image=image/nphoton
	plt.imshow(image,'gray')
	plt.colorbar()
	#plt.plot(image[imagenx/2])
	f=open("moments.dat","w")
	for i in range(ngrid):
			tauint[i]=taumax-tauint[i]
			print>>f,r[i],tauint[i],jp[i],jm[i],hp[i],hm[i],kp[i],km[i]
	print "Energy = ",energy
	#plt.plot(tauint,(jp+jm))
	x,y=np.meshgrid(np.linspace(0,imagenx-1,imagenx),np.linspace(0,imagenx-1,imagenx))
	u=deg*np.cos(angle)
	v=deg*np.sin(angle)
	#plt.quiver(x[::10,::10],y[::10,::10],u[::10,::10],v[::10,::10],color='y')
	plt.quiver(x,y,u,v,color='b')
	plt.show()

output()
