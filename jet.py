# Parallel rays through uniform blob from
# under it.
# Changes to be made :
# 1) Sphere -> Cylinder   -- This is done.
# 2) Uniform -> Non-uniform (grids) -- This is done. With unit Cartesian grid.
# 3) From under it -> Anywhere -- This can probably be done by changing the mu initial.


from math import *
from random import *
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

seed(32783)

H=1.
albedo=1.0	
nphoton=1000000
taumax=1.
mubins=20 # Number of bins = mubins-1. mubins is the number of edges of bins
R=1.	  # Radius of Sphere
a=0.75	  # Required in the polarisation matrices
muspace=np.linspace(-1,1,mubins)
xmax=H/2
ymax=R
zmax=R
resx=10
resy=10
resz=10
xface=np.linspace(-xmax,xmax,resx)
yface=np.linspace(-ymax,ymax,resy)
zface=np.linspace(-zmax,zmax,resz)

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


# isoscatt isotropically scatters a photon in all directions.
#       Parameters : None
#       Returns : 
# mu : Cosine of angle between direction of emission and +Z-axis
# phi : Azimuthal angle on XY plane from +X-axis


def isoscatt():

	mu=2*random()-1
	phi=2*pi*random()
	return mu,phi


# newphot emits a new photon from the XY plane within a radius R
# and travelling parallel to the +Z axis. from below a cylinder lying along the x axis
#       Parameters : None
#       Returns : 
# x : Randomly selected X co-ordinate between 0 and R using
#     the probability distribution P(r)=2*pi*r*dr/(pi*R**2)
# y : Y co-ordinate set to 0
# z : Z co-ordinate calculated such the (x,y,z) touches sphere


def newphot():
	r=sqrt(random())*R
	x=random()*H/2
	y=r
	z=R-sqrt(R**2-r**2)
	mu=0.5
	phi=2*pi*random()
	return x,y,z,mu,phi

def makegrid_cart():
	opac=np.ones([resx,resy,resz])
	return opac

opac=makegrid_cart()    # 3-D grid of optical density

# finds the nearest wall

def findwall_cart(x,y,z,ix,iy,iz,nx,ny,nz):
	s=sx=sy=sz=0.0
	ind=xind=yind=zind=0
	which_wall=0

	if (nx>0):
		sx=(xface[ix+1]-x)/nx
		xind=1
	elif (nx<0) :
		sx=(xface[ix]-x)/nx
		xind=-1
	else:
		sx=1000*xmax

	if (ny>0):
		sy=(yface[iy+1]-y)/ny
		yind=1
	elif(ny<0) :
		sy=(yface[iy]-y)/ny
		yind=-1
	else:
		sy=1000*ymax

	if (nz>0):
		sz=(zface[iz+1]-z)/nz
		zind=1
	elif (nz<0):
		sz=(zface[iz]-z)/nz
		zind=-1
	else:
		sz=1000*zmax

	s=min(sx,sy,sz)
	if(s==sx):
		ind=xind
		which_wall=1*ind    #1,2,3 -> x,y,z , +-
	elif(s==sy):
		ind=yind
		which_wall=2*ind
	elif(s==sz):
		ind=zind
		which_wall=3*ind
	return (s,which_wall)

# takes the photon to the next scattering or exiting location

def tauint(x,y,z,ix,iy,iz,nx,ny,nz,tau,exitflag):
	taurun=0.0
	taucell=0.0
	while (True):
		s,which_wall=findwall_cart(x,y,z,ix,iy,iz,nx,ny,nz)
		taucell=s*opac[ix,iy,iz]
		if ((taurun+taucell)>tau):
			break

		taurun=taurun+taucell
		x=x+nx*s
		y=y+ny*s
		z=z+nz*s

		if abs(which_wall)==1:
			ix+=which_wall
		elif abs(which_wall)==2:
			iy+=(which_wall/2)
		elif abs(which_wall)==3:
			iz+=(which_wall/3)
		if(ix<0) or (iy<0) or (iz<0) or (ix==resx-1) or (iy==resy-1) or (iz==resz-1):
			exitflag=1
			x=x+nx*s
			y=y+ny*s
			z=z+nz*s
			break
	if (exitflag==0):
		s=(tau-taurun)/opac[ix,iy,iz]
		x=x+nx*s
		y=y+ny*s
		z=z+nz*s
	return (x,y,z,ix,iy,iz,exitflag)


# transfer tracks the photons through the sphere and returns
#the final polarizations binned by exit angle(mu).
#       Parameters : None
#       Returns :
# poli,polq,polu,polv : The exit values of Stokes' parameters
#                       binned by exit value of mu
# energy : The number of photons in each mu bin


def transfer():

	i=1
	energy=np.zeros(mubins)
	poli=np.zeros(mubins)
	polq=np.zeros(mubins)
	polu=np.zeros(mubins)
	polv=np.zeros(mubins)
	while(i<=nphoton):
		exitflag=0
		x,y,z,mu,phi=newphot()
		nx=sqrt(1-mu**2)*cos(phi)
		ny=sqrt(1-mu**2)*sin(phi)
		nz=mu
		z=z-R   #Origin shifted to centre of circle at (0,0,R)
		ix=locate(xface,x)
		iy=locate(yface,y)
		iz=locate(zface,z)
		aflag=0
		s=np.array([1,0,0,0]) #Unpolarized
		while(True):
			tau=-log(random())
			x,y,z,ix,iy,iz,exitflag=tauint(x,y,z,ix,iy,iz,nx,ny,nz,tau,exitflag)

			#if (sqrt(y**2+z**2)>=R) or (abs(x)>H/2):
			#	exitflag=1
			#	break
			if (exitflag==1):
				break
			if (random()<albedo):
				mu,phi,s=raylscatt(mu,phi,s)
				#mu,phi=isoscatt()
			else:
				aflag=1
				break
		if (exitflag==1):
			i=i+1
			for l in range(mubins-1):
				if (muspace[l]<=mu) and (mu<muspace[l+1]):
					poli[l]=poli[l]+s[0]
					polq[l]=polq[l]+s[1]
					polu[l]=polu[l]+s[2]
					polv[l]=polv[l]+s[3]
					energy[l]=energy[l]+1
			if (i%1000==0):
				print i,"\tphotons completed."
	return (poli[:mubins-1],polq[:mubins-1],polu[:mubins-1],polv[:mubins-1],energy[:mubins-1])
 
# output calculates the degree of polarization and
# angle of polarization as a function of exit angle,
# saves it in a file and plots it.
#       Parameters : None
#       Returns : None


def output():

	poli,polq,polu,polv,energy=transfer()
	poli=poli/energy
	polq=polq/energy
	polu=polu/energy
	polv=polv/energy
	mu_plot=muspace[:mubins-1]+(1./mubins) #mu values at centre of bins
	theta=180*np.arccos(mu_plot)/np.pi
	deg=np.sqrt(polq**2+polu**2)/poli
	#np.savez("blobtau10",theta,deg)
	plt.plot(theta,deg)
	plt.axvline(x=90.)
	plt.legend()
	plt.show()

def main():
	output()

if (__name__=='__main__'):
	main()
