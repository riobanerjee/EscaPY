# Parallel rays through uniform blob from
# under it.

from math import *
from random import *
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


seed(54366)
albedo=1.0	
nphoton=10000000
taumax=1.
mubins=20 # Number of bins = mubins-1. mubins is 
			#the number of edges of bins
R=1.	  # Radius of Sphere
a=0.75	  # Required in the polarisation matrices
muspace=np.linspace(-1,1,mubins)

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
# and travelling parallel to the +Z axis.
#       Parameters : None
#       Returns : 
# x : Randomly selected X co-ordinate between 0 and R using
#     the probability distribution P(r)=2*pi*r*dr/(pi*R**2)
# y : Y co-ordinate set to 0
# z : Z co-ordinate calculated such the (x,y,z) touches sphere


def newphot():
	r=sqrt(random())*R
	x=r
	y=0.
	z=R-sqrt(R**2-r**2)
	mu=1.0
	phi=2*pi*random()
	return x,y,z,mu,phi


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
		z=z-R   #Origin shifted to centre of circle at (0,0,R)
		aflag=0
		s=np.array([1,0,0,0]) #Unpolarized
		while(True):
			tau=-log(random())
			L=tau/taumax
			x=x+L*sqrt(1-mu*mu)*cos(phi)
			y=y+L*sqrt(1-mu*mu)*sin(phi)
			z=z+L*mu
			if (sqrt(x**2+y**2+z**2)>=R):
				exitflag=1
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
	np.savez("blobtau10",theta,deg)
	plt.plot(theta,deg)
	plt.axvline(x=90.)
	plt.xlabel("Angle")
	plt.ylabel("Degree of Polarisation")
	plt.legend()
	plt.show()

if (__name__=='__main__'):
	output()
