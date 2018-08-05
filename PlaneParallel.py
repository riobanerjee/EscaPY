# Rays isotropically emitting from a point at the base of a uniform
# plane parallel atmosphere. 

from math import *
from random import *
import numpy as np
import matplotlib.pyplot as plt

seed(87654)
albedo=1.											#initial values
nlevel=10
mubins=20
nphoton=100000											#Takes a lot of time


def newphot():											#emits new photon
	mu=sqrt(random())
	phi=2*pi*random()
	xt=np.array([0.,0.,0.])
	return (mu,phi,xt)

def isoscatt():											#isotropically scatter photon
	mu=2*random()-1
	phi=2*pi*random()
	return (mu,phi)

def moments(x1,x2,cost,nlevel,jp,hp,kp,jm,hm,km,zmax):						#calculates the intensity moments
	if((x1[2]>0) and (x2[2]>0) and (int(x1[2]*nlevel/zmax)==int(x2[2]*nlevel/zmax))):	#when photon is in atmosphere but does not cross a level
		return (jp,hp,kp,jm,hm,km)
	if(cost>0):										#counting from bottom to top
		if(x1[2]<0):									#l1 is initial level and l2 is updated level
			l1=0									#level 0 when photon is below atmosphere
		else:
			l1=int(x1[2]*nlevel/zmax)
		if(x2[2]>1):
			l2=nlevel								#crossed last level(photon exits from top)
		else:
			l2=int(x2[2]*nlevel/zmax)
		for i in range(l1,l2):
			jp[i]=jp[i]+1./cost
			hp[i]=hp[i]+1.
			kp[i]=kp[i]+cost
	elif(cost<0):										#counting from top to bottom
		l1=int(x1[2]*nlevel/zmax)
		if(x2[2]<0):
			l2=0
		else:
			l2=int(x2[2]*nlevel/zmax)
		for i in range(l2,l1):
			jm[i]=jm[i]+1./abs(cost)
			hm[i]=hm[i]-1.
			km[i]=km[i]+abs(cost)
	return (jp,hp,kp,jm,hm,km)


def transfer(taumax=10,zmax=1):
	i=1
	energy=np.zeros(mubins)
	jp=np.zeros(nlevel)
	hp=np.zeros(nlevel)
	kp=np.zeros(nlevel)
	jm=np.zeros(nlevel)
	hm=np.zeros(nlevel)
	km=np.zeros(nlevel)
	while(i<=nphoton):
		
		mu,phi,xt=newphot()
		aflag=0											#aflag=1 means photon is absorbed
		
		
		while (xt[2]>=0) and (xt[2]<=zmax) :
			x1=xt
			nx=np.array([cos(phi)*sqrt(1-mu**2),sin(phi)*sqrt(1-mu**2),mu])			#direction of propagation
			tau=-log(random())
			L=tau*zmax/taumax								#length of step
			xt=xt+(L*nx)
			x2=xt
			(jp,hp,kp,jm,hm,km)=moments(x1,x2,mu,nlevel,jp,hp,kp,jm,hm,km,zmax)
			if (xt[2]<0) or (xt[2]>zmax) :							#when photon exits from bottom or top of atmosphere
				break	
			if (random()<albedo):
				mu,phi=isoscatt()							#photon scattered when not absorbed	
			else:
				aflag=1
				break
		
		if (aflag==0) and (xt[2]>zmax):								#counts photons exiting from the top of the atmosphere
			l=(int(mubins*mu))
			if (i%1000==0):
				print i,"\tphotons completed"
			i+=1										#and bins them in respective mu bins (last mu value)
			energy[l]=energy[l]+1
			
	return energy,jp,hp,kp,jm,hm,km

def output():
	tau=np.linspace(nlevel,1,nlevel)
	theta=np.zeros(mubins)
	dthet=1./mubins
	halfw=0.5*dthet
	mu=np.zeros(mubins)

	for j in range(mubins):										#mu,theta values at the middle of the mu bins
		theta[j]=acos(j*dthet+halfw)
		mu[j]=j*dthet+halfw
	
	energy,jp,hp,kp,jm,hm,km=transfer()
	j=(jp+jm)/nphoton									#moments
	h=(hp+hm)/nphoton
	k=(kp+km)/nphoton
	f=k/j											#Eddington factors
	g=h/j
	intensity=energy*mubins/(2*nphoton*mu)							#normalised intensity
	energy=energy/nphoton
	theta=theta*180/np.pi

	inte=plt.figure()
	ax1=inte.add_subplot(111)									#plots
	ax1.scatter(theta,intensity)
	ax1.set_xlabel('theta')
	ax1.set_ylabel('Normalised Intensity')
	mom=plt.figure()
	ax2=mom.add_subplot(111)
	ax2.plot(tau,j,label='J')
	ax2.plot(tau,h,label='H')
	ax2.plot(tau,k,label='K')
	ax2.set_xlabel('tau')
	ax2.invert_xaxis()
	ax2.set_ylabel('Intensity Moments')
	eddie=plt.figure()
	ax3=eddie.add_subplot(111)
	ax3.plot(tau,f,label='f')
	ax3.plot(tau,g,label='g')
	ax3.set_xlabel('tau')
	ax3.invert_xaxis()
	ax3.set_ylabel('Eddington Factor')
	ax3.legend()
	plt.show()

output()
