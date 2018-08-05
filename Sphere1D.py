# Tracks single photon through Spherical atmosphere

from math import *
from random import *
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

seed(10480)
albedo=1.0											#initial values
nlevel=11
mubins=20
nphoton=1
taumax=10.
rmin=1.						
rmax=10.
rstar=1.
limb=0.
rhoexp=-1.
ngrid=10
imagenx=99
r=np.linspace(rmin,rmax,ngrid)
									

def newphot():							#emits new photon
	mu=sqrt(random())
	phi=2*pi*random()
	xt=yt=zt=0.
	return (mu,phi,xt,yt,zt)

def newphotsph():
	mu,phi,xt,yt,zt=newphot()
	cosb=1-2*random()
	sinb=sqrt(1-cosb**2)
	if limb!=0:						
		mu=darkening(random(),mu)
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
			if(cospn>1):
				cospn=1.0
			elif(cospn<-1):
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
		else:
			costi=mu
	phi=phi+lp
	if(phi>2*pi):
		phi=phi-2*pi
	zt=cosb*rmin*1.0001
	xt=sinb*cos(lp)*rmin*1.0001
	yt=sinb*sin(lp)*rmin*1.0001
	return(mu,phi,xt,yt,zt)

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
    return opac
    
opac=makegrid()

def isoscatt():					#isotropically scatter photon
	mu=2*random()-1
	phi=2*pi*random()
	return (mu,phi)


def findWall(r1,r2,nx,ny,mu,xp,yp,zp,exitflag):
	
	
	dt=0
	root=[-999,-999,-999,-999]
	ind=[0,0,1,1]
	#if (sqrt(xp*xp+yp*yp+zp*zp)>r2) and (r2!=rmax) and (r1!=rmin):
	#	return(sqrt(xp*xp+yp*yp+zp*zp)-r2+0.001,-1,0)
	bb=(xp*nx+yp*ny+zp*mu)						
	cc=xp*xp+yp*yp+zp*zp-r1*r1
	dd1=bb*bb-cc
	ifound=True
	if (dd1>=0):
		dd1=sqrt(dd1)
		root[0]=(-bb+dd1)
		root[1]=(-bb-dd1)
	cc=xp*xp+yp*yp+zp*zp-r2*r2
	dd2=bb*bb-cc
	if (dd2>=0):
		dd2=sqrt(dd2)
		root[2]=(-bb+dd2)
		root[3]=(-bb-dd2)
	npos=0
	posroot=[0.,0.,0.,0.]
	for i in range(4):
		if(root[i]>0):
			posroot[npos]=root[i]
			npos+=1
	print r1,r2,sqrt(cc+r2*r2)
	if (npos==0):
		ifound=False
		t=0.00001
		print "But Why??",r1,r2,sqrt(cc+r2*r2),dd1,dd2
    		ioffset=1
	
	else:
		ifound=True
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
	
	if(exitflag!=1):				#-1 for r1, 1 for r2
		dt=t
		
	return(dt,ioffset,exitflag)

def tau_1d(tau,xold,yold,zold,nx,ny,mu,phi,ii,exitflag,jp,hp,kp):
	
	tsum=0.0
	xp=xold
	yp=yold
	zp=zold
	while(True):
		if(r[ii]==rmax):
			exitflag=1
			print "Exit"
			return(exitflag,xp,yp,zp,ii,jp,hp,kp)
			
		#print ii,r[ii],r[ii+1],sqrt(xp*xp+yp*yp+zp*zp)
		rtot=sqrt(xp*xp+yp*yp+zp*zp)
		ii=locate(r,rtot)
		ds,ioffset,exitflag=findWall(r[ii],r[ii+1],nx,ny,mu,xp,yp,zp,exitflag)
		if(exitflag==1):
			return(exitflag,xp,yp,zp,ii,jp,hp,kp)
		if(ds<1e-5):
			ds=1e-5
		opaci=opac[ii]
		dtau=opaci*ds
		if(ii!=0) or (ioffset!=-1):
			ii=ii+ioffset    
		
		if(tsum+dtau)>tau:
			ds=(tau-tsum)/opaci
			xp=xp+ds*nx
			yp=yp+ds*ny
			zp=zp+ds*mu
			return(exitflag,xp,yp,zp,ii,jp,hp,kp)
		
	
   		xp=xp+ds*nx*1.001           
   		yp=yp+ds*ny*1.001
   		zp=zp+ds*mu*1.001
   		tsum=tsum+dtau

		

def transfer(taumax=10,zmax=1):
	i=1
	x=[]
	y=[]
	z=[]
	energy=0
	fractx=2.*rmax/float(imagenx)
	jp=np.zeros(ngrid)
	hp=np.zeros(ngrid)
	kp=np.zeros(ngrid)
	image=np.zeros([imagenx,imagenx])
	while(i<=nphoton):

		mu,phi,xt,yt,zt=newphotsph()
		aflag=0			#aflag=1 means photon is absorbed
		exitflag=0
		ii=0
		jp[ii]=jp[ii]+1./mu					
   		hp[ii]=hp[ii]+1.
       		kp[ii]=kp[ii]+mu
		xold=yold=zold=0.
		while (exitflag==0) :
		
			nx=sqrt(1-mu**2)*cos(phi)
			ny=sqrt(1-mu**2)*sin(phi)
			aflag=0
			xold=xt
			yold=yt
			zold=zt
			x.append(xt)
			y.append(yt)
			z.append(zt)
			tau=-log(random())
			exitflag,xt,yt,zt,ii,jp,hp,kp=tau_1d(tau,xold,yold,zold,nx,ny,mu,phi,ii,exitflag,jp,hp,kp)
			if (exitflag==1) :
				i+=1						
			if (random()<albedo):
				mu,phi=isoscatt()		#photon scattered when not absorbed	
			else:
				aflag=1
				i+=1
				break
		
		if (exitflag==1) and (aflag==0):				
			energy+=1
			
			
	return x,y,z,energy

def output():
	x,y,z,en=transfer()
	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')
	ax.plot(x, y, z)
	u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
	xs = rmin*np.cos(u)*np.sin(v)
	ys = rmin*np.sin(u)*np.sin(v)
	zs = rmin*np.cos(v)
	ax.plot_wireframe(xs, ys, zs, color="g")
	u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
	xs = rmax*np.cos(u)*np.sin(v)
	ys = rmax*np.sin(u)*np.sin(v)
	zs = rmax*np.cos(v)
	ax.plot_wireframe(xs, ys, zs, color="r")
	ax.scatter(x[0],y[0],z[0],color='g',s=100,label='Start')
	ax.scatter(x[-1],y[-1],z[-1],color='r',s=100,label='end')
	print (sqrt(x[-1]**2+y[-1]**2+z[-1]**2))
	plt.legend()
	plt.show()
	#energy= float(en)/nphoton
	#jp=jp/nphoton
	#hp=hp/nphoton
	#kp=kp/nphoton
	#image=image/nphoton
	#plt.imshow(image,'gray')
	#plt.colorbar()
	#plt.show()
	#print energy,jp,hp,kp
output()
