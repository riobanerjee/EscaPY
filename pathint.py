# Program to calculate the final intensity when looking
# down at an angle of (Considering only source function)
# Source function is jnu.


import numpy as np
import matplotlib.pyplot as plt
from math import *


# jnu follows the following function of co-ordinates.
# 		Parameters:
# x,y,z - Cartesian Co-ordinates


def rule(x,y,z):
	return x*y*z


# Makes the 3-D grid of jnu by following rule(x,y,z)
# 		Parameters:
# f - rule
# resx,resy,resz - Resolution in x, y and z co-ords
#		Returns:
# jnu - Source function


def makegrid(f,resx,resy,resz):
	jnu=np.zeros([resx,resy,resz])
	for i in range(resx):
		for j in range(resy):
			for k in range(resz):
				jnu[i,j,k]=f(i,j,k)
	return jnu


# Performs the path integral at the viewing angle
# 		Parameters:
# inuinit - Initial Intensity
# jnu - Source function
# resx,resy,resz - Resolution in x, y and z co-ords
# va - Viewing angle
# dz - z-length(height) of 1 cell
#		Returns:
# inu - Final intensity


def path(inuinit,jnu,resx,resy,resz,va,dz):

	inu=inuinit
	for i in range(resx):
		for j in range(resy):
			for k in range(resz):
				inu[i,j]=jnu[i,j,k]*dz/cos(va)
	return inu


# Initializes parameters, calls all functions in order
# and plots the results.


def main():
	resx,resy,resz = 100 , 100 , 100
	xmin,ymin,zmin = 0. , 0. , 0.
	xmax,ymax,zmax = 10. , 10. , 10.

	inu=np.zeros([resx,resy])
	jnu=makegrid(rule,resx,resy,resz)
	dx,dy,dz = (xmax-xmin)/resx , (ymax-ymin)/resy , (zmax-zmin)/resz
	va=input("Enter viewing angle in radians ")   # in radians, from z axis on the xz plane
	inu=path(inu,jnu,resx,resy,resz,va,dz)
	plt.imshow(inu,'gray')
	plt.xlabel('X')
	plt.ylabel('Y')
	plt.colorbar()
	plt.title("Viewing angle = {} radians\njnu(x,y,z) = xyz".format(va))
	plt.show()


if __name__=='__main__':
	main()