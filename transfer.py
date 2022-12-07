import numpy as np
import random
import matplotlib

npart = 7 #7 points representing the particles are along each direction
N = (npart - 1)*(npart - 1) #N is the number of real particles
/*
	Explanation: array of 36*36,
	Since over a length space the first and the last points are the same, out of 7 only 6 will be real particles
	hence the array is of size 36 by 36.
*/

command 'rm -rf output/*'
command 'mkdir output'

temp = 0.1 #temperature in scaled units
box = 5.5 #total box length
a = box/(npart - 2) #how much space is available to one molecule

	#Explanation: Suppose 3 molecules occupy a length of 2.
	#Then how much length is available for each molecule?
	#This is given by 2/3=0.66. This is in fact the length of freedom for the molecule.
	#The reason of taking npart-2 (actual particles are npart-1) is an attempt to get how much length a particle does fully enjoy.

box = box + a
boxby2 = box*0.5 #half of the updated full box length
rc2 = 2.5 #cut off radius and is manually decided

if rc2 > boxby2:
	print("Cut-off radius greater than half-box length")
	exit()

rc2 = rc2*rc2 #find rcut square
ecut = 1./(rc2*rc2*rc2)
ecut = 4.*ecut*(ecut-1) # value of the LJ potential at the rcut.
delt = 1.e-3; # time step
timesteps = (int) (1./delt) # no of time steps required to reach 1 unit in steps of delt.
timesteps = timesteps + 10002 # + two more time steps
#this is baseless anyway and is done just to make the total timesteps to be equal to 1000.

#Initialise the MD program
#Initialize the positions to start the program
#i and j are the particle indices along x and y direction respectively

for i in range(1, npart):
	for j in range(1, npart):
		I = j + (i - 1)*(npart - 1)
		vx[I] = random() - 0.5
		vy[I] = random() - 0.5
#Since, initial kinetic energy must be zero, and energy is conserved
#we initiate with zero sum of velocities along x and y direction
#sumvx and sumvy are the component velocity of the centre of mass

sumvx = 0.0
sumvy = 0.0
#square sum of velocities
sumv2x = 0.0
sumv2y = 0.0
sumv2 = 0.0

for i in range(1, N+1):
	sumvx = sumvx + vx[i]
	sumvy = sumvy + vy[i]
	sumv2x = sumv2x + vx[i]*vx[i]
	sumv2y = sumv2y + vy[i]*vy[i]

#Average sum of velocities
sumvx = sumvx/N
sumvy = sumvy/N
#Average sum of square velocities
sumv2x = sumv2x/N
sumv2y = sumv2y/N
#Sum of square velocities
sumv2 = sumv2x+ sumv2y
#velocity scaling factor (kT/(0.5mv^2)
fs = sqrt(2.*temp/sumv2)

#Subtract the individual velocities from the average sum of
#velocities and scale it by the scaling factor

for i in range(1, npart):
	for j in range(1, npart):
		I = j + (i-1)*(npart-1) #particle index
		vx[I] = (vx[I] - sumvx)*fs
		vy[I] = (vy[I] - sumvy)*fs
		#Find the previous positions
		xm[I] = x[I] - vx[I]*delt
		ym[I] = y[I] - vy[I]*delt

#Start the time loop
t = 0.0




















