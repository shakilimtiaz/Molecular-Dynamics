import numpy as np
import random
import matplotlib
import os
import math

npart = 7 #7 points representing the particles are along each direction
N = (npart - 1)*(npart - 1) #N is the number of real particles

#Explanation: array of 36*36,
#Since over a length space the first and the last points are the same, out of 7 only 6 will be real particles
#hence the array is of size 36 by 36.

x = [0] * N #x array
y = [0] * N #y array
xm = [0] * N #approximation of the previous x particle position
ym = [0] * N #approximation of the previous y particle position
vx = [0] * N #x velocity array
vy = [0] * N
fx = [0] * N
fy = [0] * N

os.system('rm -rf Output')
os.system('mkdir Output')

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

rand_int = np.random.randint(0,1)


for i in range(1, npart):
	for j in range(1, npart):
		I = j + (i - 1)*(npart - 1)
		if i == 1 and j == 1:
			x[0] = a/2
			y[0] = a/2
		else:
			x[I-1] = x[0] + (i-1)*a + rand_int*0.0005
			y[I-1] = y[0] + (i-1)*a + rand_int*0.0005

for i in range(1, npart):
	for j in range(1, npart):
		I = j + (i - 1)*(npart - 1)
		vx[I-1] = float(rand_int - 0.5)
		vy[I-1] = float(rand_int - 0.5)
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
	sumvx = sumvx + vx[i-1]
	sumvy = sumvy + vy[i-1]
	sumv2x = sumv2x + vx[i-1]*vx[i-1]
	sumv2y = sumv2y + vy[i-1]*vy[i-1]

#Average sum of velocities
sumvx = sumvx/N
sumvy = sumvy/N
#Average sum of square velocities
sumv2x = sumv2x/N
sumv2y = sumv2y/N
#Sum of square velocities
sumv2 = sumv2x+ sumv2y
#velocity scaling factor (kT/(0.5mv^2)
fs = math.sqrt(2.*temp/sumv2)

#Subtract the individual velocities from the average sum of
#velocities and scale it by the scaling factor

for i in range(1, npart):
	for j in range(1, npart):
		I = j + (i-1)*(npart-1) #particle index
		vx[I-1] = (vx[I-1] - sumvx)*fs
		vy[I-1] = (vy[I-1] - sumvy)*fs
		#Find the previous positions
		xm[I-1] = x[I-1] - vx[I-1]*delt
		ym[I-1] = y[I-1] - vy[I-1]*delt

#Start the time loop
t = 0.0

for K in range(0, timesteps):
	for i in range(1, npart):
		for j in range(1, npart):
			I = j + (i - 1)*(npart - 1)
			fx[I-1] = 0.0
			fy[I-1] = 0.0

	#start with zero potential energy
	en = 0.0

#loop over all the particles:
#i should start from 1 to N-1 and j from 2 to N
#just to distinguish between the particles

	for i in range(1, N):
		for j in range(1, N+1):
				#calculate the difference between the molecules one by one*/
				xr = x[i-1] - x[j-1]
				yr = y[i-1] - y[j-1]
				#If the difference happens to be out of the box,
				#create an image particle in the box maintaining the periodicity

				if xr > boxby2:
					xr = xr-box
				elif xr < -boxby2:
					xr = xr + box

				if yr > boxby2:
					yr = yr-box
				elif yr < -boxby2:
					yr = yr + box
				
				#square the difference
				r2 = xr*xr + yr*yr

				#If the difference is more than a certain number and less than
				#the cutoff radius then only calulate the forces

				if r2 > 1.e-12 and r2 < rc2:
					r2i = 1./r2
					r6i = r2i*r2i*r2i
					ff = 48.*r2i*r6i*(r6i-0.5) #potential energy

					#By Newton's thrid law if the first particle is
					#moved to right the seond in the interaction would move left

					fx[i-1] = fx[i-1] + ff*xr
					fy[i-1] = fy[i-1] + ff*yr
					fx[j-1] = fx[j-1] - ff*xr
					fy[j-1] = fy[j-1] - ff*yr
					en = en + 4.*r6i*(r6i-1.) - ecut

#Integrate the equations of motion: push the particles to a new position
#Once again initiate the centre of mass velocities to be zero
		sumvx = 0.0
		sumvy = 0.0
		sumv2x = 0.0
		sumv2y = 0.0

		for i in range(1, N + 1):
			#Calculate the new positions from the verlet algorithm
			xx = 2*x[i-1] - xm[i-1] +delt*delt*fx[i-1]
			yy = 2*y[i-1] - ym[i-1] +delt*delt*fy[i-1]

			if xx > box:
				xx = xx - box
			elif xx < 0:
				xx = xx + box
			
			if yy > box:
				yy = yy - box
			elif yy < 0:
				yy = yy + box
			
			xr = xx - xm[i-1]
			yr = yy - ym[i-1]

			if xr > boxby2:
				xr = xr-box;
			elif xr < -boxby2:
				xr = xr + box;

			if yr > boxby2:
				yr = yr-box;
			elif yr < -boxby2:
				yr = yr + box
			
			vx[i-1] = xr/(2.*delt);
			vy[i-1] = yr/(2.*delt);

			#Update the centre of mass velocities
			sumvx = sumvx + vx[i-1]
			sumvy = sumvy + vy[i-1]
			sumv2x = sumv2x + vx[i-1]*vx[i-1]
			sumv2y = sumv2y + vy[i-1]*vy[i-1]

			#Update positions in previous time
			xm[i-1] = x[i-1]
			ym[i-1] = y[i-1]
			
			#Update positions in current time
			x[i-1] = xx
			y[i-1] = yy

		sumv2 = sumv2x + sumv2y

		#Instantaneous temperature
		temp = sumv2/(2.*N)

		#Total energy per particle
		etot = (en + 0.5*sumv2)/N

#Forward the time
t = t + delt
