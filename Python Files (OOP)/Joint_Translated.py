import numpy as np
import random
import os
import math

############################################################

npart = 7
N = (npart-1)*(npart-1)

# delete and create the folder output
os.system('rm -rf output')
os.system('mkdir output')

# define the box length 
box = 5.5
# calculate how much space is available to one molecule/atom
a = box/(npart - 2) 
# redefine box length
box = box + a 
# define half the boxlength 
boxby2 = 0.5*box 
# define the cut-off radius 
rc2 = 2.5 
# put a caution
if rc2>boxby2:
    print("cut-off radius is greater than half-box length:modify func.py")
    exit()
# define square of the rcut
rc2 = rc2*rc2 
# define ecut (cut-off energy)
ecut = 1./(rc2*rc2*rc2) 
# find the value of the LJ potential at rcut
ecut = 4.*ecut*(ecut-1)

# define time step
delt = 1.e-3
timesteps = (int)(1./delt) + 10000

# define temperature in scaled units
temp = 0.1  

# define the potential energy
en = 0.0 

rand_int = np.random.randint(0,1)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# define x and y as the same variable as in func.py
# Initialize the arrays

x = np.zeros(N)
y = np.zeros(N)
xm = np.zeros(N)
ym = np.zeros(N)
vx = np.zeros(N)
vy = np.zeros(N)
fx = np.zeros(N)
fy = np.zeros(N)

# initialize x and y: this modifies x and y variables 
# already defined to be array of zeros 

for i in range(1, npart):
	for j in range(1, npart):
			I = j + (i-1)*(npart-1)
			if i==1 and j==1:
				 x[0] = a/2
				 y[0] = a/2 
			else:
				 x[I-1] = x[0] + (i-1)*a + rand_int*0.0005
				 y[I-1] = y[0] + (j-1)*a + rand_int*0.0005
    
for i in range(1, npart):
	for j in range(1, npart):
			I = j + (i-1)*(npart-1)
			vx[I-1] = float(rand_int - 0.5)
			vy[I-1] = float(rand_int - 0.5)

# average sum of velocities
sumvx = np.average(vx)
sumvy = np.average(vy)
# average sum of square velocities
sumv2x = np.average(vx*vx)
sumv2y = np.average(vy*vy)
# sum of average square velocities
sumv2 = sumv2x + sumv2y
# velocity scaling factor (kT/(0.5mv**2)) 
fs = math.sqrt(2.*temp/sumv2)

# subtract the average sum of velocities from 
# individual velocities and scale it by the scaling factor
vx[:] = (vx[:] - sumvx)*fs
vy[:] = (vy[:] - sumvy)*fs
# find the previous positions
xm = x - vx*delt
ym = y - vy*delt     
 

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# open a file to write the initialized values of particles 

fp = open('part.txt','w+')
fm = open('forces.txt','w+')

# write the new values of x and y to the file 'part.txt' 

for i in range(len(x)):
    fp.write('%f \t %f \n' %(x[i], y[i]))

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# define the time

t = 0.0 

for k in range(timesteps):	
	for i in range(1, npart):
		for j in range(1, npart):
			I = j + (i - 1)*(npart - 1)
			fx[I-1] = 0.0
			fy[I-1] = 0.0

	#start with zero potential energy
	en = 0.0
	# -------------------------------------------------------------------------
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


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    # write the forces on each particle to a file
	for i in range(len(fx)):
	    fm.write('%f \t %f \n' %(fx[i], fy[i]))

    # integrate on the equations of motion
    # push the particles to a new position 
    # assume the centre of mass velocities to be zero    

	for i in range(1, N+1):
		# calculate the new positions from the velocity verlet algorithm
		xx = 2*x[i-1] - xm[i-1] + delt*delt*fx[i-1]
		yy = 2*y[i-1] - ym[i-1] + delt*delt*fy[i-1]

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
			 xr = xr - box
		elif xr < -boxby2:
			 xr = xr + box

		if yr > boxby2:
			 yr = yr-box
		elif yr < -boxby2:
			 yr = yr + box

		vx[i-1] = xr/(2.*delt)
		vy[i-1] = yr/(2.*delt)

		#Update the centre of mass velocities
		sumvx  = np.sum(vx)
		sumvy  = np.sum(vy)
		sumv2x = np.sum(vx*vx)
		sumv2y = np.sum(vy*vy)      

		#Update positions in previous time
		xm[i-1] = x[i-1]
		ym[i-1] = y[i-1]

		#Update positions in current time
		x[i-1] = xx
		y[i-1] = yy  

		sumv2 = sumv2x + sumv2y


# calculate the instantaneous temperature
	temp = sumv2/(2*N)
# calculate total energy per particle
	etot = (en + 0.5*sumv2)/N

	# Output the results into a file named 'output' at every 100 steps
	if(k	%100)==0:
		path = './output/'
		file =  open(str(path)+"part_%d.txt" %(k),"w+")
		for i in range(npart):
			for j in range(npart):				
				I = j + (i-1)*(npart-1)
				file.write("%f \t %f\n" %(x[I-1],y[I-1]))
		  
# forward the time
t = t + delt

fp.close()
#file.close()