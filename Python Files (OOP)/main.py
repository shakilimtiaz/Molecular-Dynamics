import numpy as np
import func
import os

# delete and create the folder output
os.system('rm -rf output')
os.system('mkdir output')


# define x and y as the same variable as in func.py
# Initialize the arrays
x = np.zeros(func.N)
y = np.zeros(func.N)
xm = np.zeros(func.N)
ym = np.zeros(func.N)
vx = np.zeros(func.N)
vy = np.zeros(func.N)
fx = np.zeros(func.N)
fy = np.zeros(func.N)

# initialize x and y: this modifies x and y variables 
# already defined to be array of zeros 
x, y, xm, ym, vx, vy = func.init(x, y, vx, vy)
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# open a file to write the initialized values of particles 
fp = open('part.txt','w+')
ff = open('forces.txt','w+')
# write the new values of x and y to the file 'part.txt' 
for i in range(len(x)):
    fp.write('%f \t %f \n' %(x[i], y[i]))
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# define the time
t = 0.0 

for k in range(func.timesteps):    # func.timesteps    
    
    # Calculate the forces
    fx, fy, en = func.forces(x, y, fx, fy, func.en)
    
    # write the forces on each particle to a file
    #for i in range(len(func.fx)):
    #   ff.write('%f \t %f \n'%(func.fx[i],func.fy[i]))

    # integrate on the equations of motion
    # push the particles to a new position 
    # assume the centre of mass velocities to be zero    
    x, y, xm, ym, sumv2 = func.integrate(x, y, xm, ym, vx, vy, fx, fy)   

    # calculate the instantaneous temperature
    temp = sumv2/(2*func.N)
    # calculate total energy per particle
    etot = (en + 0.5*sumv2)/func.N
    if (k%100)==0:
        path = './output/'
        file = open(str(path)+'part_%d.txt'%(k),'w+')
        
        for i in range(1, func.npart):
            for j in range(1, func.npart):
                I = j + (i-1)*(func.npart-1)
                file.write('%f \t %f \n' %(x[I-1],y[I-1]))        
              
    # forward the time
    t = t + func.delt

#fp.close()
#file.close()



    
    

