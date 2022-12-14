
import numpy as np
import matplotlib.pyplot as plt
import os.path
from os.path import join as pjoin

k = 0
fig = plt.figure(figsize=(6,5))
left, bottom, width, height = 0.1, 0.1, 0.8, 0.8
ax = fig.add_axes([left, bottom, width, height])

while k<=900:
	file_name = 'part_%d.txt'%(k)
	path = './output/'

	if os.path.exists(pjoin(path,file_name)):
		x, y = np.loadtxt(pjoin(path,file_name),unpack=True)
	else:
		print('No Data')
		exit()

	plt.cla() # Remove this command and you see a cluster of all previous timesteps	
	plt.scatter(x,y)
	plt.xlim([0, 6.5])
	plt.ylim([0, 6.5])
	plt.pause(0.05)
	k = k + 100
plt.show()
