import numpy as np
import matplotlib.pyplot as plt
import os.path
from os.path import join as pjoin

file_name = 'part.txt'
path = './'
if os.path.exists(pjoin(path,file_name)):
	x, y = np.loadtxt(pjoin(path,file_name),unpack=True)
else:
	print('No Data')
	exit()

fig = plt.figure(figsize=(6,5))
left, bottom, width, height = 0.1, 0.1, 0.8, 0.8
ax = fig.add_axes([left, bottom, width, height])
plt.scatter(x,y)
plt.xlim([0, 6.5])
plt.ylim([0, 6.5])
plt.show()
