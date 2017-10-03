################################################################################
################################################################################
####                                                                        ####
####  Tool to visualize the trajectory of a single particles                ####
####                                                                        ####
####  WARNING: only works with files in H5hut format generated after        #### 
####  This routine assume that the particles are always in the same order   ####
####  27th of Septemebr of 2017                                             ####
####                                                                        ####
####  Diego Gonzalez                                                        ####
####  September 2017                                                        ####
####                                                                        ####
################################################################################
################################################################################
import numpy as np
import h5py
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

################################################################################
# Input data
basename = '/home/diego/Codes/ECsim/ecsim/Test/data/TwoCoils-Partcl'
cStart = 400    # First cycle
cEnd   = 0     # Last cycle
cStep  = -10    # Step between files
spec   = 0     # Index of the specie to plot
tagR   = 0     # tagR: rank where it was born)
tagP   = 4     # tagP: particle tag

################################################################################

hisX = []
hisY = []
hisZ = []

ind = 0
for i in range(cStart,cEnd,cStep):
	filename = "%s_%06d.h5" % (basename, i)
	keytagR = "Step#0/tagR_%d" % spec
	keytagP = "Step#0/tagP_%d" % spec
	keyX = "Step#0/x_%d" % spec
	keyY = "Step#0/y_%d" % spec
	keyZ = "Step#0/z_%d" % spec
	print filename

	file  = h5py.File(filename, driver='core')
	X = file.get(keyX)
	Y = file.get(keyY)
	Z = file.get(keyZ)
	aTagR = file.get(keytagR)[()]
	aTagP = file.get(keytagP)[()]
	found = False
	for p in range(X.shape[0]):
		if (aTagR[p] == tagR and aTagP[p] == tagP):
			hisX.append(X[p])
			hisY.append(Y[p])
			hisZ.append(Z[p])
			found = True
			break
 
	if (found == False):
		break 

fig = plt.figure()
ax = fig.add_subplot(111,projection='3d')
ax.plot(hisX,hisY,hisZ,'b-')
plt.plot(hisX, hisY, hisZ, 'b.')
plt.show()

