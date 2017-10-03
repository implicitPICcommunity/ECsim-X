################################################################################
################################################################################
####                                                                        ####
####  Get the key of the particles inside a given region.                   ####
####                                                                        ####
####  WARNING: only works with files in H5hut format generated after        #### 
####  This routine assume that the particles are always in the same order   ####
####  27th of Septemebr of 2017                                             ####
####                                                                        ####
####  Diego Gonzalez                                                        ####
####  October 2017                                                          ####
####                                                                        ####
################################################################################
################################################################################
import numpy as np
import h5py
import sys
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import csv

################################################################################
# Input data
filename = '/home/diego/Codes/ECsim/Python/particles/TrackParticles/Data/part70000/TwoCoils2D-Partcl_060001.h5'
spec   = 0     # Index of the specie to plot
cyl = True
x0 = 3.4
x1 = 3.7
y0 = 24.8
y1 = 25.2

################################################################################

hisX = []
hisY = []
hisZ = []
hisW = []

keytagR = "Step#0/tagR_%d" % spec
keytagP = "Step#0/tagP_%d" % spec
keyX = "Step#0/x_%d" % spec
keyY = "Step#0/y_%d" % spec
keyZ = "Step#0/z_%d" % spec
print filename

file  = h5py.File(filename, driver='core')
X = file.get(keyX)[()]
Y = file.get(keyY)[()]
Z = file.get(keyZ)[()]
aTagR = file.get(keytagR)
aTagP = file.get(keytagP)
n = 0
fd = open("particlesID.csv","wt")
writer = csv.writer(fd)
for p in range(X.shape[0]):
	if (X[p] > x0 and X[p] < x1 and Y[p] > y0 and Y[p] < y1):
		writer.writerow((aTagR[p], aTagP[p]))
		n = n +1
fd.close()
print n

