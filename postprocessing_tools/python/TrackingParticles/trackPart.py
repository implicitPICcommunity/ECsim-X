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
import csv

################################################################################
# Input data
################################################################################

# File with particles
basename = '/home/diego/Codes/ECsim/Python/particles/TrackParticles/Data/part70000/TwoCoils2D-Partcl'
# File with the IDs of the particles to track
fileids = 'particlesID.csv'
cStart = 60010    # First cycle
cEnd   = 62600    # Last cycle
cStep  = 100      # Step between files
spec   = 0        # Index of the specie to plot
cyl = True        # True if the particles files are in cylindrical coordinates
dt = 0.5          # Time step used in the simulation
Lx = 50           # Size of the box in 3D cartesian
Ly = 50           # Size of the box in 3D cartesian
Lz = 50           # Size of the box in 3D cartesian
saveAscii = True  # True to save the particle trajectories in ASCII
autozoom  = True
################################################################################

Dt = cStep*dt

fig = plt.figure()
ax = fig.add_subplot(111,projection='3d')
if (autozoom == False):
	ax.set_xlim([0, Lx])
	ax.set_ylim([0, Ly])
	ax.set_zlim([0, Lz])

fd = open(fileids,"rt")
reader = csv.reader(fd)
listId = []
for row in reader:
  idp = []
  idp.append(row[0])
  idp.append(row[1])
  listId.append(idp)
fd.close()

for idp in listId:
	tagR = int(idp[0])
	tagP = int(idp[1])
	print("particles %d %d" %(tagR, tagP))
	hisX = []
	hisY = []
	hisZ = []
	hisW = []
	
	ind = 0
	for i in range(cStart,cEnd,cStep):
		filename = "%s_%06d.h5" % (basename, i)
		keytagR = "Step#0/tagR_%d" % spec
		keytagP = "Step#0/tagP_%d" % spec
		keyX = "Step#0/x_%d" % spec
		keyY = "Step#0/y_%d" % spec
		keyZ = "Step#0/z_%d" % spec
		keyW = "Step#0/w_%d" % spec
		print "   ",filename
	
		file  = h5py.File(filename, driver='core')
		X = file.get(keyX)
		Y = file.get(keyY)
		Z = file.get(keyZ)
		W = file.get(keyW)
		aTagR = file.get(keytagR)[()]
		aTagP = file.get(keytagP)[()]
		found = False
		for p in range(X.shape[0]):
			if (aTagR[p] == tagR and aTagP[p] == tagP):
				hisX.append(X[p])
				hisY.append(Y[p])
				hisZ.append(Z[p])
				if (cyl):
					hisW.append(W[p]);
				found = True
				break
	 
		if (found == False):
			break 
	
	if (cyl):
		theta = 0;
		for i in range(len(hisX)):
			r = hisX[i]
			z = hisY[i]
			theta += Dt*hisW[i]/r
			hisX[i] = Lx/2 + np.cos(theta)*r
			hisY[i] = Ly/2 + np.sin(theta)*r
			hisZ[i] = z
	
	#ax.plot(hisX,hisY,hisZ,'b.')
	ax.plot(hisX,hisY,hisZ,'-')
	if (saveAscii):
		partfd = "part_%d_%d.txt" % (tagR, tagP)
		fda = open(partfd,"w")
		for i in range(len(hisX)):
			msg = "%13.6e %13.6e %13.6e\n" % (hisX[i], hisY[i], hisZ[i])
			fda.write(msg)
		fda.close()
	
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
plt.show()

