################################################################################
####
####  	VISUALIZATION OF PARTICLE POSITION FORM THE RESULTS OF ECSIM 
####    (IN CYLINDRICAL COORDINATES)
####
####    File type: H5HUT
####
################################################################################
####
####  Diego Gonzalez
####  KU Leuven
####  June 2017
####
################################################################################
################################################################################
import numpy as np
import h5py
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import sys

mode = 0
Y0 = 0
Y1 = 1
if (len(sys.argv) < 3):
	print "Visualizarion of particle position from the results of ECsim in cylindrical coordinates (using H5hut)" 
	print "Use %s <1> <2> (<3> <4>)" % (sys.argv[0])
	print "\t <1>   Name of the particles file (h5hut)"
	print "\t <2>   Mode: 3D  -> three dimensional visualization in cartesian coordinates"
	print "\t             2D  -> two dimensional visualization in cartesian and cylindrical coordinates"
	print "\t             Cut -> 2D visulaization of the particles in a plane between Y0 and Y1."
	print "\t(<3>)  Y0 only used if Cut mode is selected."
	print "\t(<4>)  Y1 only used if Cut mode is selected."

	sys.exit()

elif(len(sys.argv) > 2):
	if (sys.argv[2] == '3D'):
		mode = 1
	elif (sys.argv[2] == '2D'):
		mode = 2
	elif (sys.argv[2] == 'Cut'):
		if (len(sys.argv) < 5):
			print "Mode Cut requires two more arguments (Y0, Y1)"
			sys.exit()
		Y0 = float(sys.argv[3])
		Y1 = float(sys.argv[4])
		mode = 3
	else:
		print "Mode %s is not available" % (sys.argv[2])
		sys.exit()

filename = sys.argv[1]
print filename

with h5py.File(filename,'r') as hf:
	step0 = hf.get('Step#0')
	Q = step0.get('q_0')
	R = step0.get('x_0')
	Z = step0.get('y_0')
	T = step0.get('z_0')

	n = len(R)
	theta = 2*np.pi*np.random.rand(n)
	X3d =  np.cos(theta)*R;
	Y3d = -np.sin(theta)*R
	Z3d = Z

	Xmax = np.amax(R)

	if (mode == 1):
		f = plt.figure()
		ax = f.add_subplot(111, projection='3d')
		ax.set_xlabel('x')
		ax.set_ylabel('y')
		ax.set_zlabel('z')
		ax.scatter(X3d, Y3d, Z3d, marker = 'o', s=0.5, c='r')
		plt.show()
	elif (mode == 2):
		f, axarr = plt.subplots(1,2)
		axarr[0].set_xlim(0, Xmax)
		#axarr[0].set_ylim(0, 75)
		axarr[1].set_xlim(0, Xmax)
		#axarr[1].set_ylim(0, 75)
		axarr[0].set_title('Cylindrical coordinates (2D)')
		axarr[0].set_xlabel('r')
		axarr[1].set_xlabel('x')
		axarr[0].set_ylabel('z')
		axarr[1].set_ylabel('z')
		axarr[1].set_title('Cartesian coordinates (3D)')
		axarr[0].scatter(R, Z, marker = 'o', s=0.1, c='b')
		axarr[1].scatter(X3d, Z3d, marker = 'o', s=0.1, c='r')
		plt.show()
	elif (mode == 3):
		p = len(X3d)
		nX = []
		nY = []
		nZ = []
		for i in range(0,n-1):
			if (Y3d[i] > Y0 and Y3d[i] < Y1):
				nX.append(X3d[i])
				nY.append(Y3d[i])
				nZ.append(Z3d[i])

		print "%d particles in the range y=[%g, %g]" % (len(nX), Y0, Y1)

		plt.figure()
		plt.title('3D cut. y = [%g, %g]' % (Y0, Y1))
		plt.xlabel('x')
		plt.ylabel('z')
		#plt.xlim(0, 37.5)
		#plt.ylim(0, 75)
		plt.scatter(nX, nZ, marker = 'o', s=0.1, c='b')
		plt.gca().set_aspect('equal', adjustable='box')
		plt.show()

		#fileres = "fig.png" % (ind)
		#f.savefig(fileres)


