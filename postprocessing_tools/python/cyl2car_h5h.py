################################################################################
################################################################################
####                                                                        ####
####  Convert a field file computed in 2D in cylindrical coordinates into   ####
####  a file in 3D in Cartesian coordinates.                                ####
####                                                                        ####
####  Version for H5hut files                                               ####
####                                                                        ####
####  Only the following magnitudes are converted:                          ####
####     Ex                                                                 ####
####     Ey                                                                 ####
####     Ez                                                                 ####
####     Bx                                                                 ####
####     By                                                                 ####
####     Bz                                                                 ####
####     Bx_ext                                                             ####
####     By_ext                                                             ####
####     Bz_ext                                                             ####
####     Jx_0                                                               ####
####     Jy_0                                                               ####
####     Jz_0                                                               ####
####     rho_0                                                              ####
####                                                                        ####
####                                                                        ####
####                                                                        ####
####  Diego Gonzalez                                                        ####
####  August 2017                                                           ####
####                                                                        ####
################################################################################
################################################################################

import numpy as np
import h5py
import sys

# Set the basename, range and final filename
def readField(field, file):
	key = 'Step#0/Block/'+ field + '/0'
	data  = file.get(key)
	nnx = data.shape[2]
	nny = data.shape[1]
	nnz = data.shape[0]
	
  # This has to be done to force hdf5 to laod the data into memory.
  # Otherwise each time data[i][j][k] is needed it is read from the file.
	dst = np.empty((nnz, nny, nnx))
	dst[:][:][:] = data[:][:][:] 
	return dst

def computeValuesVector(Bx_3d, By_3d, Bz_3d, Bx, By, Bz, i, j, k, i0, j0, i1, j1, w00, w01, w10, w11): 
	ValR = w00*Bx[0][j0][i0] + w01*Bx[0][j0][i1] + w10*Bx[0][j0][i1] + w11*Bx[0][j1][i1]
	ValZ = w00*By[0][j0][i0] + w01*By[0][j0][i1] + w10*By[0][j0][i1] + w11*By[0][j1][i1]
	ValT = w00*Bz[0][j0][i0] + w01*Bz[0][j0][i1] + w10*Bz[0][j0][i1] + w11*Bz[0][j1][i1]
	alpha = np.arctan2(z, x)
	cA = np.cos(alpha)
	sA = np.sin(alpha)
	Bx_3d[k][j][i] = cA*ValR - sA*ValT
	By_3d[k][j][i] = ValZ
	Bz_3d[k][j][i] = sA*ValR + cA*ValT

def computeValuesScaler(A_3d, A, i, j, k, i0, j0, i1, j1, w00, w01, w10, w11): 
	val = w00*A[0][j0][i0] + w01*A[0][j0][i1] + w10*A[0][j0][i1] + w11*A[0][j1][i1]
	A_3d[k][j][i] = val


#####################################################################################
##    	MAIN
#####################################################################################

if (len(sys.argv) < 2):
	print("")
	print("  Use: " + sys.argv[0] + " <file> (<file> ...)")
	sys.exit(0)

nf = len(sys.argv)-1

for i in range(nf):
	old_file = str(sys.argv[1])
	#new_file = base + "_3D.h5"
	if (old_file[len(old_file)-17:len(old_file)-9] != "-Fields_"):
		print("ERROR: The name of the files (" + old_file + ") does not have the usual format: SimulationName-Fields_000000.h5")
		sys.exit(0)

	new_file = old_file[:len(old_file)-17] + "_3D" + old_file[len(old_file)-17:len(old_file)]

	print "File ", old_file , " -> ", new_file

  # Input file in Cylindrical coordinates
	fileIn  = h5py.File(old_file, driver='core');


	fields = fileIn.get("/Step#0/Block")

  # Calculation of the new grid in 3D
	print("   Reading data ...")
	Bx_ext = readField("Bx_ext",fileIn)
	By_ext = readField("By_ext",fileIn)
	Bz_ext = readField("Bz_ext",fileIn)
	Bx = readField("Bx",fileIn)
	By = readField("By",fileIn)
	Bz = readField("Bz",fileIn)
	Ex = readField("Ex",fileIn)
	Ey = readField("Ey",fileIn)
	Ez = readField("Ez",fileIn)
	Jx = readField("Jx_0",fileIn)
	Jy = readField("Jy_0",fileIn)
	Jz = readField("Jz_0",fileIn)
	rho= readField("rho_0",fileIn)


	nx_2D = Bx.shape[2]
	ny_2D = Bx.shape[1]
	nz_2D = Bx.shape[0]

	nx = (Bx.shape[2]-1)*2 + 1
	ny = Bx.shape[1]
	nz = nx
	x0 = nx/2.0
	z0 = nz/2.0

	Bx_ext_3d = np.empty((nz,ny,nx))
	By_ext_3d = np.empty((nz,ny,nx))
	Bz_ext_3d = np.empty((nz,ny,nx))
	Bx_3d = np.empty((nz,ny,nx))
	By_3d = np.empty((nz,ny,nx))
	Bz_3d = np.empty((nz,ny,nx))
	Ex_3d = np.empty((nz,ny,nx))
	Ey_3d = np.empty((nz,ny,nx))
	Ez_3d = np.empty((nz,ny,nx))
	Jx_3d = np.empty((nz,ny,nx))
	Jy_3d = np.empty((nz,ny,nx))
	Jz_3d = np.empty((nz,ny,nx))
	rho_3d = np.empty((nz,ny,nx))

	ntn = nx*ny*nz
	n = 0

  # Interpolation from 2D cylindrical to 3D cartesian
	for k in range(0,nz-1):
		for j in range(0,ny-1):
			for i in range(0,nx-1):
				n = n + 1
				if (np.mod(n,100) == 0):
					prog = float(n)/ntn*100
					msg = "\r   Processing ( %.1f %% )" % (prog)
					sys.stdout.write(msg)
				x = float(i) - x0
				y = float(j)
				z = float(k) - z0;
				r = np.sqrt(z*z + x*x)
				# Interpolation
				j0 = int(y)
				j1 = j0 + 1
				i0 = int(r)
				i1 = i0 + 1
				outofdom = False

				if (r > nx_2D-1):
					#print "Out of domain"	
					outofdom = True
				else:
					#print "r=",r,j0,j1,i0,i1
					w00 = (j1 - y)*(i1 - r)
					w01 = (y - j0)*(i1 - r)
					w10 = (j1 - y)*(r - i0)
					w11 = (y - j0)*(r - i0)
					if (i1 > nx_2D - 1): i1 = i0
					if (j1 > ny_2D - 1): j1 = j0
				
				if (outofdom == False):
					computeValuesVector(Bx_3d, By_3d, Bz_3d, Bx, By, Bz, i, j, k, i0, j0, i1, j1, w00, w01, w10, w11)
					computeValuesVector(Bx_ext_3d, By_ext_3d, Bz_ext_3d, Bx_ext, By_ext, Bz_ext, i, j, k, i0, j0, i1, j1, w00, w01, w10, w11)
					computeValuesVector(Ex_3d, Ey_3d, Ez_3d, Ex, Ey, Ez, i, j, k, i0, j0, i1, j1, w00, w01, w10, w11)
					computeValuesVector(Jx_3d, Jy_3d, Jz_3d, Jx, Jy, Jz, i, j, k, i0, j0, i1, j1, w00, w01, w10, w11)
					computeValuesScaler(rho_3d, rho, i, j, k, i0, j0, i1, j1, w00, w01, w10, w11)

	msg = "\r   Processing ( 100.0 %% )\n"
	sys.stdout.write(msg)

	print("   Writting data ...")
  # Output file in cartesian geometry
	fileOut = h5py.File(new_file,"w");
	maingrp = fileOut.create_group('Step#0/Block')

	grp = maingrp.create_group("Bx") 
	grp.create_dataset('0', data=Bx_3d)

	grp = maingrp.create_group("By") 
	grp.create_dataset('0', data=By_3d)

	grp = maingrp.create_group("Bz") 
	grp.create_dataset('0', data=Bz_3d)

	grp = maingrp.create_group("Bx_ext") 
	grp.create_dataset('0', data=Bx_ext_3d)

	grp = maingrp.create_group("By_ext") 
	grp.create_dataset('0', data=By_ext_3d)

	grp = maingrp.create_group("Bz_ext") 
	grp.create_dataset('0', data=Bz_ext_3d)

	grp = maingrp.create_group("Ex") 
	grp.create_dataset('0', data=Ex_3d)

	grp = maingrp.create_group("Ey") 
	grp.create_dataset('0', data=Ey_3d)

	grp = maingrp.create_group("Ez") 
	grp.create_dataset('0', data=Ez_3d)

	grp = maingrp.create_group("Jx_0") 
	grp.create_dataset('0', data=Jx_3d)

	grp = maingrp.create_group("Jy_0") 
	grp.create_dataset('0', data=Jy_3d)

	grp = maingrp.create_group("Jz_0") 
	grp.create_dataset('0', data=Jz_3d)

	grp = maingrp.create_group("rho_0") 
	grp.create_dataset('0', data=rho_3d)

	fileOut.close()
	print("   DONE")

  

