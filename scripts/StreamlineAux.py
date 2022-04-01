import matplotlib.pyplot as plt
import numpy as np
from math import floor
import os
import open3d as o3d

def write_box(PATH,Lx,Ly,Lz):
	if 'boxPoints.ply' in os.listdir(f'{PATH}'): 
		return

	points = np.asarray([
		Lx/2,Ly/2,Lz/2,
		-Lx/2,Ly/2,Lz/2,
		Lx/2,-Ly/2,Lz/2,
		-Lx/2,-Ly/2,Lz/2,
		Lx/2,Ly/2,-Lz/2,
		-Lx/2,Ly/2,-Lz/2,
		Lx/2,-Ly/2,-Lz/2,
		-Lx/2,-Ly/2,-Lz/2	
		]).reshape((-1,3))
	pcd = o3d.geometry.PointCloud()
	pcd.points = o3d.utility.Vector3dVector(points)
	o3d.io.write_point_cloud(f"{PATH}/boxPoints.ply", pcd)

def get_fields(Lx=4,Ly=3,Lz=5):
	# Make the grid
	dx = float(Lx / 10)
	dy = float(Ly / 10)
	dz = float(Lz / 10)
	x, y, z = np.meshgrid(np.arange(-Lx/2, Lx/2, dx),
		              np.arange(-Ly/2, Ly/2, dy),
		              np.arange(-Lz/2, Lz/2, dz))
	idX = x * 0 + 1
	idY = y * 0 + 1
	idZ = z * 0 + 1

	# Make the direction data for the arrows
	u = 1 * np.sin(np.pi * x) * np.cos(np.pi * y) + 0.1
	v = 1 * -np.cos(np.pi * x) * np.sin(np.pi * y) + 0.1
	w = idZ * 0.1

	return u,v,w,x,y,z,idX,idY,idZ,Lx,Ly,Lz,dx,dy,dz

def getFieldsFromCSV(filename):
	import pandas as pd
	df = pd.read_csv(filename)
	Lx = df['x'].max() - df['x'].min()
	Ly = df['y'].max() - df['y'].min()
	Lz = df['z'].max() - df['z'].min()

	# Make the grid
	# pend create!

	return u,v,w,x,y,z,idX,idY,idZ,Lx,Ly,Lz,dx,dy,dz


def get_streamline(u,v,w,x,y,z,idX,idY,idZ,Lx,Ly,Lz,dx,dy,dz):
	F = 2
	pos = [np.random.uniform(-Lx/F,Lx/F), np.random.uniform(-Ly/F,Ly/F), np.random.uniform(-Lz/F,Lz/F)]
	vel = []
	EPOCHS = 100000
	dt = 0.02
	for _ in range(EPOCHS):
		# Compute indices
		i = floor((pos[-3] + Lx/2) / dx)
		j = floor((pos[-2] + Ly/2) / dy)
		k = floor((pos[-1] + Lz/2) / dz)
		try:
			vel += [u[i,j,k], v[i,j,k], w[i,j,k]]
		except Exception as ins:
			print(f'[I]: unknown error, {ins.args}')
			break
		update = [dt * vel[-3], dt * vel[-2], dt * vel[-1]]
		pos += [pos[-3] + update[0], pos[-2] + update[1], pos[-1] + update[2]]
		if (pos[-3] > Lx/2) or (pos[-3] < -Lx/2) or (pos[-2] > Ly/2) or (pos[-2] < -Ly/2) or (pos[-1] > Lz/2) or (pos[-1] < -Lz/2):
			print('[I]: ending the streamline because of out of bounds')
			break
	pos = pos[:-3]
	vel = np.asarray(vel).reshape(-1,3)
	pos = np.asarray(pos).reshape(-1,3)
	return pos,vel

def plot_streamline(pos, c, ax):
	ax.plot(pos[:,0], pos[:,1], pos[:,2], c=c)


if __name__ == '__main__':

	PATH = f'workdir'
	with open(f'{PATH}/Mesh.txt','r') as f:
		data = f.readlines()
	Lx = float(data[0]) # In meters
	Ly = float(data[1])
	Lz = float(data[2])

	#ax = plt.figure().add_subplot(projection='3d')
	u,v,w,x,y,z,idX,idY,idZ,Lx,Ly,Lz,dx,dy,dz =  get_fields(Lx,Ly,Lz)
	A = 0.6
	#ax.quiver(x, y, z, u, v, w, length=0.1, normalize=True, alpha=A, color='r')
	#plt.show()
	ax = plt.figure().add_subplot(projection='3d')
	ax.quiver(x, y, z, u, v, w, length=0.1, normalize=True, alpha=A, color='r')
	S = []
	for _ in range(30):
		S += [get_streamline(u,v,w,x,y,z,idX,idY,idZ,Lx,Ly,Lz,dx,dy,dz)]
		print(S[-1][0].shape, S[-1][1].shape)
		plot_streamline(S[-1][0],
				#['r','y','magenta','orange','purple'][_ % 5],
				'k',
				 ax)
	ax.set_xlabel('X')
	ax.set_ylabel('Y')
	ax.set_zlabel('Z')
	plt.show()


