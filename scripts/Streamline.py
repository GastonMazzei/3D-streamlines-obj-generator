import sys, os

PATH = f'workdir'

if 'Mesh.txt' not in os.listdir(f'{PATH}'): 
	print('[E]: Mesh not found')
	sys.exit(1)


from StreamlineAux import (
		get_fields,
		get_streamline,
		plot_streamline,
		write_box,
			)

import numpy as np

import matplotlib.pyplot as plt

from scipy.optimize import minimize
from scipy.spatial.transform import Rotation as R
import open3d as o3d




def make_v_from_two_points(alpha,x0,xf,N=100):
	return (np.linspace(x0,xf,N) - np.asarray(x0)) * (np.linspace(x0,xf,N) - np.asarray(xf)) * alpha


# Read from boxpoints_info.txt the first three numbers, i.e.
with open(f'{PATH}/Mesh.txt','r') as f:
	data = f.readlines()
Lx = float(data[0]) # In meters
Ly = float(data[1])
Lz = float(data[2])

# Other params
Offx = float(data[3])
Offy = float(data[4])
Offz = float(data[5])
Qx = float(data[6])
Qy = float(data[7])
Qz = float(data[8])
Qw = float(data[9])
Rot = R.from_quat([Qx,Qy,Qz,Qw]).as_matrix()
Off = np.asarray([Offx,Offy,Offz])

# Record box
write_box(PATH,Lx,Ly,Lz)

u,v,w,x,y,z,idX,idY,idZ,Lx,Ly,Lz,dx,dy,dz =  get_fields(Lx,Ly,Lz)
velocity_field_magnitude = np.sqrt(u.flatten() * u.flatten() + v.flatten() * v.flatten() + w.flatten() * w.flatten())
keep_trying = True
while keep_trying:
	X,Y = get_streamline(u,v,w,x,y,z,idX,idY,idZ,Lx,Ly,Lz,dx,dy,dz)
	if X.shape[0] > 0:
		keep_trying = False
DEBUG = False
if DEBUG:
	ax = plt.figure().add_subplot(projection='3d')
	ax.quiver(x, y, z, u, v, w, length=0.1, normalize=True, alpha=0.1)
	plot_streamline(X,'r', ax)
	plt.show()



# Function to fit
def f(a,b,c,d, _X, _X_2, _X_3):
	"""
	f(a,b,c,d) is a pol of deg 3 w/ offset param
	"""
	return a * _X_3 + b * _X_2 + c * _X + d
def cost(params, Y,X,X_2,X_3):
	"""
	f(X,a,b,c,d)=y' should yield Y
	"""
	a,b,c,d  = params
	return np.linalg.norm(f(a,b,c,d,X,X_2,X_3) - Y) 	


# Resample it with more density
N=200
toy_support = np.linspace(0,1,X.shape[0])
toy_support_2 = toy_support * toy_support
toy_support_3 = toy_support_2 * toy_support
denseX = np.linspace(0,1,N)
denseX_2 = denseX * denseX
denseX_3 = denseX_2 * denseX


# Smooth parametrization of domain at each point
Xpred = []
for i in range(3):
	res = minimize(lambda l: cost(l,X[:,i],toy_support,toy_support_2,toy_support_3), [1,1,1,1], method='Nelder-Mead', tol=1e-6)		
	Xpred += [f(res.x[0], res.x[1], res.x[2], res.x[3], denseX, denseX_2, denseX_3).reshape(-1,1)]
Xpred = np.concatenate(Xpred,1)
print('finished fitting the values of the points of interest')


# Smooth parametrization of arrows at each point
v = []
for i in range(3):
	res = minimize(lambda l: cost(l,Y[:,i],toy_support,toy_support_2,toy_support_3), [1,1,1,1], method='Nelder-Mead', tol=1e-6)			
	v += [f(res.x[0], res.x[1], res.x[2], res.x[3], denseX, denseX_2, denseX_3).reshape(-1,1)]
v = np.concatenate(v,1)
print('finished fitting the values of the field in the points of interest')


# Apply rotation and translation !CANCELLED, IT IS DONE INSIDE THE SPECTACLES!
#Xpred = Xpred @ Rot + Off
#denseX[:,2] *= -1

# Define a bulk operation to make the curves three dimensionals
L = 0.01
L_2 = L / 2
sorrounders = np.asarray([ # square section
	-L_2, L_2, 0,
	-L_2, -L_2, 0,
	L_2, L_2, 0,
	L_2, -L_2, 0,	
		]).reshape(-1,3)
def bulk(vec, normal_vec, alpha_dim1, alpha_dim2):
	"""
	Give some bulk to the 1D curve
	"""
	# Normalize the field directions
	normal_vec_normalized = normal_vec / np.linalg.norm(normal_vec, axis=1).reshape(-1,1)
	# Set up the initial perturbation 
	initial_pos = np.cross(np.asarray([1,1,1]), normal_vec_normalized)
	initial_pos /= np.linalg.norm(initial_pos, axis=1).reshape(-1,1)
	# Set up the second perturbation
	second_pos = np.cross(initial_pos, normal_vec_normalized)
	second_pos /= np.linalg.norm(second_pos, axis=1).reshape(-1,1)
	# Set up the third and fourth perturbations
	third_pos = -initial_pos
	fourth_pos = -second_pos
	# Four more perturbations are possible
	four_more = [initial_pos + second_pos, second_pos + third_pos, fourth_pos + third_pos, fourth_pos + initial_pos]
	four_more = [x / np.linalg.norm(x, axis=1).reshape(-1,1) * (alpha_dim1 ** 2 + alpha_dim2 ** 2) + vec for i,x in enumerate(four_more)]
	result = [
			#vec, 
			alpha_dim1 * initial_pos + vec,  #0
			alpha_dim2 * second_pos + vec,  #1
			alpha_dim1 * third_pos + vec,  #2
			alpha_dim2 * fourth_pos + vec] #3
	result += four_more # 4,5,6,7
	return np.asarray(result)

# Compute the bulk
print('bulk computing section')
sufficient = False
target_abs = 0.05 # min([Lx,Ly,Lz,4]) * 0.02
target_rel = (0.8,1.3)
alpha_dim1 = 1
alpha_dim2 = 1
while not sufficient:
	all_v = bulk(Xpred, v, alpha_dim1, alpha_dim2)
	# Criteria with the maximum
	#abs_disp = max ( [ max( np.linalg.norm(all_v[0]-all_v[2],axis=1) ), max( np.linalg.norm(all_v[1]-all_v[3],axis=1) ) ] )
	# Criteria with the average
	abs_disp = (np.mean( np.linalg.norm(all_v[0]-all_v[2],axis=1) ) + np.mean( np.linalg.norm(all_v[1]-all_v[3],axis=1) ) ) / 2
	if abs_disp > target_abs:
		alpha_dim1 *= 0.5
		alpha_dim2 *= 0.5
		print(f'alpha was not sufficient, shrinking to: {alpha_dim1}')
		sufficient = False
	else:
		sufficient = True
print(all_v.shape)


# Display
DISPLAY = False
if DISPLAY:
	fig = plt.figure()
	ax = fig.add_subplot(projection='3d')
	for j,v in enumerate(all_v):
		for i in range(len(v)-1):
			ax.plot([v[i+1,0],v[i,0]], [v[i+1,2],v[i,2]], [v[i+1,1],v[i,1]], c=['r','g','b','k','y','magenta','orange','violet','gray'][j])
	ax.set_xlim(-Lx/2,Lx/2)
	ax.set_ylim(-Lz/2,Lz/2)
	ax.set_zlim(-Ly/2,Ly/2)
	plt.show()

# Give colors
def speed_to_color(vec: np.ndarray, absmax = np.inf, absmin = -np.inf):
	"""
	Convert magnitude to color
	"""
	# Test: hardcoded red
	#COLORING = [1,0,0] # yellow
	#return [COLORING for _ in range(len(vec))]
	# end of test
	EPS = 1e-5
	magnitudes = np.linalg.norm(vec, axis=1)
	if absmin == -np.inf and absmax == np.inf:
		M,m = max(magnitudes), min(magnitudes)
	else:
		M = absmax
		m = absmin
	width = M-m
	if width<EPS: 
		# Velocity is constant, choose some coloring
		COLORING = [1,1,0] # yellow
		return [COLORING for _ in range(len(magnitudes))]
	width_2 = width / 2
	scaled_magnitudes = magnitudes - m
	result = []
	for val in scaled_magnitudes:
		if val < width_2:
			result += [[val / width_2, 1, 0]]
		elif val >= width_2:
			result += [[1, 1 - (val - width_2) / width_2, 0]]
		else:
			print(val, width)
			raise Exception("Color out of bounds!")		
	return result

absmax, absmin = np.max(velocity_field_magnitude), np.min(velocity_field_magnitude)
colors = np.asarray(speed_to_color(v, absmax, absmin))
print(colors.shape)

# Save results
total_mesh = o3d.geometry.TriangleMesh()
colors = np.concatenate([colors,colors],0)
colorsL_2 = colors.shape[0] // 2
print('entering the meshing section')
SPLITTED_INTO = 3
Nsections = int(N * 0.4) # 40
radii = [min([max([dx * 2, dy * 2, dz * 2 ]),0.4])]
for Q in [(0,4),(4,1),(1,5),(5,2),(2,6),(6,3),(3,7),(7,0)]:
	a,b = Q
	section_length = int(all_v[[a,b],:].shape[1] / Nsections)
	for split in range(Nsections):
		low,high = split * section_length, (split+1) * section_length
		if low>1:
			low -= 2
		high += 2 
		if high > all_v[[a,b],:].shape[1]:
			high = all_v[[a,b],:].shape[1]
		points = np.concatenate(all_v[[a,b], low:high, :],0)
		#print(points.shape)
		pcd = o3d.geometry.PointCloud()
		pcd.points = o3d.utility.Vector3dVector(points)
		local_colors = np.asarray(colors[low:high].tolist() + colors[low + colorsL_2 : high + colorsL_2].tolist())
		pcd.colors = o3d.utility.Vector3dVector(local_colors)
		#o3d.io.write_point_cloud(f"./data{a}-{b}-{split}.ply", pcd)
		# Pcd to Mesh
		pcd.estimate_normals(search_param=o3d.geometry.KDTreeSearchParamHybrid(radius=max([dx * 5, dy * 5, dz * 5]), max_nn=12))
		if True:
			rec_mesh = o3d.geometry.TriangleMesh.create_from_point_cloud_ball_pivoting(pcd, o3d.utility.DoubleVector(radii))
		elif False:
			rec_mesh, densities = o3d.geometry.TriangleMesh.create_from_point_cloud_poisson(
											pcd, depth=9)
		total_mesh += rec_mesh
		#o3d.io.write_triangle_mesh(f"./dataMesh{a}-{b}-{split}.obj", rec_mesh)

try:
	os.mkdir(f"{PATH}/streamlines")
except: pass

THRESHOLD = all_v.shape[0] * 0.5
if np.asarray(total_mesh.triangles).shape[0] < THRESHOLD:
	print('calling myself recursively')
	os.system(f'python3 ./scripts/Streamline.py {sys.argv[1]}')
else:
	latest_num = max([1] + [int(x.split('.')[0].split('streamLine')[1]) for x in os.listdir(f"{PATH}/streamlines")])
	o3d.io.write_triangle_mesh(f"{PATH}/streamlines/streamLine{latest_num+1}.obj", total_mesh)
	if int(sys.argv[1])>1:
		os.system(f'python3 ./scripts/Streamline.py {int(sys.argv[1]) - 1}')




