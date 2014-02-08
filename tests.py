
import matplotlib.pyplot as plt
import pyshull
import numpy as np

if __name__ == "__main__":
	
	pts = np.random.rand(50, 2)

	triangles = pyshull.PySHull(pts)

	for tri in triangles:
		tri2 = list(tri[:])
		tri2.append(tri[0])
		plt.plot(pts[tri2,0], pts[tri2,1])

	#plt.plot(pts[:,0], pts[:,1], 'x')
	plt.show()

	
