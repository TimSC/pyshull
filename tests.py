
import matplotlib.pyplot as plt
import pyshull, sys, time
import numpy as np

if __name__ == "__main__":
	
	n = 50
	if len(sys.argv) >= 2:
		n = int(sys.argv[1])

	pts = np.random.rand(n, 2)

	startTime = time.time()
	triangles = pyshull.PySHull(pts)
	print "Processed", n, "points in", time.time() - startTime, "sec"

	for tri in triangles:
		tri2 = list(tri[:])
		tri2.append(tri[0])
		plt.plot(pts[tri2,0], pts[tri2,1])

	#plt.plot(pts[:,0], pts[:,1], 'x')
	plt.show()

	
