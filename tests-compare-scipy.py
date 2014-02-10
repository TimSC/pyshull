
import matplotlib.pyplot as plt
import pyshull, sys, time, pickle, random
import numpy as np
import scipy.spatial as spatial

def CompareTriangleLists(triangles1, triangles2):
	listOfTuples1 = [tuple(tri) for tri in triangles1]
	listOfTuples2 = [tuple(tri) for tri in triangles2]
	count = 0
	for tri in listOfTuples1:
		match = tri in listOfTuples2
		if not match:
			tri2 = (tri[1],tri[2],tri[0])
			match = tri2 in listOfTuples2
		if not match:
			tri2 = (tri[2],tri[0],tri[1])
			match = tri2 in listOfTuples2
		if not match:
			tri2 = (tri[2],tri[1],tri[0])
			match = tri2 in listOfTuples2
		if not match:
			tri2 = (tri[1],tri[0],tri[2])
			match = tri2 in listOfTuples2
		if not match:
			tri2 = (tri[0],tri[2],tri[1])
			match = tri2 in listOfTuples2

		if match:
			count += 1
	return float(count) / float(len(triangles1))

if __name__ == "__main__":

	n = 50
	if len(sys.argv) >= 2:
		n = int(sys.argv[1])
	problemCount = 0	

	while 1:
	
		#Generate random points
		pts = np.random.rand(n, 2)

		#plt.plot(pts[:,0], pts[:,1], 'x')
		#plt.show()

		#Align some of the points
		for i in range(pts.shape[0]):
			for j in range(pts.shape[1]):
				if random.randint(0,1):
					pts[i,j] = random.randint(0,10) / 10.

		startTime = time.time()
		triangles = pyshull.PySHull(pts)
		print "pyqhull Processed", n, "points in", time.time() - startTime, "sec"

		startTime = time.time()
		triangles2 = spatial.Delaunay(pts).simplices
		print "scipy Processed", n, "points in", time.time() - startTime, "sec"

		print "Problems found", problemCount
		
		#print triangles
		compare = CompareTriangleLists(triangles, triangles2)
		
		if compare < 1.:
			print "Problem detected", compare
			plt.clf()
			plt.subplot(211)
			for tri in triangles:
				tri2 = list(tri[:])
				tri2.append(tri[0])
				plt.plot(pts[tri2,0], pts[tri2,1])

			plt.plot(pts[:,0], pts[:,1], 'x')

			plt.subplot(212)
			for tri in triangles2:
				tri2 = list(tri[:])
				tri2.append(tri[0])
				plt.plot(pts[tri2,0], pts[tri2,1])

			plt.plot(pts[:,0], pts[:,1], 'x')
			plt.savefig("problem{0}.svg".format(problemCount))
			
			pickle.dump(pts, open("problem{0}.dat".format(problemCount),"wb"), protocol=-1)

			problemCount += 1

