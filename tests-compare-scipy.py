from __future__ import print_function
import matplotlib.pyplot as plt
import pyshull, sys, time, pickle, random
import numpy as np
import scipy.spatial as spatial

def CompareTriangleLists(triangles1, triangles2):
	listOfTuples1 = [tuple(tri) for tri in triangles1]
	listOfTuples2 = [tuple(tri) for tri in triangles2]
	count = 0
	probIndex = []
	for triNum, tri in enumerate(listOfTuples1):
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
		else:
			probIndex.append(triNum)
	return float(count) / float(len(triangles1)), probIndex

def HeronsFormula(pts, tri):

	a = pyshull.CalcDist(pts[tri[0]], pts[tri[1]])
	b = pyshull.CalcDist(pts[tri[1]], pts[tri[2]])
	c = pyshull.CalcDist(pts[tri[2]], pts[tri[0]])

	#https://en.wikipedia.org/wiki/Heron%27s_formula#Numerical_stability
	x1 = (a+(b+c))
	x2 = (c-(a-b))
	x3 = (c+(a-b))
	x4 = (a+(b-c))
	area = 0.25 * ((x1*x2*x3*x4) ** 0.5)
	return area

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

		#Remove duplicates
		pts = pyshull.RemoveDuplicatePoints(pts)
		pts = np.array(pts)

		startTime = time.time()
		triangles = pyshull.PySHull(pts)
		print("pyqhull Processed", n, "points in", time.time() - startTime, "sec")

		startTime = time.time()
		triangles2 = spatial.Delaunay(pts).simplices
		print("scipy Processed", n, "points in", time.time() - startTime, "sec")
		
		for tri in triangles:
			area = HeronsFormula(pts, tri)
			if area == 0.:
				print("Problem: Zero size triangle")
				pickle.dump(pts, open("problem{0}.dat".format(problemCount),"wb"), protocol=-1)
				problemCount += 1

		#print(triangles)
		compare, probIndex = CompareTriangleLists(triangles, triangles2)
		compare2, probIndex2 = CompareTriangleLists(triangles2, triangles)
		
		if compare + compare2 < 2.:
			print("Problem detected", compare, compare2, len(triangles), len(triangles2))
			plt.clf()
			plt.subplot(211)
			plt.xlim([-0.1, 1.1])
			plt.ylim([-0.1, 1.1])
			for triNum, tri in enumerate(triangles):
				tri2 = list(tri[:])
				tri2.append(tri[0])
				col = 'g-'
				if triNum in probIndex: 
					print("p1", triNum, tri, HeronsFormula(pts, tri))
					#print("z1", pts[tri2,0], pts[tri2,1])
					col = "r-"
				plt.plot(pts[tri2,0], pts[tri2,1], col)

			plt.plot(pts[:,0], pts[:,1], 'x')

			plt.subplot(212)
			plt.xlim([-0.1, 1.1])
			plt.ylim([-0.1, 1.1])
			for triNum, tri in enumerate(triangles2):
				tri2 = list(tri[:])
				tri2.append(tri[0])
				col = 'g-'
				if triNum in probIndex2: 
					print("p2", triNum, tri, HeronsFormula(pts, tri))
					#print("z2", pts[tri2,0], pts[tri2,1])
					col = "r-"

				plt.plot(pts[tri2,0], pts[tri2,1], col)

			plt.plot(pts[:,0], pts[:,1], 'x')
			plt.savefig("problem{0}.svg".format(problemCount))
			
			pickle.dump(pts, open("problem{0}.dat".format(problemCount),"wb"), protocol=-1)

			problemCount += 1

		print("Problems found", problemCount)

