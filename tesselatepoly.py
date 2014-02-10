
import pyshull, math
import scipy.spatial as spatial

def PolyToTriangles(poly, holes = []):
	
	pass

def CheckIfFlippable(pts, triangles, triNum1, triNum2):
	tri1 = triangles[triNum1]
	tri2 = triangles[triNum2]
	commonEdge = pyshull.HasCommonEdge(tri1, tri2)
	if commonEdge is None:
		return False
	triInd1, triInd2 = commonEdge
	#print "original ind", tri1, tri2

	#Reorder nodes so the common edge is the first two verticies
	triOrdered1 = (tri1[triInd1[0]], tri1[triInd1[1]], tri1[triInd1[2]])
	triOrdered2 = (tri2[triInd2[0]], tri2[triInd2[1]], tri2[triInd2[2]])

	#print "reordered ind", triOrdered1, triOrdered2

	flipTri1 = (triOrdered1[2], triOrdered2[2], triOrdered1[1])
	flipTri2 = (triOrdered2[2], triOrdered1[2], triOrdered1[0])

	#print "flipped", flipTri1, flipTri2

	rh1 = pyshull.RightHandedCheck(pts, flipTri1[1], flipTri1[2], flipTri1[0])
	rh2 = pyshull.RightHandedCheck(pts, flipTri2[1], flipTri2[2], flipTri2[0])

	if rh1 <= 0. or rh2 <= 0.:
		return False

	return True	

def PerformFlipOfTriPair(pts, triangles, triNum1, triNum2):
	tri1 = triangles[triNum1]
	tri2 = triangles[triNum2]
	commonEdge = pyshull.HasCommonEdge(tri1, tri2)
	if commonEdge is None:
		return False
	triInd1, triInd2 = commonEdge
	#print "original ind", tri1, tri2

	#Reorder nodes so the common edge is the first two verticies
	triOrdered1 = (tri1[triInd1[0]], tri1[triInd1[1]], tri1[triInd1[2]])
	triOrdered2 = (tri2[triInd2[0]], tri2[triInd2[1]], tri2[triInd2[2]])

	#print "reordered ind", triOrdered1, triOrdered2

	flipTri1 = (triOrdered1[2], triOrdered2[2], triOrdered1[1])
	flipTri2 = (triOrdered2[2], triOrdered1[2], triOrdered1[0])

	triangles[triNum1] = flipTri1
	triangles[triNum2] = flipTri2

if __name__=="__main__":
	import matplotlib.pyplot as plt
	import numpy as np
	#pts = [(2.,1.),(4.,5.),(2.,0.),(0.,5.)]
	pts = [(0.,0.),(1.,0.),(1.,1.),(2.,1.),(2.,2.),(4.,2.),(4.,4.),(2.,4.),(2.,3.),(1.,3.),(1.,2.),(0.,2.)]
	try:
		triangles = pyshull.PySHull(pts)
	except Exception as err:
		print err
		triangles = err.triangles

	if 1:
		#Flip pairs of triangles while avoiding overlaps
		for i in range(len(triangles)):
			for j in range(len(triangles)):
				if j >= i: continue
				flipable = CheckIfFlippable(pts, triangles, i, j)
				#print i, j, flipable
				if flipable:
					PerformFlipOfTriPair(pts, triangles, i, j)

	#triangles = pyshull.FlipTriangles(pts, triangles)
	#triangles = spatial.Delaunay(pts).simplices

	pts = np.array(pts)
	#plt.plot(pts[:,0], pts[:,1])
	#plt.show()

	for tri in triangles:
		tri2 = list(tri[:])
		tri2.append(tri[0])
		plt.plot(pts[tri2,0], pts[tri2,1])
	plt.show()

