import math

def CalcTriangleAng(pts, angleCache, pt1, pt2, pt3):

	angId = (pt1, pt2, pt3)
	if angId in angleCache:
		return angleCache[angId]

	#Angle is computed on pt3. pt1 and pt2 define the side opposite the angle
	pt1v = pts[pt1]
	pt2v = pts[pt2]
	pt3v = pts[pt3]
	v31 = (pt1v[0] - pt3v[0], pt1v[1] - pt3v[1])
	v32 = (pt2v[0] - pt3v[0], pt2v[1] - pt3v[1])
	mv31 = (v31[0]**2. + v31[1]**2.) ** 0.5
	mv32 = (v32[0]**2. + v32[1]**2.) ** 0.5
	if mv31 == 0. or mv32 == 0.:
		raise RuntimeError("Angle not defined for zero area triangles")

	v31n = [c / mv31 for c in v31]
	v32n = [c / mv32 for c in v32]
	crossProd = - v31n[0] * v32n[1] + v31n[1] * v32n[0]
	dotProd = v31n[0] * v32n[0] + v31n[1] * v32n[1]
	
	#Limit to valid range
	if dotProd > 1.: dotProd = 1.
	if dotProd < -1.: dotProd = -1.

	#print crossProd < 0., crossProd
	#print math.asin(crossProd), math.acos(dotProd), cosAng
	if crossProd < 0.:
		#Reflex angle detected
		trigAng = 2. * math.pi - math.acos(dotProd)
	else:
		#Acute or obtuse angle
		trigAng = math.acos(dotProd)

	angleCache[angId] = trigAng
	return trigAng

def RightHandedCheck(pts, pt1, pt2, pt3):
	vec21 = (pts[pt1][0] - pts[pt2][0], pts[pt1][1] - pts[pt2][1])
	vec23 = (pts[pt3][0] - pts[pt2][0], pts[pt3][1] - pts[pt2][1])
	return vec21[0] * vec23[1] - vec21[1] * vec23[0]

def EarClipping(poly):
	#Based on Triangulation by Ear Clipping by David Eberly

	workingPoly = range(len(poly))
	angleCache = {}
	triangles = []

	while len(workingPoly) > 3:
		workingNodes = len(workingPoly)
		for nodeNum in range(workingNodes):
			prevNode = (nodeNum - 1) % workingNodes
			nextNode = (nodeNum + 1) % workingNodes
			
			#Ignore this vertex if it has a reflex angle
			ang = CalcTriangleAng(pts, angleCache, workingPoly[prevNode], workingPoly[nextNode], workingPoly[nodeNum])
			if ang >= math.pi:
				continue
			#print prevNode, nodeNum, nextNode, ang

			#Check if nodes are in this ear
			foundNode = False
			for nodeNum2 in range(workingNodes):
				if nodeNum2 in [prevNode, nodeNum, nextNode]: continue
				chk1 = RightHandedCheck(pts, workingPoly[prevNode], workingPoly[nodeNum], workingPoly[nodeNum2])
				chk2 = RightHandedCheck(pts, workingPoly[nodeNum], workingPoly[nextNode], workingPoly[nodeNum2])
				chk3 = RightHandedCheck(pts, workingPoly[nextNode], workingPoly[prevNode], workingPoly[nodeNum2])
				if chk1 <= 0. and chk2 <= 0. and chk3 <= 0.:
					foundNode = True

				#print "chk", nodeNum2, chk1, chk2, chk3
			if foundNode:
				continue
		
			#print "Found ear at node", nodeNum

			#Store ear
			triangles.append((workingPoly[prevNode], workingPoly[nodeNum], workingPoly[nextNode]))

			#Remove ear from working poly
			workingPoly.pop(nodeNum)

			if 0:
				import matplotlib.pyplot as plt
				import numpy as np
				ptsArr = np.array(pts)
				plt.clf()
				for tri in triangles:
					triTemp = list(tri[:])
					triTemp.append(tri[0])
					print triTemp
					plt.plot(ptsArr[triTemp,0], ptsArr[triTemp,1],'r-')
				plt.plot(ptsArr[workingPoly,0], ptsArr[workingPoly,1],'g-')
				plt.show()

			break

	triangles.append(workingPoly)
	return triangles

if __name__=="__main__":
	import matplotlib.pyplot as plt
	import numpy as np
	#pts = [(2.,1.),(4.,5.),(2.,0.),(0.,5.)]

	#Polygon is defined in an anti-clockwise order
	pts = [(0.,0.),(1.,0.),(1.,1.),(2.,1.),(2.,2.),(4.,2.),(4.,4.),(2.,4.),(2.,3.),(1.,3.),(1.,2.),(0.,2.)]

	triangles = EarClipping(pts)

	#import pyshull
	#triangles = pyshull.FlipTriangles(pts, triangles)

	print triangles

	ptsArr = np.array(pts)
	plt.clf()
	for tri in triangles:
		triTemp = list(tri[:])
		triTemp.append(tri[0])
		plt.plot(ptsArr[triTemp,0], ptsArr[triTemp,1],'r-')
	plt.plot(ptsArr[:,0], ptsArr[:,1],'g-')
	plt.show()



