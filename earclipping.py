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

def MergeHoleIntoOuter(workingPoly, pts, outerInd, hole, holeInd):
	#Outer polygon before cut
	filterWorkingPoly = workingPoly[:outerInd+1]
	filteredPts = pts[:]

	#Reorder hole
	reorderedHole = hole[holeInd:]
	reorderedHole.extend(hole[:holeInd])

	#Insert hole
	holdStartInd = len(filteredPts)
	filteredPts.extend(reorderedHole)
	filterWorkingPoly.extend(range(holdStartInd, holdStartInd+len(reorderedHole)))

	#Close hole
	filterWorkingPoly.append(holdStartInd)
	
	#Outer polygon after cut
	filterWorkingPoly.extend(workingPoly[outerInd:])

	return filterWorkingPoly, filteredPts

def line(p1, p2):
	p1 = map(float, p1)
	p2 = map(float, p2)
	A = (p1[1] - p2[1])
	B = (p2[0] - p1[0])
	C = (p1[0]*p2[1] - p2[0]*p1[1])
	return A, B, -C

def InfiniteLineIntersection(L1in, L2in):
	#Based on http://stackoverflow.com/a/20679579
	L1 = line(*L1in)
	L2 = line(*L2in)
	
	D  = L1[0] * L2[1] - L1[1] * L2[0]
	Dx = L1[2] * L2[1] - L1[1] * L2[2]
	Dy = L1[0] * L2[2] - L1[2] * L2[0]
	if D != 0:
		x = Dx / D
		y = Dy / D
		return x,y
	else:
		return False

def IsPointInSegment(L1in, pt):
	#Check if intersection is in L1 segment
	vecL1 = L1in[1][0] - L1in[0][0], L1in[1][1] - L1in[0][1]
	vecAtoI = pt[0] - L1in[0][0], pt[1] - L1in[0][1]
	magVecL1 = (vecL1[0]**2. + vecL1[1]**2.) ** 0.5
	#magVecAtoI = (vecAtoI[0]**2. + vecAtoI[1]**2.) ** 0.5
	if magVecL1 == 0.:
		raise RuntimeError("Zero length input line")
	#if magVecAtoI == 0.:
	#	return True

	vecL1n = (vecL1[0] / magVecL1, vecL1[1] / magVecL1)
	dotProd = vecL1n[0] * vecAtoI[0] + vecL1n[1] * vecAtoI[1]
	if dotProd >= 0. and dotProd < magVecL1:
		return True
	return False

def LineSegmentIntersection(L1in, L2in):

	infIntersect = InfiniteLineIntersection(L1in, L2in)
	if infIntersect is False:
		return False

	chk1 = IsPointInSegment(L1in, infIntersect)
	chk2 = IsPointInSegment(L2in, infIntersect)	

	return chk1 and chk2

def PointVisibility(pts, poly, holeInd, holeNum, holes):
	visiblePoints = []
	#print "holeShape", holeShape
	ptCoord = holes[holeNum][holeInd]

	#Check each point
	for ptIndex, ptNum in enumerate(poly):
		#See if any line segments block
		blocked = False
		for edgeStart, edgeStartPt in enumerate(poly):
			edgeEnd = (edgeStart + 1) % len(poly)
			if poly[edgeStart] == ptNum: continue
			if poly[edgeEnd] == ptNum: continue
			ret = LineSegmentIntersection((ptCoord, pts[ptNum]), (pts[poly[edgeStart]], pts[poly[edgeEnd]]))
			if ret is not False:
				blocked=True
				break

		#Check if the hole self blocks
		holeShape = holes[holeNum]
		for holePtNum, holeChkCoord in enumerate(holeShape):
			if blocked: 
				break
			nextPtNum = (holePtNum + 1) % len(holeShape)
			if holePtNum == holeInd: continue
			if nextPtNum == holeInd: continue
			ret = LineSegmentIntersection((ptCoord, pts[ptNum]), (holeShape[holePtNum], holeShape[nextPtNum]))
			#print ptIndex, holeInd, holePtNum, nextPtNum, ret
			if ret is not False:
				#print (ptCoord, pts[ptNum]), (holeShape[holePtNum], holeShape[nextPtNum])
				blocked=True
		
		#Check if it would be blocked by a future fole
		for holeNumChk, holeShape in enumerate(holes):
			if blocked:
				break
			if holeNumChk == holeNum: continue #Already done self collisions

			for holePtNum, holeChkCoord in enumerate(holeShape):
				if blocked: 
					break
				nextPtNum = (holePtNum + 1) % len(holeShape)
				if holePtNum == holeInd: continue
				if nextPtNum == holeInd: continue
				ret = LineSegmentIntersection((ptCoord, pts[ptNum]), (holeShape[holePtNum], holeShape[nextPtNum]))
				#print ptIndex, holeInd, holePtNum, nextPtNum, ret
				if ret is not False:
					#print (ptCoord, pts[ptNum]), (holeShape[holePtNum], holeShape[nextPtNum])
					blocked=True
			
		#print ptNum, blocked
		if not blocked:
			visiblePoints.append(ptIndex)

	return visiblePoints


def EarClipping(poly, holes):
	#Based on Triangulation by Ear Clipping by David Eberly
	
	angleCache = {}
	triangles = []

	workingPoly = range(len(poly))
	pts = poly[:]

	for holeNum, hole in enumerate(holes):
		#Find place to make cut
		foundCut = None
		for holdPtNum, holeCoord in enumerate(hole):

			visible = PointVisibility(pts, workingPoly, holdPtNum, holeNum, holes)
			#print "vis", holeCoord, visible
			if len(visible) > 0:
				foundCut = (visible[0], holdPtNum)
				break

		if foundCut is None:
			raise RuntimeError("Failed to join hole to other polygon")

		workingPoly, pts = MergeHoleIntoOuter(workingPoly, pts, foundCut[0], hole, foundCut[1])
		#print "wp", workingPoly
		#print "pts", pts

	if 1:
		import matplotlib.pyplot as plt
		import numpy as np
		ptsArr = np.array(pts)
		plt.clf()
		plt.plot(ptsArr[workingPoly,0], ptsArr[workingPoly,1],'g-')
		plt.show()

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
				if workingPoly[nodeNum2] in [workingPoly[prevNode], workingPoly[nodeNum], workingPoly[nextNode]]: continue
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
				print ptsArr[workingPoly,:]
				plt.clf()
				for tri in triangles:
					triTemp = list(tri[:])
					triTemp.append(tri[0])
					#print triTemp
					plt.plot(ptsArr[triTemp,0], ptsArr[triTemp,1],'r-')
				plt.plot(ptsArr[workingPoly,0], ptsArr[workingPoly,1],'g-')
				plt.show()

			break

	triangles.append(workingPoly)
	return pts, triangles

if __name__=="__main__":
	import matplotlib.pyplot as plt
	import numpy as np
	#pts = [(2.,1.),(4.,5.),(2.,0.),(0.,5.)]

	#Polygon is defined in an anti-clockwise order
	outer = [(0.,0.),(1.,0.),(1.,1.),(2.,1.),(2.,2.),(4.,2.),(4.,4.),(2.,4.),(2.,3.),(1.,3.),(1.,2.),(0.,2.)]

	#Holes are specified clockwise
	holes = [[(0.25,0.25),(0.25,0.75),(0.75,0.75),(0.75,0.25)],
		[(3.25,3.25),(3.25,3.75),(3.75,3.75),(3.75,3.25)],
		[(2.1,2.5),(2.5,2.9),(2.9,2.5),(2.5,2.1),]]

	pts, triangles = EarClipping(outer, holes)

	#import pyshull
	#triangles = pyshull.FlipTriangles(pts, triangles)

	print triangles

	ptsArr = np.array(pts)
	outerArr = np.array(outer)
	plt.clf()
	for tri in triangles:
		triTemp = list(tri[:])
		triTemp.append(tri[0])
		plt.plot(ptsArr[triTemp,0], ptsArr[triTemp,1],'r-')
	plt.plot(outerArr[:,0], outerArr[:,1],'g-')
	for hole in holes:
		holeTmp = list(hole[:])
		holeTmp.append(hole[0])
		holeArr = np.array(holeTmp)
		plt.plot(holeArr[:,0], holeArr[:,1],'g-')
	
	plt.show()



