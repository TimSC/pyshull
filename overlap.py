
#Polygon, point and line overlap functions
#by Tim Sheerman-Chase, 2014-2016
from __future__ import print_function

def PointInSideTriangle(pt, tri, winding = None):
	if winding is None:
		winding = GetWindingDirection(tri)

	for side in range(3):
		sideStart = tri[side]
		sideEnd = tri[(side+1)%3]
		sideVec = (sideEnd[0] - sideStart[0], sideEnd[1] - sideStart[1])
		ptVec = (pt[0] - sideStart[0], pt[1] - sideStart[1])

		crossProd = (sideVec[0]*ptVec[1] - ptVec[0]*sideVec[1])

		if winding > 0. and crossProd < 0.: return False
		if winding < 0. and crossProd > 0.: return False

	return True

def Check1DOverlap(range1, range2):
	min1 = min(range1)
	max1 = max(range1)
	min2 = min(range2)
	max2 = max(range2)
	#print(min1, max1, min2, max2)
	if min1 >= min2 and min1 <= max2: return True
	if max1 >= min2 and max1 <= max2: return True
	if min1 <= min2 and max1 >= max2: return True
	return False

def GetWindingDirection(tri2):
	sideVec1 = (tri2[1][0] - tri2[0][0], tri2[1][1] - tri2[0][1])
	sideVec2 = (tri2[2][0] - tri2[1][0], tri2[2][1] - tri2[1][1])
	crossProd = (sideVec1[0]*sideVec2[1] - sideVec2[0]*sideVec1[1])
	if crossProd != 0.: return crossProd #Zero is ambiguous

	sideVec1 = (tri2[2][0] - tri2[1][0], tri2[2][1] - tri2[1][1])
	sideVec2 = (tri2[0][0] - tri2[2][0], tri2[0][1] - tri2[2][1])
	crossProd = (sideVec1[0]*sideVec2[1] - sideVec2[0]*sideVec1[1])
	return crossProd

def CheckFirstTriangleIsContained(tri1, tri2):
	crossProd = GetWindingDirection(tri2)
	#print("winding", crossProd)
	
	r1 = PointInSideTriangle(tri1[0], tri2, crossProd)
	if r1 is False: return False
	r2 = PointInSideTriangle(tri1[1], tri2, crossProd)
	if r2 is False: return False
	r3 = PointInSideTriangle(tri1[2], tri2, crossProd)
	return r3

def CheckResult(expected, actual, description):
	if expected == actual:
		print("Test OK", description)
	else:
		print("Test Failed:", description)
	return expected == actual

def ReorderTriangleThenTest(tri1, tri2, swap, expected, description):
	if not swap:
		result = DoTrianglesCollide(tri1, tri2)
		CheckResult(expected, result, description)
	else:
		result = DoTrianglesCollide(tri1, tri2[::-1])
		CheckResult(expected, result, description)

def RunTriangleTestBattery(tri1, tri2, expected, description):
	ReorderTriangleThenTest(tri1, tri2, 0, expected, description)
	ReorderTriangleThenTest(tri1, tri2, 1, expected, description)
	ReorderTriangleThenTest(tri2, tri1, 0, expected, description)
	ReorderTriangleThenTest(tri2, tri1, 1, expected, description)
	
#*************** Line segment collisions *********************

def line(p1, p2):
	p1 = list(map(float, p1))
	p2 = list(map(float, p2))
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

	xvals = [p[0] for p in L1in]	

	if pt[0] < min(xvals): return False
	if pt[0] > max(xvals): return False
	yvals = [p[1] for p in L1in]
	if pt[1] < min(yvals): return False
	if pt[1] > max(yvals): return False
	return True

#*************** Main public functions ***********************

def DoPointCollideTriangle(pt, tri, winding = None):
	#Do bounding box check
	trix = [p[0] for p in tri]
	if pt[0] < min(trix) or pt[0] > max(trix): return False

	triy = [p[1] for p in tri]
	if pt[1] < min(triy) or pt[1] > max(triy): return False

	#Do full check
	return PointInSideTriangle(pt, tri, winding)

def DoTrianglesCollide(tri1, tri2):
	#Do bounding box check
	tri1x = [p[0] for p in tri1]
	tri2x = [p[0] for p in tri2]
	overlapx = Check1DOverlap(tri1x, tri2x)
	if not overlapx: return False

	tri1y = [p[1] for p in tri1]
	tri2y = [p[1] for p in tri2]
	overlapy = Check1DOverlap(tri1y, tri2y)
	if not overlapy: return False	

	#Check for line overlaps
	for i in range(3):
		for j in range(3):
			crossing = LineSegmentIntersection((tri1[i],tri1[(i+1)%3]),(tri2[i],tri2[(j+1)%3]))
			if crossing: return True

	#Check for entirely contained triangle
	contained = CheckFirstTriangleIsContained(tri1, tri2)
	#print("contained1", contained)
	if contained: return True

	contained = CheckFirstTriangleIsContained(tri2, tri1)
	#print("contained2", contained)
	if contained: return True

	return False

def DoPolyPolyCollision(polyAverts, polyAtris, polyBverts, polyBtris):

	for triA in polyAtris:
		triApos = []
		if len(triA) != 3:
			raise ValueError("Triangle has wrong number of points")
		for ptNum in triA:
			triApos.append(polyAverts[ptNum])

		for triB in polyBtris:
			triBpos = []
			if len(triB) != 3:
				raise ValueError("Triangle has wrong number of points")
			for ptNum in triB:
				triBpos.append(polyBverts[ptNum])

			if DoTrianglesCollide(triApos, triBpos):
				return True

	return False

def LineSegmentIntersection(L1in, L2in):

	#Check if bounding boxes overlap
	L1x = [p[0] for p in L1in]
	L1y = [p[1] for p in L1in]
	L2x = [p[0] for p in L2in]
	L2y = [p[1] for p in L2in]

	#Added end of line comparisons
	#This improves stability but is it the correct approach?
	if L1in[0] == L2in[0]: return True
	if L1in[0] == L2in[1]: return True
	if L1in[1] == L2in[0]: return True
	if L1in[1] == L2in[1]: return True

	if Check1DOverlap(L1x, L2x) is False:
		#print("fail x")
		return False
	if Check1DOverlap(L1y, L2y) is False:
		#print("fail y")
		return False

	#Find intersection assuming lines are infinitely long
	infIntersect = InfiniteLineIntersection(L1in, L2in)
	if infIntersect is False:
		return False
	chk1 = IsPointInSegment(L1in, infIntersect)

	if chk1 is False: return False
	chk2 = IsPointInSegment(L2in, infIntersect)	

	return chk1 and chk2

# **************** Test functions *******************

if __name__ == "__main__":
	#Unit tests

	#Identical triangle
	RunTriangleTestBattery(((0.,0.),(1.,0.),(0.5,1.)),((0.,0.),(1.,0.),(0.5,1.)), 
		True, "Identical triangle")

	#Overlapping triangles, crossing edge, point on edge
	RunTriangleTestBattery(((0.5,0.),(1.5,0.),(1.5,1.)),((0.,0.),(1.,0.),(0.5,1.)), 
		True, "Overlapping triangles, crossing edge, point on edge")

	#Overlapping triangles, crossing edge, point is interior
	RunTriangleTestBattery(((0.0,0.5),(1.,0.),(1.,1.)),((0.5,0.5),(1.5,0.),(1.5,1.)),
		True, "Overlapping triangles, crossing edge, point is interior")

	#Non overlapping triangles
	RunTriangleTestBattery(((0.,0.),(1.,0.),(0.5,1.)),((10.,0.),(11.,0.),(10.5,1.)),
		False, "Non overlapping triangles")

	#Nearby, non overlapping triangles
	RunTriangleTestBattery(((0.,0.),(1.,0.),(0.5,1.)),((0.6,1.),(1.1,0.),(1.6,1.)),
		False, "Nearby, non overlapping triangles")

	#Common point, no crossing
	RunTriangleTestBattery(((0.,0.),(1.,0.),(0.5,1.)),((0.,10.0),(1.,10.),(0.5,1.)),
		True, "Common point, no crossing")

	#Shared edge, no crossing
	RunTriangleTestBattery(((0.,0.),(1.,0.),(0.5,1.)),((0.,0.),(1.,0.),(0.5,-1.)),
		True, "Common point, no crossing")

	#Contained triangle
	RunTriangleTestBattery(((0.,0.),(10.,0.),(5.,10.)),((4.,1.),(6.,1.),(5.,2.)),
		True, "Contained triangle")



	#result = DoTrianglesCollide(((),(),()),((),(),()))

	#result = DoTrianglesCollide(((),(),()),((),(),()))
