import earclipping

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
	#print "winding", crossProd
	
	r1 = PointInSideTriangle(tri1[0], tri2, crossProd)
	if r1 is False: return False
	r2 = PointInSideTriangle(tri1[1], tri2, crossProd)
	if r2 is False: return False
	r3 = PointInSideTriangle(tri1[2], tri2, crossProd)
	return r3

def DoTrianglesCollide(tri1, tri2):
	#Do bounding box check
	tri1x = [p[0] for p in tri1]
	tri2x = [p[0] for p in tri2]
	overlapx = earclipping.Check1DOverlap(tri1x, tri2x)
	if not overlapx: return False

	tri1y = [p[1] for p in tri1]
	tri2y = [p[1] for p in tri2]
	overlapy = earclipping.Check1DOverlap(tri1y, tri2y)
	if not overlapy: return False	

	#Check for line overlaps
	for i in range(3):
		for j in range(3):
			crossing = earclipping.LineSegmentIntersection((tri1[i],tri1[(i+1)%3]),(tri2[i],tri2[(j+1)%3]))
			if crossing: return True

	#Check for entirely contained triangle
	contained = CheckFirstTriangleIsContained(tri1, tri2)
	#print "contained1", contained
	if contained: return True

	contained = CheckFirstTriangleIsContained(tri2, tri1)
	#print "contained2", contained
	if contained: return True

	return False

def CheckResult(expected, actual, description):
	if expected == actual:
		print "Test OK", description
	else:
		print "Test Failed:", description
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
