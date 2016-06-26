from __future__ import print_function
import earclipping, overlap

def CheckResult(expected, actual, description):
	if expected == actual:
		print("Test OK")
	else:
		print("Test Failed:", description)
	return expected == actual

if __name__=="__main__":
	
	#Two horizontal lines
	val = overlap.LineSegmentIntersection(((0.,0.),(10.,0.)),((1.,1.),(9.,1.)))
	CheckResult(False, val, "Two horizontal lines")
	
	#Two vertical lines
	val = overlap.LineSegmentIntersection(((5.,0.),(5.,60.)),((12.,-5.),(12.,50.)))
	CheckResult(False, val, "Two vertical lines")

	#Simple intersection
	val = overlap.LineSegmentIntersection(((0.,0.),(10.,10.)),((0.,10.),(10.,0.)))
	CheckResult(True, val, "Simple intersection")
	
	#Node touching line
	val = overlap.LineSegmentIntersection(((0.,5.),(10.,5.)),((5.,5.),(10.,10.)))
	CheckResult(True, val, "Node touching line")

	#Nearby non-touching line
	val = overlap.LineSegmentIntersection(((0.,0.),(10.,10.)),((2.,7.),(-5.,10.)))
	CheckResult(False, val, "Nearby non-touching line")

	#Horizontal and vertical non-touching lines
	val = overlap.LineSegmentIntersection(((0.,0.),(10.,0.)),((20.,-10.),(20.,10.)))
	CheckResult(False, val, "Horizontal and vertical non-touching lines")

	#Add troublesome real example #1
	commonPt = (243.15634513052646, 200.78910711687058)
	val = overlap.LineSegmentIntersection(((247.27367794327438, 198.35831062297802), commonPt),(commonPt, (243.15231527818833, 200.37556044757366)))
	CheckResult(True, val, "Real line example #1")

