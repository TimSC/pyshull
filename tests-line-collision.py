import earclipping

def CheckResult(expected, actual, description):
	if expected == actual:
		print "Test OK"
	else:
		print "Test Failed:", description
	return expected == actual

if __name__=="__main__":
	
	#Two horizontal lines
	val = earclipping.LineSegmentIntersection(((0.,0.),(10.,0.)),((1.,1.),(9.,1.)))
	CheckResult(False, val, "Two horizontal lines")
	
	#Two vertical lines
	val = earclipping.LineSegmentIntersection(((5.,0.),(5.,60.)),((12.,-5.),(12.,50.)))
	CheckResult(False, val, "Two vertical lines")

	#Simple intersection
	val = earclipping.LineSegmentIntersection(((0.,0.),(10.,10.)),((0.,10.),(10.,0.)))
	CheckResult(True, val, "Simple intersection")
	
	#Node touching line
	val = earclipping.LineSegmentIntersection(((0.,5.),(10.,5.)),((5.,5.),(10.,10.)))
	CheckResult(True, val, "Node touching line")

	#Nearby non-touching line
	val = earclipping.LineSegmentIntersection(((0.,0.),(10.,10.)),((7.,2.),(-5.,10.)))
	CheckResult(False, val, "Nearby non-touching line")

