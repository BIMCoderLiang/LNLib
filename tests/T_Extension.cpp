#include "gtest/gtest.h"
#include "XYZ.h"
#include "XYZW.h"
#include "NurbsCurve.h"
using namespace LNLib;

TEST(Test_Extension, Curve)
{
	int degree = 2;
	std::vector<double> kv = { 0,0,0,1,2,3,3,3};
	std::vector<XYZW> cps = { XYZW(XYZ(0,0,0),1), XYZW(XYZ(1,1,0),4), XYZW(XYZ(3,2,0),1), XYZW(XYZ(4,1,0),1), XYZW(XYZ(5,-1,0),1)};
	
	LN_NurbsCurve curve;
	curve.Degree = degree;
	curve.KnotVector = kv;
	curve.ControlPoints = cps;
	
	double originLength = NurbsCurve::ApproximateLength(curve);
	LN_NurbsCurve tangentResult;
	NurbsCurve::Extend(curve, 10, true, ExtensionType::Tangent, tangentResult);
	double extensionLength = NurbsCurve::ApproximateLength(tangentResult);
	EXPECT_TRUE(extensionLength == originLength + 10);
	NurbsCurve::Extend(curve, 10, false, ExtensionType::Tangent, tangentResult);
	extensionLength = NurbsCurve::ApproximateLength(tangentResult);
	EXPECT_TRUE(extensionLength == originLength + 10);
}