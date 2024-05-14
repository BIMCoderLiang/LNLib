#include "gtest/gtest.h"
#include "XYZ.h"
#include "XYZW.h"
#include "NurbsCurve.h"
using namespace LNLib;

TEST(Test_NurbsCurve, All)
{
	int degree = 2;
	std::vector<double> kv = { 0,0,0,1,2,3,3,3};
	std::vector<XYZW> cps = { XYZW(XYZ(0,0,0),1), XYZW(XYZ(1,1,0),4), XYZW(XYZ(3,2,0),1), XYZW(XYZ(4,1,0),1), XYZW(XYZ(5,-1,0),1)};
	
	LN_NurbsCurve curve;
	curve.Degree = degree;
	curve.KnotVector = kv;
	curve.ControlPoints = cps;
	
	bool canDer = NurbsCurve::CanComputerDerivative(curve, 1.0);
	EXPECT_TRUE(canDer);

	XYZ result = NurbsCurve::GetPointOnCurve(curve, 1.0);
	EXPECT_TRUE(result.IsAlmostEqualTo(XYZ(7.0/5,6.0/5,0)));

	std::vector<double> kv1 = { 0,0,0,1,1,1 };
	std::vector<XYZW> cps1 = { XYZW(XYZ(1,0,0),1), XYZW(XYZ(1,1,0),1), XYZW(XYZ(0,1,0),2)};

	curve.KnotVector = kv1;
	curve.ControlPoints = cps1;

	std::vector<XYZ> ders = NurbsCurve::ComputeRationalCurveDerivatives(curve, 2, 0.0);
	EXPECT_TRUE(ders[1].IsAlmostEqualTo(XYZ(0, 2, 0)));
	EXPECT_TRUE(ders[2].IsAlmostEqualTo(XYZ(-4, 0, 0)));
}