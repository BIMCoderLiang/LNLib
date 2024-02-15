#include "gtest/gtest.h"
#include "NurbsCurve.h"
#include "UV.h"
#include "XYZ.h"
#include "XYZW.h"
#include "KnotVectorUtils.h"
#include "MathUtils.h"
using namespace LNLib;

TEST(Test_Clamp, All)
{
	int degree = 2;
	std::vector<double> kv = { 0,0,0,1.0 / 4,1.0 / 4,1.0 / 2,1.0 / 2,3.0 / 4,3.0 / 4,1,1,1 };
	XYZW P0 = XYZW(10, 0, 0, 1);
	XYZW P1 = XYZW(10, 10, 0, 1);
	XYZW P2 = XYZW(0, 10, 0, 1);
	XYZW P3 = XYZW(-10, 10, 0, 1);
	XYZW P4 = XYZW(-10, 0, 0, 1);
	XYZW P5 = XYZW(-10, -10, 0, 1);
	XYZW P6 = XYZW(0, -10, 0, 1);
	XYZW P7 = XYZW(10, -10, 0, 1);
	XYZW P8 = XYZW(10, 0, 0, 1);
	std::vector<XYZW> cps = { P0,P1,P2,P3,P4,P5,P6,P7,P8 };

	LN_NurbsCurve curve;
	curve.Degree = degree;
	curve.KnotVector = kv;
	curve.ControlPoints = cps;
	EXPECT_TRUE(NurbsCurve::IsClamp(curve));
	EXPECT_TRUE(NurbsCurve::IsClosed(curve));
	XYZ checkP = NurbsCurve::GetPointOnCurve(curve, 0.5);

	LN_NurbsCurve result;
	NurbsCurve::ToUnclampCurve(curve, result);
	EXPECT_TRUE(NurbsCurve::IsClosed(result));
	curve = result;
	EXPECT_FALSE(NurbsCurve::IsClamp(curve));
	XYZ newP = NurbsCurve::GetPointOnCurve(curve, 0.5);
	EXPECT_TRUE(checkP.IsAlmostEqualTo(newP));

	auto kvCopy = curve.KnotVector;
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(kvCopy[0], -1.0 / 4));
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(kvCopy[1], -1.0 / 4));
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(kvCopy[2], 0));
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(kvCopy[3], 1.0 / 4));
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(kvCopy[4], 1.0 / 4));
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(kvCopy[5], 1.0 / 2));
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(kvCopy[6], 1.0 / 2));
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(kvCopy[7], 3.0 / 4));
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(kvCopy[8], 3.0 / 4));
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(kvCopy[9], 1));
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(kvCopy[10], 5.0 / 4));
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(kvCopy[11], 5.0 / 4));

	NurbsCurve::ToClampCurve(curve, result);
	EXPECT_TRUE(NurbsCurve::IsClamp(result));
	EXPECT_TRUE(NurbsCurve::IsClosed(result));
	XYZ reCheckP = NurbsCurve::GetPointOnCurve(result, 0.5);
	EXPECT_TRUE(checkP.IsAlmostEqualTo(reCheckP));

	std::vector<double> k1 = { 0,0,0,0,1,2,3,4,4,4,4 };
	std::vector<double> k2 = { -0.5,-0.5,-0.5,1,2.5,4,4,4 };
	EXPECT_TRUE(KnotVectorUtils::IsUniform(k1));
	EXPECT_TRUE(KnotVectorUtils::IsUniform(k2));
	LN_NurbsCurve curveK1;
	curveK1.Degree = 3;
	curveK1.KnotVector = k1;// ignore control points.
	EXPECT_TRUE(NurbsCurve::IsClamp(curveK1));
	LN_NurbsCurve curveK2;
	curveK2.Degree = 2;
	curveK2.KnotVector = k2;// ignore control points.
	EXPECT_TRUE(NurbsCurve::IsClamp(curveK2));

	std::vector<double> k3 = { 0,0,0,2,3,6,7,7,7 };
	std::vector<double> k4 = { 0,0,0,0,1,2,2,3,4,5,5,5,5 };
	EXPECT_FALSE(KnotVectorUtils::IsUniform(k3));
	EXPECT_FALSE(KnotVectorUtils::IsUniform(k4));
	LN_NurbsCurve curveK3;
	curveK3.Degree = 2;
	curveK3.KnotVector = k3;// ignore control points.
	EXPECT_TRUE(NurbsCurve::IsClamp(curveK3));
	LN_NurbsCurve curveK4;
	curveK4.Degree = 3;
	curveK4.KnotVector = k4;// ignore control points.
	EXPECT_TRUE(NurbsCurve::IsClamp(curveK4));

	std::vector<double> k5 = { -3,-2,-1,0,1,2,3,4,5 };
	std::vector<double> k6 = {0,1,2,3,4,5,6,7,8,9,10,11};
	EXPECT_TRUE(KnotVectorUtils::IsUniform(k5));
	EXPECT_TRUE(KnotVectorUtils::IsUniform(k6));
	LN_NurbsCurve curveK5;
	curveK5.Degree = 2;
	curveK5.KnotVector = k5;// ignore control points.
	EXPECT_FALSE(NurbsCurve::IsClamp(curveK5));
	LN_NurbsCurve curveK6;
	curveK6.Degree = 2;
	curveK6.KnotVector = k6;// ignore control points.
	EXPECT_FALSE(NurbsCurve::IsClamp(curveK6));

	std::vector<double> k7 = { 0,0,1,2,3,4 };
	std::vector<double> k8 = { -2,-1,0,4,5,6,7 };
	EXPECT_FALSE(KnotVectorUtils::IsUniform(k7));
	EXPECT_FALSE(KnotVectorUtils::IsUniform(k8));
	LN_NurbsCurve curveK7;
	curveK7.Degree = 2;
	curveK7.KnotVector = k7;// ignore control points.
	EXPECT_FALSE(NurbsCurve::IsClamp(curveK7));
	LN_NurbsCurve curveK8;
	curveK8.Degree = 2;
	curveK8.KnotVector = k8;// ignore control points.
	EXPECT_FALSE(NurbsCurve::IsClamp(curveK8));
}
