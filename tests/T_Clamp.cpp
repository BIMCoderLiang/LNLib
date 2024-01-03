#include "gtest/gtest.h"
#include "NurbsCurve.h"
#include "UV.h"
#include "XYZ.h"
#include "XYZW.h"
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

	XYZ checkP = NurbsCurve::GetPointOnCurve(curve, 0.5);

	LN_NurbsCurve result;
	NurbsCurve::ToUnclampCurve(curve, result);
	curve = result;
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
	XYZ reCheckP = NurbsCurve::GetPointOnCurve(result, 0.5);
	EXPECT_TRUE(checkP.IsAlmostEqualTo(reCheckP));
}
