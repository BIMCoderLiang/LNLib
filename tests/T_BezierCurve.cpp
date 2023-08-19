#include "gtest/gtest.h"
#include "BezierCurve.h"
#include "XYZ.h"
#include "XYZW.h"
#include "MathUtils.h"

using namespace LNLib;

TEST(Test_Bezier, All)
{
	std::vector<XYZW> controlPoints = {XYZW(1,0,0,1),XYZW(1,1,0,1),XYZW(0,2,0,2)};
	XYZ result = BezierCurve::GetPointOnRationalCurveByBernstein(2, controlPoints, 0.5);
	EXPECT_TRUE(result.IsAlmostEqualTo(XYZ(3.0 / 5, 4.0 / 5, 0)));
	result = BezierCurve::GetPointOnRationalCurveByDeCasteljau(2, controlPoints, 0.5);
	EXPECT_TRUE(result.IsAlmostEqualTo(XYZ(3.0 / 5, 4.0 / 5, 0)));
}