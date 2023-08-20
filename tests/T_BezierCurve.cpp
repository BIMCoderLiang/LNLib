#include "gtest/gtest.h"
#include "BezierCurve.h"
#include "XYZ.h"
#include "XYZW.h"
#include "MathUtils.h"

using namespace LNLib;

TEST(Test_Bezier, All)
{
	std::vector<XYZ> controlPoints = { XYZ(1,0,0),XYZ(1,1,0),XYZ(0,2,0) };
	XYZ result = BezierCurve::GetPointOnCurveByBernstein(2, controlPoints, 0.5);
	EXPECT_TRUE(result.IsAlmostEqualTo(XYZ(0.75, 1.0, 0.0)));
	result = BezierCurve::GetPointOnCurveByDeCasteljau(2, controlPoints, 0.5);
	EXPECT_TRUE(result.IsAlmostEqualTo(XYZ(0.75, 1.0, 0.0)));

	std::vector<XYZW> weightedControlPoints = {XYZW(1,0,0,1),XYZW(1,1,0,1),XYZW(0,2,0,2)};
	XYZW weightedResult = BezierCurve::GetPointOnRationalCurveByBernstein(2, weightedControlPoints, 0.5);
	EXPECT_TRUE(weightedResult.ToXYZ(true).IsAlmostEqualTo(XYZ(3.0 / 5, 4.0 / 5, 0)));
	weightedResult = BezierCurve::GetPointOnRationalCurveByDeCasteljau(2, weightedControlPoints, 0.5);
	EXPECT_TRUE(weightedResult.ToXYZ(true).IsAlmostEqualTo(XYZ(3.0 / 5, 4.0 / 5, 0)));
}