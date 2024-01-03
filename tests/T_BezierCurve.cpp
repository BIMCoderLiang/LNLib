#include "gtest/gtest.h"
#include "BezierCurve.h"
#include "XYZ.h"
#include "XYZW.h"
#include "MathUtils.h"
#include "LNObject.h"

using namespace LNLib;

TEST(Test_BezierCurve, All)
{
	std::vector<XYZ> controlPoints = { XYZ(1,0,0),XYZ(1,1,0),XYZ(0,2,0) };

	LN_BezierCurve<XYZ> curve;
	curve.Degree = 2;
	curve.ControlPoints = controlPoints;

	XYZ result = BezierCurve::GetPointOnCurveByBernstein(curve, 0.5);
	EXPECT_TRUE(result.IsAlmostEqualTo(XYZ(0.75, 1.0, 0.0)));
	result = BezierCurve::GetPointOnCurveByDeCasteljau(curve, 0.5);
	EXPECT_TRUE(result.IsAlmostEqualTo(XYZ(0.75, 1.0, 0.0)));

	std::vector<XYZW> weightedControlPoints = {XYZW(1,0,0,1),XYZW(1,1,0,1),XYZW(0,2,0,2)};

	LN_BezierCurve<XYZW> curve1;
	curve1.Degree = 2;
	curve1.ControlPoints = weightedControlPoints;

	XYZW weightedResult = BezierCurve::GetPointOnCurveByBernstein(curve1, 0.5);
	EXPECT_TRUE(weightedResult.ToXYZ(true).IsAlmostEqualTo(XYZ(3.0 / 5, 4.0 / 5, 0)));
	weightedResult = BezierCurve::GetPointOnCurveByDeCasteljau(curve1, 0.5);
	EXPECT_TRUE(weightedResult.ToXYZ(true).IsAlmostEqualTo(XYZ(3.0 / 5, 4.0 / 5, 0)));
}