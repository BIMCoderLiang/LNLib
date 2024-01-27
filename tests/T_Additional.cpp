#include "gtest/gtest.h"
#include "NurbsCurve.h"
#include "NurbsSurface.h"
#include "XYZ.h"
#include "XYZW.h"
#include "MathUtils.h"
#include "LNObject.h"

using namespace LNLib;

TEST(Test_Addintional, All)
{
	LN_NurbsCurve result;
	NurbsCurve::CreateLine(XYZ(0, 0, 0), XYZ(1, 0, 0), result);
	EXPECT_TRUE(NurbsCurve::IsLinear(result));
	EXPECT_FALSE(NurbsCurve::IsArc(result));

	XYZ center = XYZ(0, 0, 0);
	XYZ xAxis = XYZ(1, 0, 0);
	XYZ yAxis = XYZ(0, 1, 0);
	LN_NurbsCurve curve;
	bool createArc = NurbsCurve::CreateArc(center, xAxis, yAxis, 0, 2 * Constants::Pi, 10, 10, curve);
	EXPECT_TRUE(createArc);
	EXPECT_TRUE(NurbsCurve::IsArc(curve));
	EXPECT_FALSE(NurbsCurve::IsLinear(curve));
}