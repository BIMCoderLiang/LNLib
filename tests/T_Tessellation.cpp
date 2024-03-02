#include "gtest/gtest.h"
#include "NurbsCurve.h"
#include "NurbsSurface.h"
#include "XYZ.h"
#include "XYZW.h"
#include "MathUtils.h"
#include "LNObject.h"

using namespace LNLib;

TEST(Test_Tessellation, Curve)
{
	XYZ center = XYZ(0, 0, 0);
	XYZ xAxis = XYZ(1, 0, 0);
	XYZ yAxis = XYZ(0, 1, 0);
	double radius = 100;
	LN_NurbsCurve curve;
	bool createArc = NurbsCurve::CreateArc(center, xAxis, yAxis, 0, 2 * Constants::Pi, radius, radius, curve);
	std::vector<XYZ> points = NurbsCurve::Tessellate(curve);
	EXPECT_TRUE(points.size() > 0);
}