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
	NurbsCurve::CreateLine(XYZ(0, 0, 0), XYZ(100, 0, 0), result);
	EXPECT_TRUE(NurbsCurve::IsLinear(result));
	EXPECT_FALSE(NurbsCurve::IsArc(result));
	double simpson = NurbsCurve::ApproximateLength(result, IntegratorType::Simpson);
	double gaussLegendre = NurbsCurve::ApproximateLength(result, IntegratorType::GaussLegendre);
	double chebyshev = NurbsCurve::ApproximateLength(result, IntegratorType::Chebyshev);
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(simpson, 100.0));
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(gaussLegendre, 100.0));
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(chebyshev, 100.0));

	XYZ center = XYZ(0, 0, 0);
	XYZ xAxis = XYZ(1, 0, 0);
	XYZ yAxis = XYZ(0, 1, 0);
	LN_NurbsCurve curve;
	bool createArc = NurbsCurve::CreateArc(center, xAxis, yAxis, 0, 2 * Constants::Pi, 100, 100, curve);
	EXPECT_TRUE(createArc);
	EXPECT_TRUE(NurbsCurve::IsArc(curve));
	EXPECT_FALSE(NurbsCurve::IsLinear(curve));
	simpson = NurbsCurve::ApproximateLength(curve, IntegratorType::Simpson);
	gaussLegendre = NurbsCurve::ApproximateLength(curve, IntegratorType::GaussLegendre);
	chebyshev = NurbsCurve::ApproximateLength(curve, IntegratorType::Chebyshev);
	EXPECT_FALSE(MathUtils::IsAlmostEqualTo(simpson, 2 * Constants::Pi * 100)); // not accuracy when use Simpson
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(gaussLegendre, 2 * Constants::Pi * 100));
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(chebyshev, 2 * Constants::Pi * 100));
}