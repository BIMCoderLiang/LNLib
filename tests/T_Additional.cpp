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
	double radius = 100;
	LN_NurbsCurve curve;
	bool createArc = NurbsCurve::CreateArc(center, xAxis, yAxis, 0, 2 * Constants::Pi, radius, radius, curve);
	EXPECT_TRUE(createArc);
	EXPECT_TRUE(NurbsCurve::IsArc(curve));
	EXPECT_FALSE(NurbsCurve::IsLinear(curve));
	simpson = NurbsCurve::ApproximateLength(curve, IntegratorType::Simpson);
	gaussLegendre = NurbsCurve::ApproximateLength(curve, IntegratorType::GaussLegendre);
	chebyshev = NurbsCurve::ApproximateLength(curve, IntegratorType::Chebyshev);
	EXPECT_FALSE(MathUtils::IsAlmostEqualTo(simpson, 2 * Constants::Pi * radius)); // not accuracy when use Simpson
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(gaussLegendre, 2 * Constants::Pi * radius));
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(chebyshev, 2 * Constants::Pi * radius));

	double quadraticParameter = NurbsCurve::GetParamOnCurve(curve, 0.5 * Constants::Pi * radius, IntegratorType::GaussLegendre);
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(quadraticParameter, 0.25));

	std::vector<double> arcParameters = NurbsCurve::GetParamsOnCurve(curve, 0.25 * 2 * Constants::Pi * radius);
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(arcParameters[0], 0.25));
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(arcParameters[1], 0.5));
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(arcParameters[2], 0.75));
}