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

TEST(Test_Addintional, Area)
{
	int degreeU = 3;
	int degreeV = 3;
	std::vector<double> kvU = { 0,0,0,0,0.4,0.6,1,1,1,1 };
	std::vector<double> kvV = { 0,0,0,0,0.4,0.6,1,1,1,1 };
	std::vector<std::vector<XYZW>> controlPoints(6, std::vector<XYZW>(6));

	controlPoints[0][0] = XYZW(0, 0, 0, 1);
	controlPoints[0][1] = XYZW(6.666666, 0, 4, 1);
	controlPoints[0][2] = XYZW(16.666666, 0, 22, 1);
	controlPoints[0][3] = XYZW(33.333333, 0, 22, 1);
	controlPoints[0][4] = XYZW(43.333333, 0, 4, 1);
	controlPoints[0][5] = XYZW(50, 0, 0, 1);

	controlPoints[1][0] = XYZW(0, 6.666667, 0, 1);
	controlPoints[1][1] = XYZW(6.6666667, 6.666667, 9.950068, 1);
	controlPoints[1][2] = XYZW(16.6666666, 6.666667, 9.65541838, 1);
	controlPoints[1][3] = XYZW(33.3333333, 6.666667, 47.21371742, 1);
	controlPoints[1][4] = XYZW(43.3333333, 6.666667, -11.56982167, 1);
	controlPoints[1][5] = XYZW(50, 6.6666667, 0, 1);

	controlPoints[2][0] = XYZW(0, 16.666666, 0, 1);
	controlPoints[2][1] = XYZW(6.6666667, 16.666666, -7.9001371, 1);
	controlPoints[2][2] = XYZW(16.6666666, 16.666666, 46.6891632, 1);
	controlPoints[2][3] = XYZW(33.3333333, 16.666667, -28.4274348, 1);
	controlPoints[2][4] = XYZW(43.3333333, 16.666667, 35.1396433, 1);
	controlPoints[2][5] = XYZW(50, 16.6666667, 0, 1);

	controlPoints[3][0] = XYZW(0, 33.3333333, 0, 1);
	controlPoints[3][1] = XYZW(6.6666667, 33.3333333, 29.2877911, 1);
	controlPoints[3][2] = XYZW(16.6666666, 33.3333333, -30.4644718, 1);
	controlPoints[3][3] = XYZW(33.3333333, 33.3333333, 129.1582990, 1);
	controlPoints[3][4] = XYZW(43.3333333, 33.3333333, -62.1717142, 1);
	controlPoints[3][5] = XYZW(50, 33.333333, 0, 1);

	controlPoints[4][0] = XYZW(0, 43.333333, 0, 1);
	controlPoints[4][1] = XYZW(6.6666667, 43.333333, -10.384636, 1);
	controlPoints[4][2] = XYZW(16.6666666, 43.333333, 59.21371742, 1);
	controlPoints[4][3] = XYZW(33.3333333, 43.333333, -37.7272976, 1);
	controlPoints[4][4] = XYZW(43.3333333, 43.333333, 45.1599451, 1);
	controlPoints[4][5] = XYZW(50, 43.333333, 0, 1);

	controlPoints[5][0] = XYZW(0, 50, 0, 1);
	controlPoints[5][1] = XYZW(6.6666667, 50, 0, 1);
	controlPoints[5][2] = XYZW(16.6666666, 50, 0, 1);
	controlPoints[5][3] = XYZW(33.3333333, 50, 0, 1);
	controlPoints[5][4] = XYZW(43.3333333, 50, 0, 1);
	controlPoints[5][5] = XYZW(50, 50, 0, 1);

	LN_NurbsSurface surface;
	surface.DegreeU = degreeU;
	surface.DegreeV = degreeV;
	surface.KnotVectorU = kvU;
	surface.KnotVectorV = kvV;
	surface.ControlPoints = controlPoints;

	double standardArea = 4384.255895045;

	//Need fix, to be continued...
	//simpson = NurbsSurface::ApproximateArea(surface, IntegratorType::Simpson);
	//gaussLegendre = NurbsSurface::ApproximateArea(surface, IntegratorType::GaussLegendre);
	//chebyshev = NurbsSurface::ApproximateArea(surface, IntegratorType::Chebyshev);
}