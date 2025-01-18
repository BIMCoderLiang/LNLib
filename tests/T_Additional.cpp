#include "gtest/gtest.h"
#include "NurbsCurve.h"
#include "NurbsSurface.h"
#include "XYZ.h"
#include "XYZW.h"
#include "MathUtils.h"
#include "LNObject.h"

using namespace LNLib;

TEST(Test_Additional, ApproximateLength)
{
	LN_NurbsCurve result;
	NurbsCurve::CreateLine(XYZ(0, 0, 0), XYZ(100, 0, 0), result);
	EXPECT_TRUE(NurbsCurve::IsLinear(result));
	std::vector<double> lineKnotVector = result.KnotVector;
	double startLineKnot = lineKnotVector[0];
	double endLineKnot = lineKnotVector[lineKnotVector.size() - 1];
	double nonCurvature = NurbsCurve::Curvature(result, (startLineKnot + endLineKnot) / 2.0);
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(nonCurvature, 0.0));
	LN_ArcInfo arcInfo;
	EXPECT_FALSE(NurbsCurve::IsArc(result, arcInfo));
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
	double curvature = NurbsCurve::Curvature(curve, curve.KnotVector[0]);
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(curvature, 1.0/radius));
	EXPECT_TRUE(NurbsCurve::IsArc(curve, arcInfo));
	EXPECT_FALSE(NurbsCurve::IsLinear(curve));
	simpson = NurbsCurve::ApproximateLength(curve, IntegratorType::Simpson);
	gaussLegendre = NurbsCurve::ApproximateLength(curve, IntegratorType::GaussLegendre);
	chebyshev = NurbsCurve::ApproximateLength(curve, IntegratorType::Chebyshev);
	EXPECT_FALSE(MathUtils::IsAlmostEqualTo(simpson, 2 * Constants::Pi * radius)); // not accuracy when use Simpson
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(gaussLegendre, 2 * Constants::Pi * radius));
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(chebyshev, 2 * Constants::Pi * radius));

	simpson = NurbsCurve::GetParamOnCurve(curve, 0.5 * Constants::Pi * radius, IntegratorType::Simpson);
	gaussLegendre = NurbsCurve::GetParamOnCurve(curve, 0.5 * Constants::Pi * radius, IntegratorType::GaussLegendre);
	chebyshev = NurbsCurve::GetParamOnCurve(curve, 0.5 * Constants::Pi * radius, IntegratorType::Chebyshev);
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(simpson, 0.25));
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(gaussLegendre, 0.25));
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(chebyshev, 0.25));

	std::vector<double> arcParameters = NurbsCurve::GetParamsOnCurve(curve, 0.25 * 2 * Constants::Pi * radius, IntegratorType::Simpson);
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(arcParameters[0], 0.25));
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(arcParameters[1], 0.5));
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(arcParameters[2], 0.75));
	arcParameters = NurbsCurve::GetParamsOnCurve(curve, 0.25 * 2 * Constants::Pi * radius, IntegratorType::GaussLegendre);
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(arcParameters[0], 0.25));
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(arcParameters[1], 0.5));
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(arcParameters[2], 0.75));
	arcParameters = NurbsCurve::GetParamsOnCurve(curve, 0.25 * 2 * Constants::Pi * radius, IntegratorType::Chebyshev);
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(arcParameters[0], 0.25));
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(arcParameters[1], 0.5));
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(arcParameters[2], 0.75));
}

TEST(Test_Additional, Normal)
{
	XYZ center = XYZ(0, 0, 0);
	XYZ xAxis = XYZ(1, 0, 0);
	XYZ yAxis = XYZ(0, 1, 0);

	LN_NurbsCurve curve;
	bool result = NurbsCurve::CreateArc(center, xAxis, yAxis, 0, 2 * Constants::Pi, 10, 10, curve);

	XYZ t0 = NurbsCurve::ComputeRationalCurveDerivatives(curve, 1, 0)[1].Normalize();
	XYZ n0 = NurbsCurve::Normal(curve, CurveNormal::Normal, 0.0);
	XYZ b0 = NurbsCurve::Normal(curve, CurveNormal::Binormal, 0.0);
	EXPECT_TRUE(t0.IsAlmostEqualTo(XYZ(0, 1, 0)));
	EXPECT_TRUE(n0.IsAlmostEqualTo(XYZ(-1,0, 0)));
	EXPECT_TRUE(b0.IsAlmostEqualTo(XYZ(0, 0, 1)));

	XYZ t1 = NurbsCurve::ComputeRationalCurveDerivatives(curve, 1, 0.25)[1].Normalize();
	XYZ n1 = NurbsCurve::Normal(curve, CurveNormal::Normal, 0.25);
	XYZ b1 = NurbsCurve::Normal(curve, CurveNormal::Binormal, 0.25);
	EXPECT_TRUE(t1.IsAlmostEqualTo(XYZ(-1,0, 0)));
	EXPECT_TRUE(n1.IsAlmostEqualTo(XYZ(0, -1, 0)));
	EXPECT_TRUE(b1.IsAlmostEqualTo(XYZ(0, 0, 1)));

	XYZ t2 = NurbsCurve::ComputeRationalCurveDerivatives(curve, 1, 0.5)[1].Normalize();
	XYZ n2 = NurbsCurve::Normal(curve, CurveNormal::Normal, 0.5);
	XYZ b2 = NurbsCurve::Normal(curve, CurveNormal::Binormal, 0.5);
	EXPECT_TRUE(t2.IsAlmostEqualTo(XYZ(0, -1, 0)));
	EXPECT_TRUE(n2.IsAlmostEqualTo(XYZ(1, 0,  0)));
	EXPECT_TRUE(b2.IsAlmostEqualTo(XYZ(0, 0, 1)));

	XYZ t3 = NurbsCurve::ComputeRationalCurveDerivatives(curve, 1, 0.75)[1].Normalize();
	XYZ n3 = NurbsCurve::Normal(curve, CurveNormal::Normal, 0.75);
	XYZ b3 = NurbsCurve::Normal(curve, CurveNormal::Binormal, 0.75);
	EXPECT_TRUE(t3.IsAlmostEqualTo(XYZ(1, 0, 0)));
	EXPECT_TRUE(n3.IsAlmostEqualTo(XYZ(0, 1, 0)));
	EXPECT_TRUE(b3.IsAlmostEqualTo(XYZ(0, 0, 1)));

	XYZ t4 = NurbsCurve::ComputeRationalCurveDerivatives(curve, 1, 1.0)[1].Normalize();
	XYZ n4 = NurbsCurve::Normal(curve, CurveNormal::Normal, 1.0);
	XYZ b4 = NurbsCurve::Normal(curve, CurveNormal::Binormal, 1.0);
	EXPECT_TRUE(t4.IsAlmostEqualTo(XYZ(0, 1, 0)));
	EXPECT_TRUE(n4.IsAlmostEqualTo(XYZ(-1, 0, 0)));
	EXPECT_TRUE(b4.IsAlmostEqualTo(XYZ(0, 0, 1)));

}

static void NURBSSurfaceForAreaTest(LN_NurbsSurface& surface, double& standardArea)
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

	surface.DegreeU = degreeU;
	surface.DegreeV = degreeV;
	surface.KnotVectorU = kvU;
	surface.KnotVectorV = kvV;
	surface.ControlPoints = controlPoints;

	standardArea = 4384.255895045;
}

TEST(Test_Additional, Area_Simpson)
{
	LN_NurbsSurface surface;
	double standardArea;
	NURBSSurfaceForAreaTest(surface, standardArea);
	double area = NurbsSurface::ApproximateArea(surface, IntegratorType::Simpson);
	//notice the abs_error
	EXPECT_NEAR(area, standardArea, 1e-4);
}

TEST(Test_Additional, Area_GaussLegendre)
{
	LN_NurbsSurface surface;
	double standardArea;
	NURBSSurfaceForAreaTest(surface, standardArea);
	double area = NurbsSurface::ApproximateArea(surface, IntegratorType::GaussLegendre);
	//notice the abs_error
	EXPECT_NEAR(area, standardArea, 7e-5);
}

TEST(Test_Additional, Area_ChebyShev)
{
	LN_NurbsSurface surface;
	double standardArea;
	NURBSSurfaceForAreaTest(surface, standardArea);
	double area = NurbsSurface::ApproximateArea(surface, IntegratorType::Chebyshev);
	//notice the abs_error
	EXPECT_NEAR(area, standardArea, 2e-3);
}

TEST(Test_Additional, MergeCurve)
{
	// Make line.
	const XYZ start(0, 0, 0);
	const XYZ end(5, 0, 0);
	LN_NurbsCurve line;
	NurbsCurve::CreateLine(start, end, line);

	// Make arc.
	XYZ center(10, 0, 0);
	XYZ xAxis(-1, 0, 0);
	XYZ yAxis(0, 1, 0);
	double startRad = 0;
	double endRad = Constants::Pi;
	double radius = 5;
	LN_NurbsCurve arc;
	NurbsCurve::CreateArc(center, xAxis, yAxis, startRad, endRad, radius, radius, arc);

	// Merge.
	LN_NurbsCurve merged;
	bool success = NurbsCurve::Merge(line, arc, merged);
	EXPECT_TRUE(success);
	
	// 1.0 is the joint knot of merged curve.
	bool canDer = NurbsCurve::CanComputerDerivative(merged, 1.0);
	EXPECT_FALSE(canDer);

	// Test min and max parameters for CanComputerDerivative.
	canDer = NurbsCurve::CanComputerDerivative(merged, merged.KnotVector[0]);
	EXPECT_TRUE(canDer);
	canDer = NurbsCurve::CanComputerDerivative(merged, merged.KnotVector.back());
	EXPECT_TRUE(canDer);

	// Verify length.
	double lineLength = NurbsCurve::ApproximateLength(line);
	double arcLength = NurbsCurve::ApproximateLength(arc);
	double standard = lineLength + arcLength;
	double mergedLength = NurbsCurve::ApproximateLength(merged, IntegratorType::GaussLegendre);
	EXPECT_NEAR(standard, mergedLength, Constants::DistanceEpsilon);

	// The line-arc merged curve is also a good case to test tessellation.
	auto tessPoints = NurbsCurve::Tessellate(merged);
	
	// At least, 5 knot points.
	EXPECT_GT(tessPoints.size(), 5);

	// For each adjacent 2 points, verify angle deflection.
	for(auto i=0;i<tessPoints.size()-1;++i)
	{
		auto& pt1 = tessPoints[i];
		auto& pt2 = tessPoints[i+1];

		// Skip the line-arc joint point.
		if(pt1.IsAlmostEqualTo(end) || pt2.IsAlmostEqualTo(end))
		{
			continue;
		}

		// Get parameters.
		double t1 = NurbsCurve::GetParamOnCurve(merged, pt1);
		double t2 = NurbsCurve::GetParamOnCurve(merged, pt2);
		
		// Compute tangent directions.
		const int derOrder = 1;
		XYZ v1 = NurbsCurve::ComputeRationalCurveDerivatives(merged, derOrder, t1)[1];
		XYZ v2 = NurbsCurve::ComputeRationalCurveDerivatives(merged, derOrder, t2)[1];
		
		// Verify angle.
		double angle = v1.AngleTo(v2);
		EXPECT_LT(std::fabs(angle), Constants::AngleEpsilon);
	}
}


TEST(Test_Fitting, offsetByTillerAndHanson)
{
	// Code copied from issue #25
	std::vector<XYZ> cps_interpolation = { XYZ(XYZ(0, 0, 0)), XYZ(XYZ(1, 1, 0)), XYZ(XYZ(3, 2, 0)), XYZ(XYZ(4, 1, 0)), XYZ(XYZ(5, -1, 0)) };
	LN_NurbsCurve curve, new_curve;
	NurbsCurve::GlobalInterpolation(2, cps_interpolation, curve);
	const double offsetDist = 1.0;
	NurbsCurve::Offset(curve, offsetDist, OffsetType::TillerAndHanson, new_curve);

	// Verify a random proportional point.
	double ratio = 0.345;
	double t = curve.KnotVector[0] * (1 - ratio) + curve.KnotVector[curve.ControlPoints.size()] * ratio;
	auto point = NurbsCurve::GetPointOnCurve(curve, t);
	t = new_curve.KnotVector[0] * (1 - ratio) + new_curve.KnotVector[new_curve.ControlPoints.size()] * ratio;
	auto pointNew = NurbsCurve::GetPointOnCurve(new_curve, t);
	double distance = point.Distance(pointNew);
	double diff = abs(distance - offsetDist);

	double d = diff / offsetDist;
	EXPECT_TRUE(MathUtils::IsLessThanOrEqual(d,0.1));
}

TEST(Test_Fitting, offsetByPieglAndTiller)
{
	// Code copied from issue #25
	std::vector<XYZ> cps_interpolation = { XYZ(XYZ(0, 0, 0)), XYZ(XYZ(1, 1, 0)), XYZ(XYZ(3, 2, 0)), XYZ(XYZ(4, 1, 0)), XYZ(XYZ(5, -1, 0)) };
	LN_NurbsCurve curve, new_curve;
	NurbsCurve::GlobalInterpolation(2, cps_interpolation, curve);
	const double offsetDist = 1.0;
	NurbsCurve::Offset(curve, offsetDist, OffsetType::PieglAndTiller, new_curve);

	// Verify a random proportional point.
	double ratio = 0.345;
	double t = curve.KnotVector[0] * (1 - ratio) + curve.KnotVector[curve.ControlPoints.size()] * ratio;
	auto point = NurbsCurve::GetPointOnCurve(curve, t);
	t = new_curve.KnotVector[0] * (1 - ratio) + new_curve.KnotVector[new_curve.ControlPoints.size()] * ratio;
	auto pointNew = NurbsCurve::GetPointOnCurve(new_curve, t);
	double distance = point.Distance(pointNew);
	double diff = abs(distance - offsetDist);

	double d = diff / offsetDist;
	EXPECT_TRUE(MathUtils::IsLessThanOrEqual(d, 0.1));
}
