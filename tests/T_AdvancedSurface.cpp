#include "gtest/gtest.h"
#include "NurbsCurve.h"
#include "NurbsSurface.h"
#include "KnotVectorUtils.h"
#include "UV.h"
#include "XYZ.h"
#include "XYZW.h"
#include "MathUtils.h"

using namespace LNLib;

TEST(Test_AdvancedSurface, KnotVector)
{
	std::vector<double> kv0 = { 0,0,0,1.0 / 4,1.0 / 2,3.0 / 4,1,1,1 };
	std::vector<double> kv1 = { 0,0,0,0,3.0 / 10,7.0 / 10,1,1,1,1 };

	std::vector<std::vector<double>> kvs;
	kvs.emplace_back(kv0);
	kvs.emplace_back(kv1);
	auto result = KnotVectorUtils::GetInsertedKnotElements(kvs);
	EXPECT_TRUE(result[0].size() == 4);
	EXPECT_TRUE(result[1].size() == 3);
}

TEST(Test_AdvancedSurface, CoonsSurface)
{
	XYZ l1 = XYZ(0, 100, 10);
	XYZ l2 = XYZ(0, 70, 20);
	XYZ l3 = XYZ(0, 40, -10);
	XYZ l4 = XYZ(0, 0, 0);

	auto leftCps = {XYZW(l1,1),XYZW(l2,1),XYZW(l3,1),XYZW(l4,1)};

	LN_NurbsCurve leftCurve;
	leftCurve.Degree = 3;
	leftCurve.ControlPoints = leftCps;
	leftCurve.KnotVector = { 0,0,0,0,1,1,1,1 };

	XYZ b1 = XYZ(0, 0, 0);
	XYZ b2 = XYZ(20, 0, 10);
	XYZ b3 = XYZ(40, 0, 30);
	XYZ b4 = XYZ(60, 0, 50);

	auto bottomCps = { XYZW(b1,1),XYZW(b2,1),XYZW(b3,1),XYZW(b4,1) };

	LN_NurbsCurve bottomCurve;
	bottomCurve.Degree = 3;
	bottomCurve.ControlPoints = bottomCps;
	bottomCurve.KnotVector = { 0,0,0,0,1,1,1,1 };

	XYZ r1 = XYZ(60, 0,  50);
	XYZ r2 = XYZ(60, 40, 10);
	XYZ r3 = XYZ(60, 70, 50);
	XYZ r4 = XYZ(60, 100,30);

	auto rightCps = { XYZW(r1,1),XYZW(r2,1),XYZW(r3,1),XYZW(r4,1) };

	LN_NurbsCurve rightCurve;
	rightCurve.Degree = 3;
	rightCurve.ControlPoints = rightCps;
	rightCurve.KnotVector = { 0,0,0,0,1,1,1,1 };

	XYZ t1 = XYZ(60, 100, 30);
	XYZ t2 = XYZ(40, 100, 10);
	XYZ t3 = XYZ(20, 100, 20);
	XYZ t4 = XYZ(0, 100, 10);

	auto topCps = { XYZW(t1,1),XYZW(t2,1),XYZW(t3,1),XYZW(t4,1) };

	LN_NurbsCurve topCurve;
	topCurve.Degree = 3;
	topCurve.ControlPoints = topCps;
	topCurve.KnotVector = { 0,0,0,0,1,1,1,1 };

	LN_NurbsSurface coonsSurface;
	NurbsSurface::CreateCoonsSurface(leftCurve, bottomCurve, rightCurve, topCurve, coonsSurface);

	XYZ check0 = NurbsSurface::GetPointOnSurface(coonsSurface, UV(0, 0));
	XYZ check1 = NurbsSurface::GetPointOnSurface(coonsSurface, UV(1, 0));
	XYZ check2 = NurbsSurface::GetPointOnSurface(coonsSurface, UV(0, 1));
	XYZ check3 = NurbsSurface::GetPointOnSurface(coonsSurface, UV(1, 1));

	EXPECT_TRUE(check0.IsAlmostEqualTo(b1));
	EXPECT_TRUE(check1.IsAlmostEqualTo(l1));
	EXPECT_TRUE(check2.IsAlmostEqualTo(r1));
	EXPECT_TRUE(check3.IsAlmostEqualTo(t1));
}

TEST(Test_AdvancedSurface, ProjectNormal)
{
	LN_NurbsCurve circle0;
	XYZ center = XYZ(0, 0, 0);
	XYZ xAxis = XYZ(1, 0, 0);
	XYZ yAxis = XYZ(0, 1, 0);
	NurbsCurve::CreateArc(center, xAxis, yAxis, 0, Constants::Pi, 50, 50, circle0);

	auto blist = NurbsCurve::ProjectNormal(circle0);
	std::vector<double> kv = circle0.KnotVector;
	for (int i = 0; i < kv.size(); i++)
	{
		XYZ ti = NurbsCurve::ComputeRationalCurveDerivatives(circle0, 1, kv[i])[1].Normalize();
		XYZ bi = blist[i];

		EXPECT_TRUE(MathUtils::IsAlmostEqualTo(ti.DotProduct(bi), 0));
	}
}

TEST(Test_AdvancedSurface, LoftSurface)
{
	//arc0
	LN_NurbsCurve circle0;
	XYZ center = XYZ(0, 0, 0);
	XYZ xAxis = XYZ(1, 0, 0);
	XYZ yAxis = XYZ(0, 1, 0);
	NurbsCurve::CreateArc(center, xAxis, yAxis, 0, Constants::Pi, 50, 50, circle0);
	//arc1
	LN_NurbsCurve circle1;
	XYZ center1 = XYZ(0, 0, 50);
	XYZ xAxis1 = XYZ(1, 0, 0);
	XYZ yAxis1 = XYZ(0, 1, 0);
	NurbsCurve::CreateArc(center1, xAxis1, yAxis1, 0, Constants::Pi, 50, 50, circle1);
	//arc2
	LN_NurbsCurve circle2;
	XYZ center2 = XYZ(0, 0, 100);
	XYZ xAxis2 = XYZ(1, 0, 0);
	XYZ yAxis2 = XYZ(0, 1, 0);
	NurbsCurve::CreateArc(center2, xAxis2, yAxis2, 0, Constants::Pi, 50, 50, circle2);
	std::vector<LN_NurbsCurve> arc_profiles;
	arc_profiles.push_back(circle0);
	arc_profiles.push_back(circle1);
	arc_profiles.push_back(circle2);
	LN_NurbsSurface arc_surface;
	NurbsSurface::CreateLoftSurface(arc_profiles, arc_surface);
	auto knotVectorU = arc_surface.KnotVectorU;
	auto knotVectorV = arc_surface.KnotVectorV;

	double uStart = knotVectorU[0];
	double uEnd = knotVectorU[knotVectorU.size() - 1];
	double uHalf = (uStart + uEnd) / 2.0;
	double vStart = knotVectorV[0];
	double vEnd = knotVectorV[knotVectorV.size() - 1];
	double vHalf = (vStart + vEnd) / 2.0;

	XYZ point1 = NurbsSurface::GetPointOnSurface(arc_surface, UV(uStart, vStart));
	EXPECT_TRUE(point1.IsAlmostEqualTo(XYZ(50, 0, 0)));
	XYZ point2 = NurbsSurface::GetPointOnSurface(arc_surface, UV(uStart, vEnd));
	EXPECT_TRUE(point2.IsAlmostEqualTo(XYZ(50, 0, 100)));
	XYZ point3 = NurbsSurface::GetPointOnSurface(arc_surface, UV(uEnd, vStart));
	EXPECT_TRUE(point3.IsAlmostEqualTo(XYZ(-50, 0, 0)));
	XYZ point4 = NurbsSurface::GetPointOnSurface(arc_surface, UV(uEnd, vEnd));
	EXPECT_TRUE(point4.IsAlmostEqualTo(XYZ(-50, 0, 100)));
	XYZ point5 = NurbsSurface::GetPointOnSurface(arc_surface, UV(uHalf, vHalf));
	EXPECT_TRUE(point5.IsAlmostEqualTo(XYZ(0, 50, 50)));

}

TEST(Test_AdvancedSurface, CreateGeneralizedTranslationalSweepSurface)
{
	// Make circular profile.
	XYZ center(0, 0, 0);
	XYZ xAxis(1, 0, 0);
	XYZ yAxis(0, 1, 0);
	double startRad = 0;
	double endRad = Constants::Pi * 2;
	double radius = 5;
	LN_NurbsCurve profile;
	NurbsCurve::CreateArc(center, xAxis, yAxis, startRad, endRad, radius, radius, profile);

	// Make linear trajectory.
	LN_NurbsCurve trajectory;
	XYZ start(0, 0, 0);
	double height = 7;
	XYZ end(0, 0, height);
	NurbsCurve::CreateLine(start, end, trajectory);

	// Make translational swept surface.
	LN_NurbsSurface surface;
	NurbsSurface::CreateGeneralizedTranslationalSweepSurface(profile, trajectory, surface);

	// Verify area.
	double expectedArea = (endRad - startRad) * radius * height;
	double area = NurbsSurface::ApproximateArea(surface);
	EXPECT_NEAR(area, expectedArea, 1e-4);


	NurbsSurface::CreateSweepSurface(profile, trajectory, 5, surface);
	area = NurbsSurface::ApproximateArea(surface);
	EXPECT_NEAR(area, expectedArea, 1e-4);
}
