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