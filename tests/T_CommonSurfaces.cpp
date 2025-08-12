#include "gtest/gtest.h"
#include "NurbsCurve.h"
#include "NurbsSurface.h"
#include "UV.h"
#include "XYZ.h"
#include "XYZW.h"
#include "MathUtils.h"
#include "LNObject.h"
using namespace LNLib;

TEST(Test_CommonSurfaces, CreateBilinearSurface)
{
	//see The NURBS Book 2nd Edition Page335
	XYZ P00 = XYZ(100, 0, 100);
	XYZ P01 = XYZ(0, 0, 0);
	XYZ P10 = XYZ(100, 100, 0);
	XYZ P11 = XYZ(0, 100, 100);

	LN_NurbsSurface surface;
	NurbsSurface::CreateBilinearSurface(P11, P10, P01, P00, surface);
	XYZ mid = NurbsSurface::GetPointOnSurface(surface, UV(0.5, 0.5));
	EXPECT_TRUE(mid.IsAlmostEqualTo(XYZ(50, 50, 50)));

	XYZ uv00 = NurbsSurface::GetPointOnSurface(surface, UV(0, 0));
	XYZ uv01 = NurbsSurface::GetPointOnSurface(surface, UV(0, 1));
	XYZ uv10 = NurbsSurface::GetPointOnSurface(surface, UV(1, 0));
	XYZ uv11 = NurbsSurface::GetPointOnSurface(surface, UV(1, 1));
	EXPECT_TRUE(uv00.IsAlmostEqualTo(P11));
	EXPECT_TRUE(uv01.IsAlmostEqualTo(P10));
	EXPECT_TRUE(uv10.IsAlmostEqualTo(P01));
	EXPECT_TRUE(uv11.IsAlmostEqualTo(P00));
}

TEST(Test_CommonSurfaces, CreateCylindricalSurface)
{
	XYZ origin = XYZ(0, 0, 0);
	XYZ xAxis = XYZ(1, 0, 0);
	XYZ yAxis = XYZ(0, 1, 0);
	double startRad = 0.0;
	double endRad = 2 * Constants::Pi;
	double radius = 5;
	double height = 5;
	LN_NurbsSurface surface;
	bool result = NurbsSurface::CreateCylindricalSurface(origin, xAxis, yAxis, startRad, endRad, radius, height, surface);
	EXPECT_TRUE(result);
	XYZ C0 = NurbsSurface::GetPointOnSurface(surface, UV(0.5, 0.5));
	EXPECT_TRUE(C0.IsAlmostEqualTo(XYZ(-radius, 0, radius / 2.0)));
	XYZ C1 = NurbsSurface::GetPointOnSurface(surface, UV(0, 0));
	EXPECT_TRUE(C1.IsAlmostEqualTo(XYZ(radius, 0, height)));
	XYZ C2 = NurbsSurface::GetPointOnSurface(surface, UV(1, 0));
	EXPECT_TRUE(C2.IsAlmostEqualTo(XYZ(radius, 0, 0)));
	XYZ C3 = NurbsSurface::GetPointOnSurface(surface, UV(0, 1));
	EXPECT_TRUE(C3.IsAlmostEqualTo(XYZ(radius, 0, height)));
}

TEST(Test_CommonSurfaces, CreateRuledSurface)
{
	int degree0 = 2;
	std::vector<double> kv0 = { 0,0,0,1.0 / 4,1.0 / 2,3.0 / 4,1,1,1 };
	std::vector<XYZW> cp0 = { XYZW(20,0,10,1), XYZW(15,0,5,1), XYZW(10,0,5,1),XYZW(5,0,10,1),XYZW(0,0,10,1),XYZW(-5,0,5,1) };

	LN_NurbsCurve curve0;
	curve0.Degree = degree0;
	curve0.KnotVector = kv0;
	curve0.ControlPoints = cp0;

	XYZ C00 = NurbsCurve::GetPointOnCurve(curve0, 0);
	XYZ C01 = NurbsCurve::GetPointOnCurve(curve0, 1);

	int degree1 = 3;
	std::vector<double> kv1 = { 0,0,0,0,3.0 / 10,7.0 / 10,1,1,1,1 };
	std::vector<XYZW> cp1 = { XYZW(20,10,5,1), XYZW(15,10,10,1), XYZW(10,10,10,1),XYZW(5,10,5,1),XYZW(0,10,5,1),XYZW(-5,10,10,1) };

	LN_NurbsCurve curve1;
	curve1.Degree = degree1;
	curve1.KnotVector = kv1;
	curve1.ControlPoints = cp1;

	XYZ C10 = NurbsCurve::GetPointOnCurve(curve1, 0);
	XYZ C11 = NurbsCurve::GetPointOnCurve(curve1, 1);

	LN_NurbsSurface surface;
	NurbsSurface::CreateRuledSurface(curve0, curve1, surface);
	XYZ S00 = NurbsSurface::GetPointOnSurface(surface, UV(0, 0));
	EXPECT_TRUE(S00.IsAlmostEqualTo(C00));
	XYZ S01 = NurbsSurface::GetPointOnSurface(surface, UV(1, 0));
	EXPECT_TRUE(S01.IsAlmostEqualTo(C01));
	XYZ S10 = NurbsSurface::GetPointOnSurface(surface, UV(0, 1));
	EXPECT_TRUE(S10.IsAlmostEqualTo(C10));
	XYZ S11 = NurbsSurface::GetPointOnSurface(surface, UV(1, 1));
	EXPECT_TRUE(S11.IsAlmostEqualTo(C11));
}

TEST(Test_CommonSurfaces, CreateRevolvedSurface)
{
	int degreeV = 1;
	std::vector<double> kvV = { 0,0,1,1 };
	std::vector<XYZW> cpsV = { XYZW(0,0,1,1), XYZW(1,0,0,1) };
	XYZ origin = XYZ(0, 0, 0);
	XYZ axis = XYZ(0, 0, 1);
	double rad = Constants::Pi / 2.0;

	LN_NurbsCurve profile;
	profile.Degree = degreeV;
	profile.KnotVector = kvV;
	profile.ControlPoints = cpsV;

	LN_NurbsSurface surface;
	bool result = NurbsSurface::CreateRevolvedSurface(origin, axis, rad, profile, surface);
	EXPECT_TRUE(result);
	XYZ C0 = NurbsSurface::GetPointOnSurface(surface, UV(0.5, 0.5));
	EXPECT_TRUE(C0.IsAlmostEqualTo(XYZ(sqrt(2) / 4.0, sqrt(2) / 4.0, 0.5)));
}

TEST(Test_CommonSurfaces, RevolveArcToSphere)
{
	// make arc
	LN_NurbsCurve arc;
	XYZ center(0, 0, 0);
	XYZ xAxis(1, 0, 0);
	XYZ yAxis(0, 1, 0);
	double startRad = 0.0;
	double endRad = Constants::Pi;
	double radius = 1.0;
	bool created = NurbsCurve::CreateArc(center, xAxis, yAxis, startRad, endRad, 
		radius, radius, arc);
	EXPECT_TRUE(created);

	// revolve arc to sphere
	XYZ origin(0, 0, 0);
	XYZ axis(1, 0, 0);
	double angle = 2 * Constants::Pi;
	LN_NurbsSurface sphere;
	bool revolved = NurbsSurface::CreateRevolvedSurface(origin, axis, angle, 
		arc, sphere);
	EXPECT_TRUE(revolved);

	// verify area
	double area = NurbsSurface::ApproximateArea(sphere);
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(area, 4 * Constants::Pi * radius * radius));
}
