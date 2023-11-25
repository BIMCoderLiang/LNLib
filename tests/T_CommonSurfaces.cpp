#include "gtest/gtest.h"
#include "NurbsCurve.h"
#include "NurbsSurface.h"
#include "UV.h"
#include "XYZ.h"
#include "XYZW.h"
#include "MathUtils.h"
#include "LNObject.h"
using namespace LNLib;

TEST(Test_CommonSurfaces, All)
{
	{
		XYZ P00 = XYZ(0, 100, 0);
		XYZ P01 = XYZ(100, 100, 100);
		XYZ P11 = XYZ(100, 0, 0);
		XYZ P10 = XYZ(0, 0, 100);

		LN_Surface surface;
		NurbsSurface::CreateBilinearSurface(P00, P01, P11, P10, surface);
		XYZ result = NurbsSurface::GetPointOnSurface(surface, UV(0.5, 0.5));
		EXPECT_TRUE(result.IsAlmostEqualTo(XYZ(50, 50, 50)));
	}

	{
		XYZ origin = XYZ(0, 0, 0);
		XYZ xAxis = XYZ(1, 0, 0);
		XYZ yAxis = XYZ(0, 1, 0);
		double startRad = 0.0;
		double endRad = 2 * Constants::Pi;
		double radius = 5;
		double height = 5;
		LN_Surface surface;
		bool result = NurbsSurface::CreateCylindricalSurface(origin, xAxis, yAxis, startRad, endRad, radius, height, surface);
		EXPECT_TRUE(result);
		XYZ C0 = NurbsSurface::GetPointOnSurface(surface, UV(0.5, 0.5));
		EXPECT_TRUE(C0.IsAlmostEqualTo(XYZ(-radius,0,radius/2.0)));
		XYZ C1 = NurbsSurface::GetPointOnSurface(surface, UV(0, 0));
		EXPECT_TRUE(C1.IsAlmostEqualTo(XYZ(radius, 0, height)));
		XYZ C2 = NurbsSurface::GetPointOnSurface(surface, UV(1, 0));
		EXPECT_TRUE(C2.IsAlmostEqualTo(XYZ(radius, 0, 0)));
		XYZ C3 = NurbsSurface::GetPointOnSurface(surface, UV(0, 1));
		EXPECT_TRUE(C3.IsAlmostEqualTo(XYZ(radius, 0, height)));
	}

	{
		int degree0 = 2;
		std::vector<double> kv0 = { 0,0,0,1.0 / 4,1.0 / 2,3.0 / 4,1,1,1 };
		std::vector<XYZW> cp0 = { XYZW(20,0,10,1), XYZW(15,0,5,1), XYZW(10,0,5,1),XYZW(5,0,10,1),XYZW(0,0,10,1),XYZW(-5,0,5,1) };
		XYZ C00 = NurbsCurve::GetPointOnCurve(degree0, kv0, 0, cp0);
		XYZ C01 = NurbsCurve::GetPointOnCurve(degree0, kv0, 1, cp0);

		int degree1 = 3;
		std::vector<double> kv1 = { 0,0,0,0,3.0/10,7.0/10,1,1,1,1 };
		std::vector<XYZW> cp1 = { XYZW(20,10,5,1), XYZW(15,10,10,1), XYZW(10,10,10,1),XYZW(5,10,5,1),XYZW(0,10,5,1),XYZW(-5,10,10,1) };
		XYZ C10 = NurbsCurve::GetPointOnCurve(degree1, kv1, 0, cp1);
		XYZ C11 = NurbsCurve::GetPointOnCurve(degree1, kv1, 1, cp1);

		LN_Curve curve0;
		curve0.Degree = degree0;
		curve0.KnotVector = kv0;
		curve0.ControlPoints = cp0;

		LN_Curve curve1;
		curve1.Degree = degree1;
		curve1.KnotVector = kv1;
		curve1.ControlPoints = cp1;

		LN_Surface surface;
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

	{
		int degreeV = 1;
		std::vector<double> kvV = { 0,0,1,1 };
		std::vector<XYZW> cpsV = { XYZW(0,0,1,1), XYZW(1,0,0,1) };
		XYZ origin = XYZ(0, 0, 0);
		XYZ axis = XYZ(0, 0, 1);
		double rad = Constants::Pi / 2.0;

		LN_Curve profile;
		profile.Degree = degreeV;
		profile.KnotVector = kvV;
		profile.ControlPoints = cpsV;

		LN_Surface surface;
		bool result = NurbsSurface::CreateRevolvedSurface(origin, axis, rad, profile, surface);
		EXPECT_TRUE(result);
		XYZ C0 = NurbsSurface::GetPointOnSurface(surface, UV(0.5, 0.5));
		EXPECT_TRUE(C0.IsAlmostEqualTo(XYZ(sqrt(2)/4.0, sqrt(2)/ 4.0, 0.5)));
	}
}