#include "gtest/gtest.h"
#include "NurbsCurve.h"
#include "NurbsSurface.h"
#include "UV.h"
#include "XYZ.h"
#include "XYZW.h"
#include "MathUtils.h"
using namespace LNLib;

TEST(Test_CommonSurfaces, All)
{
	{
		XYZ P00 = XYZ(0, 100, 0);
		XYZ P01 = XYZ(100, 100, 100);
		XYZ P11 = XYZ(100, 0, 0);
		XYZ P10 = XYZ(0, 0, 100);

		int degreeU, degreeV;
		std::vector<double> kvU, kvV;
		std::vector<std::vector<XYZW>> cps;
		NurbsSurface::CreateBilinearSurface(P00, P01, P11, P10, degreeU, degreeV, kvU, kvV, cps);
		XYZ result = NurbsSurface::GetPointOnSurface(degreeU, degreeV, kvU, kvV, UV(0.5, 0.5), cps);
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
		int degreeU, degreeV;
		std::vector<double> kvU, kvV;
		std::vector<std::vector<XYZW>> cps;
		bool result = NurbsSurface::CreateCylindricalSurface(origin, xAxis, yAxis, startRad, endRad, radius, height, degreeU, degreeV, kvU, kvV, cps);
		EXPECT_TRUE(result);
		XYZ C0 = NurbsSurface::GetPointOnSurface(degreeU, degreeV, kvU, kvV, UV(0.5, 0.5), cps);
		EXPECT_TRUE(C0.IsAlmostEqualTo(XYZ(-radius,0,radius/2.0)));
		XYZ C1 = NurbsSurface::GetPointOnSurface(degreeU, degreeV, kvU, kvV, UV(0, 0), cps);
		EXPECT_TRUE(C1.IsAlmostEqualTo(XYZ(radius, 0, height)));
		XYZ C2 = NurbsSurface::GetPointOnSurface(degreeU, degreeV, kvU, kvV, UV(1, 0), cps);
		EXPECT_TRUE(C2.IsAlmostEqualTo(XYZ(radius, 0, 0)));
		XYZ C3 = NurbsSurface::GetPointOnSurface(degreeU, degreeV, kvU, kvV, UV(0, 1), cps);
		EXPECT_TRUE(C3.IsAlmostEqualTo(XYZ(radius, 0, height)));
	}
}