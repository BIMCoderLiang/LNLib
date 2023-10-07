#include "gtest/gtest.h"
#include "NurbsCurve.h"
#include "NurbsSurface.h"
#include "UV.h"
#include "XYZ.h"
#include "XYZW.h"
#include "MathUtils.h"
#include "Interpolation.h"

using namespace LNLib;

TEST(Test_Fitting, Interpolation)
{
	{
		int degree = 3;
		std::vector<XYZ> Q = { XYZ(0,0,0),XYZ(3,4,0),XYZ(-1,4,0),XYZ(-4,0,0),XYZ(-4,-3,0) };
		auto params = Interpolation::GetChordParameterization(Q);
		auto result = Interpolation::ComputeKnotVector(degree, params);
		EXPECT_TRUE(MathUtils::IsAlmostEqualTo(result[4], 28.0/51));
	}

	{
		int degree = 3;
		std::vector<XYZ> Q = { XYZ(0,0,0),XYZ(3,4,0),XYZ(-1,4,0),XYZ(-4,0,0),XYZ(-4,-3,0) };
		std::vector<double> kv;
		std::vector<XYZW> cps;
		NurbsCurve::GlobalInterpolation(degree, Q, kv, cps);
		EXPECT_TRUE(cps.size() == Q.size());
		EXPECT_TRUE(cps[0].ToXYZ(true).IsAlmostEqualTo(Q[0]));
		EXPECT_TRUE(cps[4].ToXYZ(true).IsAlmostEqualTo(Q[4]));
	}

	{
		int degree = 2;
		std::vector<XYZ> Q = { XYZ(100,0,0),XYZ(0,100,0),XYZ(-100,0,0),XYZ(0,-100,0)};
		std::vector<XYZ> T = { XYZ(0,1,0),XYZ(-1,0,0),XYZ(0,-1,0),XYZ(1,0,0)};
		std::vector<double> kv;
		std::vector<XYZW> cps;
		NurbsCurve::GlobalInterpolation(degree, Q, T, 1, kv, cps);
		XYZ C0 = NurbsCurve::GetPointOnCurve(degree, kv, 0.0, cps);
		EXPECT_TRUE(C0.IsAlmostEqualTo(XYZ(100, 0, 0)));
		XYZ C1 = NurbsCurve::GetPointOnCurve(degree, kv, 1.0, cps);
		EXPECT_TRUE(C1.IsAlmostEqualTo(XYZ(0, -100, 0)));
	}
}

TEST(Test_Fitting, Approximation)
{
	
}