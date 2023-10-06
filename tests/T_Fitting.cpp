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
	}
}

TEST(Test_Fitting, Approximation)
{
	
}