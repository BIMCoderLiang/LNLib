#include "gtest/gtest.h"
#include "NurbsCurve.h"
#include "XYZ.h"
#include "XYZW.h"
#include "MathUtils.h"

using namespace LNLib;

TEST(Test_Advanced, All)
{
	{
		int degree = 2;
		std::vector<double> kv = { 0,0,0,1,2,3,3,3 };
		std::vector<XYZW> cps = { XYZW(XYZ(0,0,0),1), XYZW(XYZ(1,1,0),4), XYZW(XYZ(3,2,0),1), XYZW(XYZ(4,1,0),1), XYZW(XYZ(5,-1,0),1) };
		XYZ result = NurbsCurve::GetPointOnCurve(degree, kv, 1.0, cps);
		double param = NurbsCurve::GetParamOnCurve(degree, kv, cps, result);
		EXPECT_TRUE(MathUtils::IsAlmostEqualTo(1.0,param));
	}
}