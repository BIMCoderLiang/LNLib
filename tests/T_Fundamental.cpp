#include "gtest/gtest.h"
#include "XYZ.h"
#include "XYZW.h"
#include "NurbsCurve.h"
#include "MathUtils.h"
using namespace LNLib;

TEST(Test_Fundamental, All)
{
	{
		int degree = 3;
		std::vector<double> kv = { 0,0,0,0,1,2,3,4,5,5,5,5 };
		double insertKnot = 2.0;

		XYZW P0 = XYZW(XYZ(0, 0, 0), 1);
		XYZW P1 = XYZW(XYZ(1, 5, 0), 1);
		XYZW P2 = XYZW(XYZ(2, 10, 0), 1);
		XYZW P3 = XYZW(XYZ(3, 15, 0), 1);
		XYZW P4 = XYZW(XYZ(4, 20, 0), 1);
		XYZW P5 = XYZW(XYZ(5, 15, 0), 1);
		XYZW P6 = XYZW(XYZ(6, 10, 0), 1);
		XYZW P7 = XYZW(XYZ(7, 5, 0), 1);

		std::vector<XYZW> cp = { P0,P1,P2,P3,P4,P5,P6,P7 };
		std::vector<double> newKv;
		std::vector<XYZW> newCp;
		NurbsCurve::InsertKnot(degree, kv, cp, insertKnot, 1, newKv, newCp);
		EXPECT_TRUE(newKv.size() == kv.size() + 1 &&
			MathUtils::IsAlmostEqualTo(newKv[5], newKv[6]) &&
			MathUtils::IsAlmostEqualTo(newKv[5], 2.0));
		EXPECT_TRUE(newCp[3].IsAlmostEqualTo(2.0 / 3 * P3 + 1.0 / 3 * P2));
		EXPECT_TRUE(newCp[4].IsAlmostEqualTo(1.0 / 3 * P4 + 2.0 / 3 * P3));
		EXPECT_TRUE(newCp[5].IsAlmostEqualTo(P4));
	}
	
	{
		int degree = 2;
		std::vector<double> kv = { 0,0,0,1,2,3,3,3 };
		std::vector<XYZW> cps = { XYZW(XYZ(0,0,0),1), XYZW(XYZ(1,1,0),4), XYZW(XYZ(3,2,0),1), XYZW(XYZ(4,1,0),1), XYZW(XYZ(5,-1,0),1) };
		XYZ result = NurbsCurve::GetPointOnCurveByCornerCut(degree, kv, 1.0, cps);
		EXPECT_TRUE(result.IsAlmostEqualTo(XYZ(7.0 / 5, 6.0 / 5, 0)));
	}
	
}