#include "gtest/gtest.h"
#include "XYZ.h"
#include "XYZW.h"
#include "NurbsCurve.h"
#include "NurbsSurface.h"
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

	{
		int degree = 3;
		std::vector<double> kv = { 0,0,0,0,1,2,3,4,5,5,5,5 };
		double insertKnot = 2.0;

		XYZW P00 = XYZW(XYZ(0, 0, 0), 1);
		XYZW P10 = XYZW(XYZ(1, 5, 0), 1);
		XYZW P20 = XYZW(XYZ(2, 10, 0), 1);
		XYZW P30 = XYZW(XYZ(3, 15, 0), 1);
		XYZW P40 = XYZW(XYZ(4, 20, 0), 1);
		XYZW P50 = XYZW(XYZ(5, 15, 0), 1);
		XYZW P60 = XYZW(XYZ(6, 10, 0), 1);
		XYZW P70 = XYZW(XYZ(7, 5, 0), 1);

		XYZW P01 = XYZW(XYZ(10, 0, 0), 1);
		XYZW P11 = XYZW(XYZ(11, 5, 0), 1);
		XYZW P21 = XYZW(XYZ(12, 10, 0), 1);
		XYZW P31 = XYZW(XYZ(13, 15, 0), 1);
		XYZW P41 = XYZW(XYZ(14, 20, 0), 1);
		XYZW P51 = XYZW(XYZ(15, 15, 0), 1);
		XYZW P61 = XYZW(XYZ(16, 10, 0), 1);
		XYZW P71 = XYZW(XYZ(17, 5, 0), 1);

		std::vector<std::vector<XYZW>> cp = { 
			
			{P00,P01},
			{P10,P11},
			{P20,P21},
			{P30,P31},
			{P40,P41},
			{P50,P51},
			{P60,P61},
			{P70,P71},
		
		};
		std::vector<double> newKv;
		std::vector<std::vector<XYZW>> newCp;
		NurbsSurface::InsertKnot(degree, kv, cp, insertKnot, 1, true, newKv, newCp);
		EXPECT_TRUE(newKv.size() == kv.size() + 1 &&
					MathUtils::IsAlmostEqualTo(newKv[5], newKv[6]) &&
					MathUtils::IsAlmostEqualTo(newKv[5], 2.0));
		EXPECT_TRUE(newCp[3][0].IsAlmostEqualTo(2.0 / 3 * P30 + 1.0 / 3 * P20));
		EXPECT_TRUE(newCp[4][0].IsAlmostEqualTo(1.0 / 3 * P40 + 2.0 / 3 * P30));
		EXPECT_TRUE(newCp[5][0].IsAlmostEqualTo(P40));
	}

	{
		int degree = 3;
		std::vector<double> kv = { 0,0,0,0,1,2,3,4,5,5,5,5 };
		double insertKnot = 2.0;

		XYZW P00 = XYZW(XYZ(0, 0, 0), 1);
		XYZW P01 = XYZW(XYZ(1, 5, 0), 1);
		XYZW P02 = XYZW(XYZ(2, 10, 0), 1);
		XYZW P03 = XYZW(XYZ(3, 15, 0), 1);
		XYZW P04 = XYZW(XYZ(4, 20, 0), 1);
		XYZW P05 = XYZW(XYZ(5, 15, 0), 1);
		XYZW P06 = XYZW(XYZ(6, 10, 0), 1);
		XYZW P07 = XYZW(XYZ(7, 5, 0), 1);

		XYZW P10 = XYZW(XYZ(10, 0, 0), 1);
		XYZW P11 = XYZW(XYZ(11, 5, 0), 1);
		XYZW P12 = XYZW(XYZ(12, 10, 0), 1);
		XYZW P13 = XYZW(XYZ(13, 15, 0), 1);
		XYZW P14 = XYZW(XYZ(14, 20, 0), 1);
		XYZW P15 = XYZW(XYZ(15, 15, 0), 1);
		XYZW P16 = XYZW(XYZ(16, 10, 0), 1);
		XYZW P17 = XYZW(XYZ(17, 5, 0), 1);

		std::vector<std::vector<XYZW>> cp = {

			{P00,P01,P02,P03,P04,P05,P06,P07},
			{P10,P11,P12,P13,P14,P15,P16,P17},

		};

		std::vector<double> newKv;
		std::vector<std::vector<XYZW>> newCp;
		NurbsSurface::InsertKnot(degree, kv, cp, insertKnot, 1, false, newKv, newCp);
		EXPECT_TRUE(newKv.size() == kv.size() + 1 &&
					MathUtils::IsAlmostEqualTo(newKv[5], newKv[6]) &&
					MathUtils::IsAlmostEqualTo(newKv[5], 2.0));
		EXPECT_TRUE(newCp[0][3].IsAlmostEqualTo(2.0 / 3 * P03 + 1.0 / 3 * P02));
		EXPECT_TRUE(newCp[0][4].IsAlmostEqualTo(1.0 / 3 * P04 + 2.0 / 3 * P03));
		EXPECT_TRUE(newCp[0][5].IsAlmostEqualTo(P04));
	}
}