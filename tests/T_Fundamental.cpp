#include "gtest/gtest.h"
#include "XYZ.h"
#include "XYZW.h"
#include "NurbsCurve.h"
#include "NurbsSurface.h"
#include "ValidationUtils.h"
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

	{
		int degree = 4;
		std::vector<double> kv = { 0,0,0,0,0,0.2,0.4,0.6,0.8,1,1,1,1,1 };
		std::vector<double> ike = { 0.1,0.1,0.1,0.1,0.2,0.2,0.2,0.3,0.3,0.3,0.3,0.4,0.4,0.4,0.5,0.5,0.5,0.5,
									0.6,0.6,0.6,0.7,0.7,0.7,0.7,0.8,0.8,0.8,0.9,0.9,0.9,0.9};

		XYZW P0 = XYZW(XYZ(5, 10, 0), 1);
		XYZW P1 = XYZW(XYZ(15, 25, 0), 1);
		XYZW P2 = XYZW(XYZ(30, 30, 0), 1);
		XYZW P3 = XYZW(XYZ(45, 5, 0), 1);
		XYZW P4 = XYZW(XYZ(55, 5, 0), 1);
		XYZW P5 = XYZW(XYZ(70, 40, 0), 1);
		XYZW P6 = XYZW(XYZ(60, 60, 0), 1);
		XYZW P7 = XYZW(XYZ(35, 60, 0), 1);
		XYZW P8 = XYZW(XYZ(20, 40, 0), 1);

		std::vector<XYZW> cp = { P0,P1,P2,P3,P4,P5,P6,P7,P8 };
		std::vector<double> newKv;
		std::vector<XYZW> newCp;
		NurbsCurve::RefineKnotVector(degree, kv, cp, ike, newKv, newCp);
		EXPECT_TRUE(newKv.size() == 46 && newCp.size() == 41);
		EXPECT_TRUE(newCp[20].ToXYZ(true).IsAlmostEqualTo(XYZ(55.9157986, 12.17447916,0)));
	}

	{
		int degreeU = 3;
		int degreeV = 3;
		std::vector<double> kvU = { 0,0,0,0,1,2,3,3,3,3 };
		std::vector<double> kvV = { 0,0,0,0,1,2,3,3,3,3 };

		std::vector < std::vector<XYZW>> cps = {

			{XYZW(25,-25,0,1),XYZW(15,-25,0,1),XYZW(5,-25,0,1),XYZW(-5,-25,0,1),XYZW(-15,-25,0,1),XYZW(-25,-25,0,1)},
			{XYZW(25,-15,0,1),XYZW(15,-15,0,1),XYZW(5,-15,0,1),XYZW(-5,-15,0,1),XYZW(-15,-15,0,1),XYZW(-25,-15,0,1)},
			{XYZW(25,-5,5,1),XYZW(15,-5,5,1),XYZW(5,-5,5,1),XYZW(-5,-5,5,1),XYZW(-15,-5,5,1),XYZW(-25,-5,5,1)},
			{XYZW(25,5,5,1),XYZW(15,5,5,1),XYZW(5,5,5,1),XYZW(-5,5,5,1),XYZW(-15,5,5,1),XYZW(-25,5,5,1)},
			{XYZW(25,15,0,1),XYZW(15,15,0,1),XYZW(5,15,5,1),XYZW(-5,15,5,1),XYZW(-15,15,0,1),XYZW(-25,15,0,1)},
			{XYZW(25,25,0,1),XYZW(15,25,0,1),XYZW(5,25,5,1),XYZW(-5,25,5,1),XYZW(-15,25,0,1),XYZW(-25,25,0,1)},
			
		};
		std::vector<double> newKu;
		std::vector<double> newKv;
		std::vector < std::vector<XYZW>> newCps;
		std::vector<double> ike = { 0.5,0.5,0.5,1.0,1.0,1.5,1.5,1.5,2,2,2.5,2.5,2.5 };
		NurbsSurface::RefineKnotVector(degreeU, degreeV, kvU, kvV, cps, ike, true, newKu, newKv, newCps);
		EXPECT_TRUE(newKu.size() == 23 && newKv.size() == kvV.size());
		EXPECT_TRUE(newCps[4][0].ToXYZ(true).IsAlmostEqualTo(XYZ(25, -10.208333, 2.1875)));
	}

	{
		int degree = 3;
		std::vector<double> kv = { 0,0,0,0,0.33,0.66,1.0,1.0,1.0,1.0 };
		std::vector<XYZW> cps = {
			XYZW(5,5,0,1),XYZW(10,10,0,1),XYZW(20,15,0,1),XYZW(35,15,0,1),XYZW(45,10,0,1),XYZW(50,5,0,1)
		};
		std::vector<double> newKv;
		std::vector<XYZW> newCps;
		NurbsCurve::RemoveKnot(degree, kv, cps, 0.66, 1, newKv, newCps);
		EXPECT_TRUE(newKv.size() == kv.size()-1);
		EXPECT_TRUE(MathUtils::IsAlmostEqualTo(newKv[4], 0.33));
		EXPECT_TRUE(MathUtils::IsAlmostEqualTo(newKv[5], 1.0));
		EXPECT_TRUE(newCps.size() == 5);
		EXPECT_TRUE(newCps[0].IsAlmostEqualTo(XYZW(5, 5, 0, 1)));
		EXPECT_TRUE(newCps[1].IsAlmostEqualTo(XYZW(10, 10, 0, 1)));
		EXPECT_TRUE(newCps[2].IsAlmostEqualTo(XYZW(25.15151515, 17.5757575757, 0, 1)));
		EXPECT_TRUE(newCps[3].IsAlmostEqualTo(XYZW(45, 10, 0, 1)));
		EXPECT_TRUE(newCps[4].IsAlmostEqualTo(XYZW(50, 5, 0, 1)));
	}

	{
		int degree = 3;
		std::vector<double> kvU = { 0,0,0,0,0.33,0.66,1.0,1.0,1.0,1.0 };
		std::vector<double> kvV = { 0,0,0,0,0.33,0.66,1.0,1.0,1.0,1.0 };
		std::vector<std::vector<XYZW>> cps = 
		{
			{XYZW(5,5,0,1),XYZW(10,10,0,1),XYZW(20,15,0,1),XYZW(35,15,0,1),XYZW(45,10,0,1),XYZW(50,5,0,1)},
			{XYZW(5,15,0,1),XYZW(10,20,0,1),XYZW(20,25,0,1),XYZW(35,25,0,1),XYZW(45,20,0,1),XYZW(50,15,0,1)},
			{XYZW(5,25,0,1),XYZW(10,30,0,1),XYZW(20,35,0,1),XYZW(35,35,0,1),XYZW(45,30,0,1),XYZW(50,25,0,1)},
			{XYZW(5,35,0,1),XYZW(10,40,0,1),XYZW(20,45,0,1),XYZW(35,45,0,1),XYZW(45,40,0,1),XYZW(50,35,0,1)},
			{XYZW(5,45,0,1),XYZW(10,50,0,1),XYZW(20,55,0,1),XYZW(35,55,0,1),XYZW(45,50,0,1),XYZW(50,45,0,1)},
			{XYZW(5,55,0,1),XYZW(10,60,0,1),XYZW(20,65,0,1),XYZW(35,65,0,1),XYZW(45,60,0,1),XYZW(50,55,0,1)},
		};
		std::vector<double> restkvU;
		std::vector<double> restkvV;
		std::vector<std::vector<XYZW>> newCps;
		NurbsSurface::RemoveKnot(degree, degree, kvU, kvV, cps, 0.66, 1, false, restkvU, restkvV, newCps);
		EXPECT_TRUE(restkvV.size() == kvV.size() - 1);
		EXPECT_TRUE(MathUtils::IsAlmostEqualTo(restkvV[4], 0.33));
		EXPECT_TRUE(MathUtils::IsAlmostEqualTo(restkvV[5], 1.0));
		EXPECT_TRUE(newCps[0].size() == 5);
		EXPECT_TRUE(newCps[0][0].IsAlmostEqualTo(XYZW(5, 5, 0, 1)));
		EXPECT_TRUE(newCps[0][1].IsAlmostEqualTo(XYZW(10, 10, 0, 1)));
		EXPECT_TRUE(newCps[0][2].IsAlmostEqualTo(XYZW(25.15151515, 17.5757575757, 0, 1)));
		EXPECT_TRUE(newCps[0][3].IsAlmostEqualTo(XYZW(45, 10, 0, 1)));
		EXPECT_TRUE(newCps[0][4].IsAlmostEqualTo(XYZW(50, 5, 0, 1)));
	}

	{
		int degree = 3;
		std::vector<double> kv = { 0,0,0,0,1,2,3,4,5,5,5,5 };

		XYZW P0 = XYZW(XYZ(0, 0, 0), 1);
		XYZW P1 = XYZW(XYZ(1, 5, 0), 1);
		XYZW P2 = XYZW(XYZ(2, 10, 0), 1);
		XYZW P3 = XYZW(XYZ(3, 15, 0), 1);
		XYZW P4 = XYZW(XYZ(4, 20, 0), 1);
		XYZW P5 = XYZW(XYZ(5, 15, 0), 1);
		XYZW P6 = XYZW(XYZ(6, 10, 0), 1);
		XYZW P7 = XYZW(XYZ(7, 5, 0), 1);

		std::vector<XYZW> cps = { P0,P1,P2,P3,P4,P5,P6,P7 };
		std::vector<double> updatedKv;
		std::vector<XYZW> updatedCps;
		NurbsCurve::ElevateDegree(degree, kv, cps, 1, updatedKv, updatedCps);
		EXPECT_TRUE(ValidationUtils::IsValidNurbs(degree + 1, updatedKv.size(), updatedCps.size()));
		EXPECT_TRUE(MathUtils::IsAlmostEqualTo(kv[0], updatedKv[0]));
		EXPECT_TRUE(MathUtils::IsAlmostEqualTo(kv[kv.size()-1], updatedKv[updatedKv.size()-1]));
		EXPECT_TRUE(cps[0].IsAlmostEqualTo(updatedCps[0]));
		EXPECT_TRUE(cps[7].IsAlmostEqualTo(updatedCps[updatedCps.size()-1]));
	}

	{
		int degree = 2;
		std::vector<double> kv = { 0,0,0,1,2,3,3,3 };
		std::vector<XYZW> cps = { XYZW(XYZ(0,0,0),1), XYZW(XYZ(1,1,0),4), XYZW(XYZ(3,2,0),1), XYZW(XYZ(4,1,0),1), XYZW(XYZ(5,-1,0),1) };
		std::vector<double> newKv;
		std::vector<XYZ> newCps;

		NurbsCurve::EquallyTessellate(degree, kv, cps, newCps, newKv);
		EXPECT_TRUE(newCps[0].IsAlmostEqualTo(XYZ(0, 0, 0)));
		EXPECT_TRUE(newCps[300].IsAlmostEqualTo(XYZ(5, -1, 0)));
	}

	{
		int degree = 2;
		std::vector<double> kv = { 0,0,0,1,2,3,3,3 };
		std::vector<XYZW> cps = { XYZW(XYZ(0,0,0),1), XYZW(XYZ(1,1,0),4), XYZW(XYZ(3,2,0),1), XYZW(XYZ(4,1,0),1), XYZW(XYZ(5,-1,0),1) };
		XYZ result = NurbsCurve::GetPointOnCurve(degree, kv, 1.0, cps);
		double param = NurbsCurve::GetParamOnCurve(degree, kv, cps, result);
		EXPECT_TRUE(MathUtils::IsAlmostEqualTo(1.0,param));
	}
}