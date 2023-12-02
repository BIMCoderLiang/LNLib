#include "gtest/gtest.h"
#include "NurbsCurve.h"
#include "NurbsSurface.h"
#include "LNObject.h"
#include "UV.h"
#include "XYZ.h"
#include "XYZW.h"
#include "MathUtils.h"

using namespace LNLib;

TEST(Test_AdvancedGeometric, All)
{
	{
		int degree = 2;
		std::vector<double> kv = { 0,0,0,1,2,3,3,3 };
		std::vector<XYZW> cps = { XYZW(XYZ(0,0,0),1), XYZW(XYZ(1,1,0),4), XYZW(XYZ(3,2,0),1), XYZW(XYZ(4,1,0),1), XYZW(XYZ(5,-1,0),1) };

		LN_Curve curve;
		curve.Degree = degree;
		curve.KnotVector = kv;
		curve.ControlPoints = cps;

		XYZ result = NurbsCurve::GetPointOnCurve(curve, 1.0);
		double param = NurbsCurve::GetParamOnCurve(curve, result);
		EXPECT_TRUE(MathUtils::IsAlmostEqualTo(1.0,param));
	}

	{
		int degreeU = 2;
		int degreeV = 2;
		std::vector<double> kvU = { 0,0,0,1,2,3,4,4,5,5,5 };
		std::vector<double> kvV = { 0,0,0,1,2,3,3,3 };
		UV uv = UV(5.0 / 2, 1);

		XYZW P20 = XYZW(-1, 2, 4, 1);
		XYZW P21 = XYZW(0, 2, 4, 1);
		XYZW P22 = XYZW(0, 6, 4, 2);
		XYZW P23 = XYZW(0, 2, 0, 1);
		XYZW P24 = XYZW(1, 2, 0, 1);

		XYZW P10 = 0.9 * P20;
		XYZW P11 = 0.9 * P21;
		XYZW P12 = 0.9 * P22;
		XYZW P13 = 0.9 * P23;
		XYZW P14 = 0.9 * P24;

		XYZW P00 = 0.9 * P10;
		XYZW P01 = 0.9 * P11;
		XYZW P02 = 0.9 * P12;
		XYZW P03 = 0.9 * P13;
		XYZW P04 = 0.9 * P14;

		XYZW P30 = XYZW(3, 6, 8, 2);
		XYZW P31 = XYZW(4, 6, 8, 2);
		XYZW P32 = XYZW(12, 24, 12, 6);
		XYZW P33 = XYZW(4, 6, 0, 2);
		XYZW P34 = XYZW(5, 6, 0, 2);

		XYZW P40 = XYZW(3, 2, 4, 1);
		XYZW P41 = XYZW(4, 2, 4, 1);
		XYZW P42 = XYZW(8, 6, 4, 2);
		XYZW P43 = XYZW(4, 2, 0, 1);
		XYZW P44 = XYZW(5, 2, 0, 1);

		XYZW P50 = 1.5 * P40;
		XYZW P51 = 1.5 * P41;
		XYZW P52 = 1.5 * P42;
		XYZW P53 = 1.5 * P43;
		XYZW P54 = 1.5 * P44;

		XYZW P60 = 1.5 * P50;
		XYZW P61 = 1.5 * P51;
		XYZW P62 = 1.5 * P52;
		XYZW P63 = 1.5 * P53;
		XYZW P64 = 1.5 * P54;

		XYZW P70 = 1.5 * P60;
		XYZW P71 = 1.5 * P61;
		XYZW P72 = 1.5 * P62;
		XYZW P73 = 1.5 * P63;
		XYZW P74 = 1.5 * P64;

		std::vector<std::vector<XYZW>> cps = {

			{P00, P01, P02, P03, P04},
			{P10, P11, P12, P13, P14},
			{P20, P21, P22, P23, P24},
			{P30, P31, P32, P33, P34},
			{P40, P41, P42, P43, P44},
			{P50, P51, P52, P53, P54},
			{P60, P61, P62, P63, P64},
			{P70, P71, P72, P73, P74},

		};
		LN_Surface surface;
		surface.DegreeU = degreeU;
		surface.DegreeV = degreeV;
		surface.KnotVectorU = kvU;
		surface.KnotVectorV = kvV;
		surface.ControlPoints = cps;
		XYZ result = NurbsSurface::GetPointOnSurface(surface, uv);
		UV param = NurbsSurface::GetParamOnSurface(surface, result);
		EXPECT_TRUE(param.IsAlmostEqualTo(uv));
	}

	{
		int degree = 3;
		std::vector<double> kv = { 0,0,0,0,3.0/10,7.0/10,1,1,1,1 };
		XYZW P0 = XYZW(XYZ(0, 0, 0), 1);
		XYZW P1 = XYZW(XYZ(10, 10, 0), 1);
		XYZW P2 = XYZW(XYZ(20, 20, 0), 1);
		XYZW P3 = XYZW(XYZ(30, 30, 0), 1);
		XYZW P4 = XYZW(XYZ(40, 20, 0), 1);
		XYZW P5 = XYZW(XYZ(50, 10, 0), 1);
		std::vector<XYZW> cps = {P0,P1,P2,P3,P4,P5};

		double alpha = 2;
		double beta = 1;
		double gamma = 3;
		double delta = 2;

		LN_Curve curve;
		curve.Degree = degree;
		curve.KnotVector = kv;
		curve.ControlPoints = cps;

		LN_Curve newtc;
		NurbsCurve::Reparametrize(curve, alpha, beta, gamma, delta, newtc);

		std::vector<double> newKv = newtc.KnotVector;
		std::vector<XYZW> newCps = newtc.ControlPoints;

		EXPECT_TRUE(MathUtils::IsAlmostEqualTo(newKv[0], 0.5));
		EXPECT_TRUE(MathUtils::IsAlmostEqualTo(newKv[1], 0.5));
		EXPECT_TRUE(MathUtils::IsAlmostEqualTo(newKv[2], 0.5));
		EXPECT_TRUE(MathUtils::IsAlmostEqualTo(newKv[3], 0.5));
		EXPECT_TRUE(MathUtils::IsAlmostEqualTo(newKv[4], 0.551724));
		EXPECT_TRUE(MathUtils::IsAlmostEqualTo(newKv[5], 0.585365));
		EXPECT_TRUE(MathUtils::IsAlmostEqualTo(newKv[6], 0.6));
		EXPECT_TRUE(MathUtils::IsAlmostEqualTo(newKv[7], 0.6));
		EXPECT_TRUE(MathUtils::IsAlmostEqualTo(newKv[8], 0.6));
		EXPECT_TRUE(MathUtils::IsAlmostEqualTo(newKv[9], 0.6));

		EXPECT_TRUE(MathUtils::IsAlmostEqualTo(newCps[0].GetW(), 0.125));
		EXPECT_TRUE(MathUtils::IsAlmostEqualTo(newCps[1].GetW(), 0.0862));
		EXPECT_TRUE(MathUtils::IsAlmostEqualTo(newCps[2].GetW(), 0.04205));
		EXPECT_TRUE(MathUtils::IsAlmostEqualTo(newCps[3].GetW(), 0.01682));
		EXPECT_TRUE(MathUtils::IsAlmostEqualTo(newCps[4].GetW(), 0.00975));
		EXPECT_TRUE(MathUtils::IsAlmostEqualTo(newCps[5].GetW(), 0.008));
	}
}