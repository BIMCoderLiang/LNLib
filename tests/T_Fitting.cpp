#include "gtest/gtest.h"
#include "NurbsCurve.h"
#include "NurbsSurface.h"
#include "UV.h"
#include "XYZ.h"
#include "XYZW.h"
#include "MathUtils.h"
#include "ValidationUtils.h"
#include "Interpolation.h"

using namespace LNLib;

TEST(Test_Fitting, Interpolation)
{
	{
		int degree = 3;
		std::vector<XYZ> Q = { XYZ(0,0,0),XYZ(3,4,0),XYZ(-1,4,0),XYZ(-4,0,0),XYZ(-4,-3,0) };
		auto params = Interpolation::GetChordParameterization(Q);
		auto result = Interpolation::AverageKnotVector(degree, params);
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

	{
		int degreeU = 2;
		int degreeV = 2;

		XYZ P20 = XYZ(-1, 2, 4);
		XYZ P21 = XYZ(0, 2, 4);
		XYZ P22 = XYZ(0, 6, 4);
		XYZ P23 = XYZ(0, 2, 0);
		XYZ P24 = XYZ(1, 2, 0);

		XYZ P10 = 0.9 * P20;
		XYZ P11 = 0.9 * P21;
		XYZ P12 = 0.9 * P22;
		XYZ P13 = 0.9 * P23;
		XYZ P14 = 0.9 * P24;

		XYZ P00 = 0.9 * P10;
		XYZ P01 = 0.9 * P11;
		XYZ P02 = 0.9 * P12;
		XYZ P03 = 0.9 * P13;
		XYZ P04 = 0.9 * P14;

		XYZ P30 = XYZ(3, 6, 8);
		XYZ P31 = XYZ(4, 6, 8);
		XYZ P32 = XYZ(12, 24, 12);
		XYZ P33 = XYZ(4, 6, 0);
		XYZ P34 = XYZ(5, 6, 0);

		XYZ P40 = XYZ(3, 2, 4);
		XYZ P41 = XYZ(4, 2, 4);
		XYZ P42 = XYZ(8, 6, 4);
		XYZ P43 = XYZ(4, 2, 0);
		XYZ P44 = XYZ(5, 2, 0);

		XYZ P50 = 1.5 * P40;
		XYZ P51 = 1.5 * P41;
		XYZ P52 = 1.5 * P42;
		XYZ P53 = 1.5 * P43;
		XYZ P54 = 1.5 * P44;

		XYZ P60 = 1.5 * P50;
		XYZ P61 = 1.5 * P51;
		XYZ P62 = 1.5 * P52;
		XYZ P63 = 1.5 * P53;
		XYZ P64 = 1.5 * P54;

		XYZ P70 = 1.5 * P60;
		XYZ P71 = 1.5 * P61;
		XYZ P72 = 1.5 * P62;
		XYZ P73 = 1.5 * P63;
		XYZ P74 = 1.5 * P64;

		std::vector<std::vector<XYZ>> Q = {

			{P00, P01, P02, P03, P04},
			{P10, P11, P12, P13, P14},
			{P20, P21, P22, P23, P24},
			{P30, P31, P32, P33, P34},
			{P40, P41, P42, P43, P44},
			{P50, P51, P52, P53, P54},
			{P60, P61, P62, P63, P64},
			{P70, P71, P72, P73, P74},

		};
		std::vector<double> kvU, kvV;
		std::vector<std::vector<XYZW>> cps;
		NurbsSurface::GlobalInterpolation(Q, degreeU, degreeV, kvU, kvV, cps);
		XYZ C0 = NurbsSurface::GetPointOnSurface(degreeU, degreeV, kvU, kvV, UV(kvU[0], kvV[0]), cps);
		EXPECT_TRUE(C0.IsAlmostEqualTo(P00));
		XYZ C1 = NurbsSurface::GetPointOnSurface(degreeU, degreeV, kvU, kvV, UV(kvU[0], kvV[kvV.size()-1]), cps);
		EXPECT_TRUE(C1.IsAlmostEqualTo(P04));
		XYZ C2 = NurbsSurface::GetPointOnSurface(degreeU, degreeV, kvU, kvV, UV(kvU[kvU.size()-1], kvV[0]), cps);
		EXPECT_TRUE(C2.IsAlmostEqualTo(P70));
		XYZ C3 = NurbsSurface::GetPointOnSurface(degreeU, degreeV, kvU, kvV, UV(kvU[kvU.size() - 1], kvV[kvV.size() - 1]), cps);
		EXPECT_TRUE(C3.IsAlmostEqualTo(P74));
	}

	{
		std::vector<XYZ> Q = { XYZ(0,0,0),XYZ(3,4,0),XYZ(-1,4,0),XYZ(-4,0,0),XYZ(-4,-3,0) };
		std::vector<double> kv;
		std::vector<XYZW> cps;
		bool result = NurbsCurve::CubicLocalInterpolation(Q, kv, cps);
		EXPECT_TRUE(result);
		EXPECT_TRUE(cps[0].ToXYZ(true).IsAlmostEqualTo(Q[0]));
		EXPECT_TRUE(cps[cps.size()-1].ToXYZ(true).IsAlmostEqualTo(Q[4]));
	}
}

TEST(Test_Fitting, Approximation)
{
	{
		XYZ P0 = XYZ(20, 20, 0);
		XYZ P1 = XYZ(20, 80, 0);
		XYZ P2 = XYZ(20, 120, 0);
		XYZ P3 = XYZ(20, 160, 0);
		XYZ P4 = XYZ(80, 200, 0);
		XYZ P5 = XYZ(120, 200, 0);
		XYZ P6 = XYZ(160, 160, 0);
		XYZ P7 = XYZ(160, 120, 0);
		XYZ P8 = XYZ(120, 80, 0);
		XYZ P9 = XYZ(80, 80, 0);

		std::vector<XYZ> points = { P0,P1,P2,P3,P4,P5,P6,P7,P8,P9 };
		std::vector<double> kv;
		std::vector<XYZW> cps;
		bool result = NurbsCurve::LeastSquaresApproximation(3, points, 5, kv, cps);
		EXPECT_TRUE(result);
		EXPECT_TRUE(ValidationUtils::IsValidNurbs(3,kv.size(),cps.size()));
		kv.clear();
		cps.clear();
		result = NurbsCurve::LeastSquaresApproximation(3, points, 8, kv, cps);
		EXPECT_TRUE(result);
		EXPECT_TRUE(ValidationUtils::IsValidNurbs(3, kv.size(), cps.size()));
		
	}
	{
		XYZ P0 = XYZ(20, 20, 0);
		XYZ P1 = XYZ(20, 80, 0);
		XYZ P2 = XYZ(20, 120, 0);
		XYZ P3 = XYZ(20, 160, 0);
		XYZ P4 = XYZ(80, 200, 0);
		XYZ P5 = XYZ(120, 200, 0);
		XYZ P6 = XYZ(160, 160, 0);
		XYZ P7 = XYZ(160, 120, 0);
		XYZ P8 = XYZ(120, 80, 0);
		XYZ P9 = XYZ(80, 80, 0);
		XYZ P10 = XYZ(60, 60, 0);

		std::vector<XYZ> points = { P0,P1,P2,P3,P4,P5,P6,P7,P8,P9,P10 };
		int m = points.size();
		int mm = m - 1;
		std::vector<int> Indices(2);
		std::vector<XYZ> D(2);
		std::vector<double> wd(2), wp(m);
		wp[0] = wp[mm] = -1.0;
		for (int i = 1; i < mm; i++)
		{
			wp[i] = 1.0;
		}
			
		wd[0] = wd[1] = -1.0;
		Indices[0] = 0;
		Indices[1] = mm;
		D[0] = XYZ(1,0,0);
		D[1] = XYZ(0,1,0);

		double length = 0;
		for (int i = 1; i < m; i++)
		{
			length += points[i - 1].Distance(points[i]);
		}

		D[0] = D[0].Normalize();
		D[1] = D[1].Normalize();

		std::vector<double> kv;
		std::vector<XYZW> cps;
		bool result = NurbsCurve::WeightedAndContrainedLeastSquaresApproximation(3, points, wp, D, Indices, wd, 5, kv, cps);
		EXPECT_TRUE(result);
		EXPECT_TRUE(ValidationUtils::IsValidNurbs(3, kv.size(), cps.size()));
	}
	{
		int degree = 3;
		std::vector<XYZ> Q = { XYZ(0,0,0),XYZ(3,4,0),XYZ(-1,4,0),XYZ(-4,0,0),XYZ(-4,-3,0) };
		std::vector<double> kv;
		std::vector<XYZW> cps;
		NurbsCurve::GlobalCurveApproximationByErrorBound(degree, Q, 1.5, kv, cps);
		EXPECT_TRUE(ValidationUtils::IsValidNurbs(degree, kv.size(), cps.size()));
	}
}