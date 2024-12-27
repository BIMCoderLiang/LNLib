#include "gtest/gtest.h"
#include "NurbsCurve.h"
#include "NurbsSurface.h"
#include "UV.h"
#include "XYZ.h"
#include "XYZW.h"
#include "MathUtils.h"
#include "ValidationUtils.h"
#include "Interpolation.h"
#include "Intersection.h"
#include "LNObject.h"

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

		LN_NurbsCurve curve;
		NurbsCurve::GlobalInterpolation(degree, Q, curve);
		auto cps = curve.ControlPoints;
		EXPECT_TRUE(cps.size() == Q.size());
		EXPECT_TRUE(cps[0].ToXYZ(true).IsAlmostEqualTo(Q[0]));
		EXPECT_TRUE(cps[4].ToXYZ(true).IsAlmostEqualTo(Q[4]));
	}
	{
		int degree = 3;
		std::vector<XYZ> Q = 
		{
			{282.9136, 219.7436, 0}, {323.8285, 216.8496, 0}, {359.7841, 219.0878, 0},
			{390.926,  225.7011, 0}, {417.4904, 235.9736, 0}, {439.7787, 249.2801, 0},
			{458.0961, 265.1101, 0}, {472.6844, 283.0457, 0}, {483.6853, 302.7034, 0},
			{491.1519, 323.6775, 0}, {495.0864, 345.5107, 0}, {495.4815, 367.6969, 0},
			{492.3456, 389.7038, 0}, {485.7093, 411.0023, 0}, {475.6128, 431.0941, 0},
			{462.0758, 449.5273, 0}, {445.0582, 465.8904, 0}, {424.4307, 479.7753, 0}
		};

		LN_NurbsCurve curve;
		NurbsCurve::GlobalInterpolation(degree, Q, curve);
		NurbsCurve::Check(curve);
	}
	{
		int degree = 2;
		std::vector<XYZ> Q = { XYZ(100,0,0),XYZ(0,100,0),XYZ(-100,0,0),XYZ(0,-100,0)};
		std::vector<XYZ> T = { XYZ(0,1,0),XYZ(-1,0,0),XYZ(0,-1,0),XYZ(1,0,0)};

		LN_NurbsCurve curve;
		NurbsCurve::GlobalInterpolation(degree, Q, T, 1, curve);
		XYZ C0 = NurbsCurve::GetPointOnCurve(curve, 0.0);
		EXPECT_TRUE(C0.IsAlmostEqualTo(XYZ(100, 0, 0)));
		XYZ C1 = NurbsCurve::GetPointOnCurve(curve, 1.0);
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
		LN_NurbsSurface surface;
		NurbsSurface::GlobalInterpolation(Q, degreeU, degreeV, surface);
		XYZ C0 = NurbsSurface::GetPointOnSurface(surface, UV(surface.KnotVectorU[0], surface.KnotVectorV[0]));
		EXPECT_TRUE(C0.IsAlmostEqualTo(P00));
		XYZ C1 = NurbsSurface::GetPointOnSurface(surface, UV(surface.KnotVectorU[0], surface.KnotVectorV[surface.KnotVectorV.size()-1]));
		EXPECT_TRUE(C1.IsAlmostEqualTo(P04));
		XYZ C2 = NurbsSurface::GetPointOnSurface(surface, UV(surface.KnotVectorU[surface.KnotVectorU.size()-1], surface.KnotVectorV[0]));
		EXPECT_TRUE(C2.IsAlmostEqualTo(P70));
		XYZ C3 = NurbsSurface::GetPointOnSurface(surface, UV(surface.KnotVectorU[surface.KnotVectorU.size() - 1], surface.KnotVectorV[surface.KnotVectorV.size() - 1]));
		EXPECT_TRUE(C3.IsAlmostEqualTo(P74));
	}

	{
		std::vector<XYZ> Q = { XYZ(0,0,0),XYZ(3,4,0),XYZ(-1,4,0),XYZ(-4,0,0),XYZ(-4,-3,0) };

		LN_NurbsCurve curve;
		bool result = NurbsCurve::CubicLocalInterpolation(Q, curve);
		EXPECT_TRUE(result);

		auto cps = curve.ControlPoints;
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
		
		LN_NurbsCurve curve;
		bool result = NurbsCurve::LeastSquaresApproximation(3, points, 5, curve);
		EXPECT_TRUE(result);
		auto kv = curve.KnotVector;
		auto cps = curve.ControlPoints;
		EXPECT_TRUE(ValidationUtils::IsValidNurbs(3,kv.size(),cps.size()));
		result = NurbsCurve::LeastSquaresApproximation(3, points, 8, curve);
		EXPECT_TRUE(result);
		kv = curve.KnotVector;
		cps = curve.ControlPoints;
		EXPECT_TRUE(ValidationUtils::IsValidNurbs(3, kv.size(), cps.size()));
		
	}
	{
		XYZ P0 = XYZ(0, 0, 0);
		XYZ P1 = XYZ(10, 10, 0);
		XYZ P2 = XYZ(20, 20, 0);
		XYZ P3 = XYZ(30, 30, 0);
		XYZ P4 = XYZ(40, 40, 0);
		XYZ P5 = XYZ(50, 40, 0);
		XYZ P6 = XYZ(60, 40, 0);
		XYZ P7 = XYZ(70, 40, 0);
		XYZ P8 = XYZ(80, 30, 0);
		XYZ P9 = XYZ(90, 20, 0);
		XYZ P10 = XYZ(100, 10, 0);
		XYZ P11 = XYZ(110, 0, 0);
		XYZ P12 = XYZ(120, -10, 0);
		XYZ P13 = XYZ(130, -20, 0);
		XYZ P14 = XYZ(140, -20, 0);
		XYZ P15 = XYZ(150, -20, 0);
		XYZ P16 = XYZ(160, -20, 0);
		XYZ P17 = XYZ(170, -10, 0);
		XYZ P18 = XYZ(180, 0, 0);
		XYZ P19 = XYZ(190, 10, 0);
		XYZ P20 = XYZ(200, 20, 0);

		std::vector<XYZ> points = { P0,P1,P2,P3,P4,P5,P6,P7,P8,P9,P10,
									P11,P12,P13,P14,P15,P16,P17,P18,P19,P20};
		int m = points.size();

		std::vector<double> pointWeights(m);
		for (int i = 0; i < m; i++)
		{
			pointWeights[i] = 1;
		}
		pointWeights[0] = -1;
		pointWeights[20] = -1;

		std::vector<XYZ> tangents(4);
		tangents[0] = XYZ(0, 1, 0);
		tangents[1] = XYZ(1, -1, 0);
		tangents[2] = XYZ(1, -1, 0);
		tangents[3] = XYZ(0, 1, 0);

		std::vector<int> indices(4);
		indices[0] = 0;
		indices[1] = 7;
		indices[2] = 13;
		indices[3] = m - 1;

		std::vector<double> tangentWeights(4);
		tangentWeights[0] = 1;
		tangentWeights[1] = 1;
		tangentWeights[2] = 1;
		tangentWeights[3] = 1;

		LN_NurbsCurve curve;
		bool result = NurbsCurve::WeightedAndContrainedLeastSquaresApproximation(3, points, pointWeights, tangents, indices, tangentWeights, 4, curve);
		EXPECT_TRUE(result);
		auto kv = curve.KnotVector;
		auto cps = curve.ControlPoints;
		EXPECT_TRUE(ValidationUtils::IsValidNurbs(3, kv.size(), cps.size()));
	}
	{
		int degree = 3;
		std::vector<XYZ> Q = { XYZ(0,0,0),XYZ(3,4,0),XYZ(-1,4,0),XYZ(-4,0,0),XYZ(-4,-3,0) };
		LN_NurbsCurve curve;
		NurbsCurve::GlobalApproximationByErrorBound(degree, Q, 1.5, curve);
		auto kv = curve.KnotVector;
		auto cps = curve.ControlPoints;
		EXPECT_TRUE(ValidationUtils::IsValidNurbs(degree, kv.size(), cps.size()));
	}
	{
		XYZ P00 = XYZ(0,  0, 0);
		XYZ P01 = XYZ(10, 0, 0);
		XYZ P02 = XYZ(20, 0, 0);
		XYZ P03 = XYZ(30, 0, 0);
		XYZ P04 = XYZ(40, 0, 0);
		XYZ P05 = XYZ(50, 0, 0);
		XYZ P06 = XYZ(60, 0, 0);

		XYZ P10 = XYZ(0,  10, 10);
		XYZ P11 = XYZ(10, 10, 20);
		XYZ P12 = XYZ(20, 10, 20);
		XYZ P13 = XYZ(30, 10, 10);
		XYZ P14 = XYZ(40, 10, 30);
		XYZ P15 = XYZ(50, 10, 20);
		XYZ P16 = XYZ(60, 10, 40);

		XYZ P20 = XYZ(0,  20, 20);
		XYZ P21 = XYZ(10, 20, 10);
		XYZ P22 = XYZ(20, 20, 30);
		XYZ P23 = XYZ(30, 20, 20);
		XYZ P24 = XYZ(40, 20, 40);
		XYZ P25 = XYZ(50, 20, 10);
		XYZ P26 = XYZ(60, 20, 20);

		XYZ P30 = XYZ(0,  30, 10);
		XYZ P31 = XYZ(10, 30, 20);
		XYZ P32 = XYZ(20, 30, 40);
		XYZ P33 = XYZ(30, 30, 10);
		XYZ P34 = XYZ(40, 30, 20);
		XYZ P35 = XYZ(50, 30, 30);
		XYZ P36 = XYZ(60, 30, 50);

		XYZ P40 = XYZ(0,  40, 20);
		XYZ P41 = XYZ(10, 40, 30);
		XYZ P42 = XYZ(20, 40, 20);
		XYZ P43 = XYZ(30, 40, 50);
		XYZ P44 = XYZ(40, 40, 20);
		XYZ P45 = XYZ(50, 40, 10);
		XYZ P46 = XYZ(60, 40, 40);

		XYZ P50 = XYZ(0,  50, 40);
		XYZ P51 = XYZ(10, 50, 20);
		XYZ P52 = XYZ(20, 50, 30);
		XYZ P53 = XYZ(30, 50, 10);
		XYZ P54 = XYZ(40, 50, 40);
		XYZ P55 = XYZ(50, 50, 20);
		XYZ P56 = XYZ(60, 50, 30);

		XYZ P60 = XYZ(0,  60, 10);
		XYZ P61 = XYZ(10, 60, 30);
		XYZ P62 = XYZ(20, 60, 20);
		XYZ P63 = XYZ(30, 60, 30);
		XYZ P64 = XYZ(40, 60, 20);
		XYZ P65 = XYZ(50, 60, 10);
		XYZ P66 = XYZ(60, 60, 30);

		int degreeU = 3;
		int degreeV = 3;
		std::vector<std::vector<XYZ>> Q = {
		
			{P00, P01, P02, P03, P04, P05, P06},
			{P10, P11, P12, P13, P14, P15, P16},
			{P20, P21, P22, P23, P24, P25, P26},
			{P30, P31, P32, P33, P34, P35, P36},
			{P40, P41, P42, P43, P44, P45, P46},
			{P50, P51, P52, P53, P54, P55, P56},
		    {P60, P61, P62, P63, P64, P65, P66},
		};
		
		LN_NurbsSurface surface;
		NurbsSurface::GlobalApproximation(Q, degreeU, degreeV, 4, 4, surface);
		EXPECT_TRUE(ValidationUtils::IsValidNurbs(degreeU, surface.KnotVectorU.size(), surface.ControlPoints.size()));
		EXPECT_TRUE(ValidationUtils::IsValidNurbs(degreeV, surface.KnotVectorV.size(), surface.ControlPoints[0].size()));
	}
	{
		XYZ P00 = XYZ(0, 0, 0);
		XYZ P01 = XYZ(10, 0, 0);
		XYZ P02 = XYZ(20, 0, 0);
		XYZ P03 = XYZ(30, 0, 0);
		XYZ P04 = XYZ(40, 0, 0);
		XYZ P05 = XYZ(50, 0, 0);
		XYZ P06 = XYZ(60, 0, 0);
		XYZ P07 = XYZ(70, 0, 0);
		XYZ P08 = XYZ(80, 0, 0);
		XYZ P09 = XYZ(90, 0, 0);

		XYZ P10 = XYZ(0, 10, 10);
		XYZ P11 = XYZ(10, 10, 20);
		XYZ P12 = XYZ(20, 10, 20);
		XYZ P13 = XYZ(30, 10, 10);
		XYZ P14 = XYZ(40, 10, 30);
		XYZ P15 = XYZ(50, 10, 20);
		XYZ P16 = XYZ(60, 10, 40);
		XYZ P17 = XYZ(70, 10, 30);
		XYZ P18 = XYZ(80, 10, 20);
		XYZ P19 = XYZ(90, 10, 10);

		XYZ P20 = XYZ(0, 20, 20);
		XYZ P21 = XYZ(10, 20, 10);
		XYZ P22 = XYZ(20, 20, 30);
		XYZ P23 = XYZ(30, 20, 20);
		XYZ P24 = XYZ(40, 20, 40);
		XYZ P25 = XYZ(50, 20, 10);
		XYZ P26 = XYZ(60, 20, 20);
		XYZ P27 = XYZ(70, 20, 30);
		XYZ P28 = XYZ(80, 20, 30);
		XYZ P29 = XYZ(90, 20, 20);

		XYZ P30 = XYZ(0, 30, 10);
		XYZ P31 = XYZ(10, 30, 20);
		XYZ P32 = XYZ(20, 30, 40);
		XYZ P33 = XYZ(30, 30, 10);
		XYZ P34 = XYZ(40, 30, 20);
		XYZ P35 = XYZ(50, 30, 30);
		XYZ P36 = XYZ(60, 30, 50);
		XYZ P37 = XYZ(70, 30, 40);
		XYZ P38 = XYZ(80, 30, 20);
		XYZ P39 = XYZ(90, 30, 20);

		XYZ P40 = XYZ(0, 40, 20);
		XYZ P41 = XYZ(10, 40, 30);
		XYZ P42 = XYZ(20, 40, 20);
		XYZ P43 = XYZ(30, 40, 50);
		XYZ P44 = XYZ(40, 40, 20);
		XYZ P45 = XYZ(50, 40, 10);
		XYZ P46 = XYZ(60, 40, 40);
		XYZ P47 = XYZ(70, 40, 30);
		XYZ P48 = XYZ(80, 40, 10);
		XYZ P49 = XYZ(90, 40, 20);

		XYZ P50 = XYZ(0, 50, 40);
		XYZ P51 = XYZ(10, 50, 20);
		XYZ P52 = XYZ(20, 50, 30);
		XYZ P53 = XYZ(30, 50, 10);
		XYZ P54 = XYZ(40, 50, 40);
		XYZ P55 = XYZ(50, 50, 20);
		XYZ P56 = XYZ(60, 50, 30);
		XYZ P57 = XYZ(70, 50, 20);
		XYZ P58 = XYZ(80, 50, 10);
		XYZ P59 = XYZ(90, 50, 20);

		XYZ P60 = XYZ(0, 60, 10);
		XYZ P61 = XYZ(10, 60, 30);
		XYZ P62 = XYZ(20, 60, 20);
		XYZ P63 = XYZ(30, 60, 30);
		XYZ P64 = XYZ(40, 60, 20);
		XYZ P65 = XYZ(50, 60, 10);
		XYZ P66 = XYZ(60, 60, 30);
		XYZ P67 = XYZ(70, 60, 20);
		XYZ P68 = XYZ(80, 60, 30);
		XYZ P69 = XYZ(90, 60, 40);

		XYZ P70 = XYZ(0, 70, 10);
		XYZ P71 = XYZ(10,70, 30);
		XYZ P72 = XYZ(20, 70, 20);
		XYZ P73 = XYZ(30, 70, 30);
		XYZ P74 = XYZ(40, 70, 20);
		XYZ P75 = XYZ(50, 70, 10);
		XYZ P76 = XYZ(60, 70, 10);
		XYZ P77 = XYZ(70, 70, 20);
		XYZ P78 = XYZ(80, 70, 20);
		XYZ P79 = XYZ(90, 70, 40);

		XYZ P80 = XYZ(0, 80, 10);
		XYZ P81 = XYZ(10, 80, 30);
		XYZ P82 = XYZ(20, 80, 20);
		XYZ P83 = XYZ(30, 80, 30);
		XYZ P84 = XYZ(40, 80, 20);
		XYZ P85 = XYZ(50, 80, 10);
		XYZ P86 = XYZ(60, 80, 30);
		XYZ P87 = XYZ(70, 80, 10);
		XYZ P88 = XYZ(80, 80, 30);
		XYZ P89 = XYZ(90, 80, 40);

		XYZ P90 = XYZ(0, 90, 10);
		XYZ P91 = XYZ(10, 90, 30);
		XYZ P92 = XYZ(20, 90, 20);
		XYZ P93 = XYZ(30, 90, 30);
		XYZ P94 = XYZ(40, 90, 20);
		XYZ P95 = XYZ(50, 90, 10);
		XYZ P96 = XYZ(60, 90, 30);
		XYZ P97 = XYZ(70, 90, 20);
		XYZ P98 = XYZ(80, 90, 20);
		XYZ P99 = XYZ(90, 90, 40);

		std::vector<std::vector<XYZ>> Q = {

			{P00, P01, P02, P03, P04, P05, P06, P07, P08, P09},
			{P10, P11, P12, P13, P14, P15, P16, P17, P18, P19},
			{P20, P21, P22, P23, P24, P25, P26, P27, P28, P29},
			{P30, P31, P32, P33, P34, P35, P36, P37, P38, P39},
			{P40, P41, P42, P43, P44, P45, P46, P47, P48, P49},
			{P50, P51, P52, P53, P54, P55, P56, P57, P58, P59},
			{P60, P61, P62, P63, P64, P65, P66, P67, P68, P69},
			{P70, P71, P72, P73, P74, P75, P76, P77, P78, P79},
			{P80, P81, P82, P83, P84, P85, P86, P87, P88, P89},
			{P90, P91, P92, P93, P94, P95, P96, P97, P98, P99},
		};
		LN_NurbsSurface surface;
		bool result = NurbsSurface::BicubicLocalInterpolation(Q, surface);
		EXPECT_TRUE(result);
		EXPECT_TRUE(ValidationUtils::IsValidNurbs(3, surface.KnotVectorU.size(), surface.ControlPoints.size()));
		EXPECT_TRUE(ValidationUtils::IsValidNurbs(3, surface.KnotVectorV.size(), surface.ControlPoints[0].size()));
	}
	{
		XYZ P;
		Intersection::ComputeLineAndPlane(XYZ(0, 0, 1), XYZ(0, 0, 0), XYZ(0, 0, 10), XYZ(0, 0, 1), P);
		EXPECT_TRUE(P.IsAlmostEqualTo(XYZ(0, 0, 0)));
	}
	{
		XYZ P0 = XYZ(0, 0, 0);
		XYZ P1 = XYZ(2, 0, 0);
		std::vector<XYZ> Q = { P0, P1};
		std::vector<XYZW> middlePoints;
		bool result = NurbsCurve::FitWithConic(Q, 0, 1, XYZ(1, 1, 0), XYZ(1, -1, 0), 0.1, middlePoints);
		EXPECT_TRUE(result);
	}
	{
		XYZ Px = XYZ(5, -5, 0);
		XYZ P0 = XYZ(10, 0, 0);
		XYZ P3 = XYZ(0, 0, 10);
		std::vector<XYZ> Q = { Px, P0, P3 };
		std::vector<XYZW> middlePoints;
		bool result = NurbsCurve::FitWithCubic(Q, 1, 2, XYZ(0, 1, 0), XYZ(0, -1, 1), 0.1, middlePoints);
		EXPECT_TRUE(result);
	}
}

TEST(Test_Fitting, offset)
{
	// Code copied from issue #25
    std::vector<XYZ> cps_interpolation = {XYZ(XYZ(0, 0, 0)), XYZ(XYZ(1, 1, 0)), XYZ(XYZ(3, 2, 0)), XYZ(XYZ(4, 1, 0)), XYZ(XYZ(5, -1, 0))};
    LN_NurbsCurve curve, new_curve;
    NurbsCurve::GlobalInterpolation(2, cps_interpolation, curve);
	const double offsetDist = 0.1;
    NurbsCurve::Offset(curve, offsetDist, new_curve);

	// Verify a random proportional point.
	double ratio = 0.345;
	double t = curve.KnotVector[0] * (1 - t) + curve.KnotVector[curve.ControlPoints.size()] * t;
	auto point = NurbsCurve::GetPointOnCurve(curve, t);
	t = new_curve.KnotVector[0] * (1 - t) + new_curve.KnotVector[new_curve.ControlPoints.size()] * t;
	auto pointNew = NurbsCurve::GetPointOnCurve(new_curve, t);
	double distance = point.Distance(pointNew);
	EXPECT_NEAR(distance, offsetDist, Constants::DistanceEpsilon);
}
