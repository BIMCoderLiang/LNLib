#include "gtest/gtest.h"
#include "BsplineSurface.h"
#include "XYZ.h"
#include "XYZW.h"
#include "MathUtils.h"
#include "LNObject.h"

using namespace LNLib;

TEST(Test_BsplineSurface, All)
{
	int degreeU = 2;
	int degreeV = 2;
	std::vector<double> kvU = { 0,0,0,2.0/5,3.0/5,1,1,1 };
	std::vector<double> kvV = { 0,0,0,1.0/5,1.0/2,4.0/5,1,1,1 };

	XYZ P00 = XYZ(0, 0, 0);
	XYZ P01 = XYZ(0, 1, 0);
	XYZ P02 = XYZ(0, 2, 0);
	XYZ P03 = XYZ(0, 3, 0);
	XYZ P04 = XYZ(0, 4, 0);
	XYZ P05 = XYZ(0, 5, 0);
	XYZ P10 = XYZ(1, 0, 0);
	XYZ P11 = XYZ(1, 1, 0);
	XYZ P12 = XYZ(1, 2, 0);
	XYZ P13 = XYZ(1, 3, 0);
	XYZ P14 = XYZ(1, 4, 0);
	XYZ P15 = XYZ(1, 5, 0);
	XYZ P20 = XYZ(2, 0, 0);
	XYZ P21 = XYZ(2, 1, 0);
	XYZ P22 = XYZ(2, 2, 0);
	XYZ P23 = XYZ(2, 3, 0);
	XYZ P24 = XYZ(2, 4, 0);
	XYZ P25 = XYZ(2, 5, 0);
	XYZ P30 = XYZ(3, 0, 0);
	XYZ P31 = XYZ(3, 1, 0);
	XYZ P32 = XYZ(3, 2, 0);
	XYZ P33 = XYZ(3, 3, 0);
	XYZ P34 = XYZ(3, 4, 0);
	XYZ P35 = XYZ(3, 5, 0);
	XYZ P40 = XYZ(4, 0, 0);
	XYZ P41 = XYZ(4, 1, 0);
	XYZ P42 = XYZ(4, 2, 0);
	XYZ P43 = XYZ(4, 3, 0);
	XYZ P44 = XYZ(4, 4, 0);
	XYZ P45 = XYZ(4, 5, 0);

	std::vector<std::vector<XYZ>> controlPoints = {

		{P00, P01, P02, P03, P04, P05},
		{P10, P11, P12, P13, P14, P15},
		{P20, P21, P22, P23, P24, P25},
		{P30, P31, P32, P33, P34, P35},
		{P40, P41, P42, P43, P44, P45},

	};
	
	double u = 1.0 / 5;
	double v = 3.0 / 5;

	LN_BsplineSurface<XYZ> bsplineSurface;
	bsplineSurface.DegreeU = degreeU;
	bsplineSurface.DegreeV = degreeV;
	bsplineSurface.KnotVectorU = kvU;
	bsplineSurface.KnotVectorV = kvV;
	bsplineSurface.ControlPoints = controlPoints;

	XYZ result = BsplineSurface::GetPointOnSurface(bsplineSurface,UV(u,v));

	double N02 = Polynomials::OneBasisFunction(0, degreeU, kvU, u);
	double N12 = Polynomials::OneBasisFunction(1, degreeU, kvU, u);
	double N22 = Polynomials::OneBasisFunction(2, degreeU, kvU, u);

	double N22_ = Polynomials::OneBasisFunction(2, degreeV, kvV, v);
	double N32_ = Polynomials::OneBasisFunction(3, degreeV, kvV, v);
	double N42_ = Polynomials::OneBasisFunction(4, degreeV, kvV, v);

	XYZ a = P02 * N22_ + P03 * N32_ + P04 * N42_;
	XYZ b = P12 * N22_ + P13 * N32_ + P14 * N42_;
	XYZ c = P22 * N22_ + P23 * N32_ + P24 * N42_;

	XYZ checked = N02 * a + N12 * b + N22 * c;
	EXPECT_TRUE(result.IsAlmostEqualTo(checked));

	int derivative = 2;

	std::vector<std::vector<XYZ>> ders = BsplineSurface::ComputeDerivatives(bsplineSurface, derivative, UV(u, v));
	std::vector<std::vector<XYZ>> checkedders = BsplineSurface::ComputeDerivativesByAllBasisFunctions(bsplineSurface, derivative, UV(u, v));
	EXPECT_TRUE(ders[0][0].IsAlmostEqualTo(checkedders[0][0]) &&
				ders[1][0].IsAlmostEqualTo(checkedders[1][0]) &&
				ders[2][0].IsAlmostEqualTo(checkedders[2][0]));
}