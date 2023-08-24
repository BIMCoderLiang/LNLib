#include "gtest/gtest.h"
#include "BsplineCurve.h"
#include "XYZ.h"
#include "XYZW.h"
#include "MathUtils.h"

using namespace LNLib;

TEST(Test_BsplineCurve, All)
{
	int degree = 2;
	std::vector<double> knotVector = { 0,0,0,1,2,3,4,4,5,5,5 };
	double paramT = 5.0 / 2;
	XYZ P2 = XYZ(5, 6, 7);
	XYZ P3 = XYZ(6, 7, 8);
	XYZ P4 = XYZ(8, 9, 10);
	std::vector<XYZ> controlPoints = { XYZ(1,0,0), XYZ(2,3,4), P2, P3, P4,
										XYZ(9,10,11), XYZ(10,11,12), XYZ(11,12,13)};
	XYZ result = BsplineCurve::GetPointOnCurve(degree, knotVector, paramT, controlPoints);
	EXPECT_TRUE(result.IsAlmostEqualTo(1.0/8 * P2 + 6.0/8 * P3 + 1.0/8 * P4));
	std::vector<XYZ> ders = BsplineCurve::ComputeDerivatives(degree, 1, knotVector, paramT, controlPoints);
	EXPECT_TRUE(ders[1].IsAlmostEqualTo(-0.5 * P2 + 0.5 * P4));
	ders = BsplineCurve::ComputeDerivativesByAllBasisFunctions(degree, 1, knotVector, paramT, controlPoints);
	EXPECT_TRUE(ders[1].IsAlmostEqualTo(-0.5 * P2 + 0.5 * P4));
}