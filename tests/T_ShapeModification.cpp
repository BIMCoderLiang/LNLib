#include "gtest/gtest.h"
#include "NurbsCurve.h"
#include "NurbsSurface.h"
#include "UV.h"
#include "XYZ.h"
#include "XYZW.h"
#include "MathUtils.h"
#include "ValidationUtils.h"
#include "Projection.h"

using namespace LNLib;

TEST(Test_Advanced, All)
{
	{
		int degree = 2;
		std::vector<double> kv = { 0,0,0,1,2,3,3,3 };
		std::vector<XYZW> cps = { XYZW(XYZ(0,0,0),1), XYZW(XYZ(1,1,0),4), XYZW(XYZ(3,2,0),1), XYZW(XYZ(4,1,0),1), XYZW(XYZ(5,-1,0),1) };
		XYZ direction = XYZ(1, 1, 0);
		double distance = 5;
		int moveIndex = 3;
		std::vector<XYZW> updatedCps;
		bool result = NurbsCurve::ControlPointReposition(degree, kv, cps, 2, moveIndex, direction, distance, updatedCps);
		EXPECT_TRUE(result);
	}
	{
		int degree = 2;
		std::vector<double> kv = { 0,0,0,1,2,3,3,3 };
		std::vector<XYZW> cps = { XYZW(XYZ(0,0,0),1), XYZW(XYZ(1,1,0),4), XYZW(XYZ(3,2,0),1), XYZW(XYZ(4,1,0),1), XYZW(XYZ(5,-1,0),1) };
		XYZ direction = XYZ(1, 1, 0);
		double distance = 5;
		std::vector<XYZW> updatedCps;
		int moveIndex = 3;
		NurbsCurve::WeightModification(degree, kv, cps, 2, moveIndex, distance, updatedCps);
		EXPECT_TRUE(updatedCps[moveIndex].IsAlmostEqualTo(cps[moveIndex]));
		EXPECT_FALSE(MathUtils::IsAlmostEqualTo(updatedCps[moveIndex].GetW(), cps[moveIndex].GetW()));
	}
	{
		int degree = 2;
		std::vector<double> kv = { 0,0,0,1,2,3,3,3 };
		std::vector<XYZW> cps = { XYZW(XYZ(0,0,0),1), XYZW(XYZ(1,1,0),4), XYZW(XYZ(3,2,0),1), XYZW(XYZ(4,1,0),1), XYZW(XYZ(5,-1,0),1) };
		XYZ direction = XYZ(1, 1, 0);
		double distance = 5;
		std::vector<XYZW> updatedCps;
		int moveIndex = 3;
		NurbsCurve::NeighborWeightsModification(degree, kv, cps, 2, moveIndex, distance, 2.0, updatedCps);
		EXPECT_FALSE(MathUtils::IsAlmostEqualTo(updatedCps[moveIndex].GetW(), cps[moveIndex].GetW()));
		EXPECT_FALSE(MathUtils::IsAlmostEqualTo(updatedCps[moveIndex + 1].GetW(), cps[moveIndex + 1].GetW()));
	}
	{
		std::vector<double> kv = { 0,0,0,0,1.0/5,2.0/5,3.0/5,4.0/5,1,1,1,1 };
		std::vector<XYZW> cps = { XYZW(0,0,0,1), XYZW(10,0,0,1), XYZW(15,10,0,1), XYZW(12,20,0,1), XYZW(30,25,0,1), XYZW(40,18,0,1), XYZW(36,0,0,1), XYZW(28,0,0,1) };
		std::vector<double> constraintParams = { 0.5 };
		std::vector<XYZ> derivativeConstraint = { XYZ(1,-2,0)};
		std::vector<int> dr = {0};
		std::vector<int> dk = {0};
		std::vector<int> fixed;
		std::vector<XYZW> result = NurbsCurve::ConstraintBasedModification(3, kv, cps, constraintParams, derivativeConstraint, dr, dk, fixed);
		EXPECT_TRUE(ValidationUtils::IsValidNurbs(3, kv.size(), result.size()));
	}
	{
		XYZ origin = XYZ(0, 0, 0);
		XYZ dir = XYZ(1, 0, 0);
		XYZ result = Projection::PointToRay(origin, dir, XYZ(0, 10, 0));
		EXPECT_TRUE(result.IsAlmostEqualTo(origin));
		result = Projection::PointToRay(origin, dir, XYZ(10, 10, 0));
		EXPECT_TRUE(result.IsAlmostEqualTo(XYZ(10,0,0)));
	}

}