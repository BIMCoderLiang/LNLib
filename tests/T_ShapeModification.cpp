#include "gtest/gtest.h"
#include "NurbsCurve.h"
#include "NurbsSurface.h"
#include "UV.h"
#include "XYZ.h"
#include "XYZW.h"
#include "MathUtils.h"
#include "ValidationUtils.h"

using namespace LNLib;

TEST(Test_Advanced, All)
{
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
}