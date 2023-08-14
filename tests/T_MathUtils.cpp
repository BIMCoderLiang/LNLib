#include "gtest/gtest.h"
#include "MathUtils.h"
using namespace LNLib;

TEST(Test_MathUtils, Matrix)
{
	std::vector<std::vector<double>> a = { {1,1,0},{2,1,3},{3,1,1} };
	std::vector<std::vector<double>> lower;
	std::vector<std::vector<double>> upper;
	bool result = MathUtils::LUDecomposition(a, lower, upper);
	EXPECT_TRUE(result);
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(lower[0][0], 1) && 
				MathUtils::IsAlmostEqualTo(lower[0][1], 0) &&
				MathUtils::IsAlmostEqualTo(lower[0][2], 0)&&
				MathUtils::IsAlmostEqualTo(lower[1][0], 2)&&
				MathUtils::IsAlmostEqualTo(lower[1][1], -1)&&
				MathUtils::IsAlmostEqualTo(lower[1][2], 0)&&
				MathUtils::IsAlmostEqualTo(lower[2][0], 3)&&
				MathUtils::IsAlmostEqualTo(lower[2][1], -2)&&
				MathUtils::IsAlmostEqualTo(lower[2][2], -5));
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(upper[0][0], 1) &&
				MathUtils::IsAlmostEqualTo(upper[0][1], 1) &&
				MathUtils::IsAlmostEqualTo(upper[0][2], 0) &&
				MathUtils::IsAlmostEqualTo(upper[1][0], 0) &&
				MathUtils::IsAlmostEqualTo(upper[1][1], 1) &&
				MathUtils::IsAlmostEqualTo(upper[1][2], -3) &&
				MathUtils::IsAlmostEqualTo(upper[2][0], 0) &&
				MathUtils::IsAlmostEqualTo(upper[2][1], 0) &&
				MathUtils::IsAlmostEqualTo(upper[2][2], 1));
}