#include "gtest/gtest.h"
#include "MathUtils.h"
using namespace LNLib;

TEST(Test_MathUtils, Compare)
{
	double a = -10;
	EXPECT_FALSE(MathUtils::IsAlmostEqualTo(a, 0.0));
	EXPECT_TRUE(MathUtils::IsLessThanOrEqual(a, Constants::DoubleEpsilon));
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(a, a + Constants::DoubleEpsilon));
	EXPECT_TRUE(MathUtils::IsGreaterThanOrEqual(a, a + Constants::DoubleEpsilon));
	EXPECT_TRUE(MathUtils::IsLessThanOrEqual(a, a + Constants::DoubleEpsilon));
	EXPECT_TRUE(MathUtils::IsGreaterThan(a, -20));
	EXPECT_TRUE(MathUtils::IsLessThan(a, 0));
}

TEST(Test_MathUtils, Matrix)
{
	std::vector<std::vector<double>> a = { {1,2,4},{3,7,2},{2,3,3} };
	std::vector<std::vector<double>> lower;
	std::vector<std::vector<double>> upper;
	bool result = MathUtils::LUDecomposition(a, lower, upper);
	EXPECT_TRUE(result);
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(lower[0][0], 1) && 
				MathUtils::IsAlmostEqualTo(lower[0][1], 0) &&
				MathUtils::IsAlmostEqualTo(lower[0][2], 0)&&
				MathUtils::IsAlmostEqualTo(lower[1][0], 3)&&
				MathUtils::IsAlmostEqualTo(lower[1][1], 1)&&
				MathUtils::IsAlmostEqualTo(lower[1][2], 0)&&
				MathUtils::IsAlmostEqualTo(lower[2][0], 2)&&
				MathUtils::IsAlmostEqualTo(lower[2][1], -1)&&
				MathUtils::IsAlmostEqualTo(lower[2][2], 1));
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(upper[0][0], 1) &&
				MathUtils::IsAlmostEqualTo(upper[0][1], 2) &&
				MathUtils::IsAlmostEqualTo(upper[0][2], 4) &&
				MathUtils::IsAlmostEqualTo(upper[1][0], 0) &&
				MathUtils::IsAlmostEqualTo(upper[1][1], 1) &&
				MathUtils::IsAlmostEqualTo(upper[1][2], -10) &&
				MathUtils::IsAlmostEqualTo(upper[2][0], 0) &&
				MathUtils::IsAlmostEqualTo(upper[2][1], 0) &&
				MathUtils::IsAlmostEqualTo(upper[2][2], -15));
}

TEST(Test_MathUtils, Binomial)
{
	double bi = MathUtils::Binomial(5, 3);
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(bi, 10.0));
}
