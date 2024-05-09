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

TEST(Test_MathUtils, Binomial)
{
	double bi = MathUtils::Binomial(5, 3);
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(bi, 10.0));
}
