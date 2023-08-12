#include "gtest/gtest.h"
#include "UV.h"
#include "MathUtils.h"
using namespace LNLib;

TEST(Test_UV, Construct)
{
	UV uv = UV();
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(uv.GetU(), 0.0) &&
				MathUtils::IsAlmostEqualTo(uv.GetV(), 0.0));
	EXPECT_TRUE(uv.IsZero());

	uv = UV(1.5,2.5);
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(uv.GetU(), 1.5) &&
				MathUtils::IsAlmostEqualTo(uv.GetV(), 2.5));
}

TEST(Test_UV, GetSet)
{
	UV uv = UV(1.5, 2.5);
	uv.SetU(2.5);
	uv.SetV(1.5);
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(uv.GetU(), 2.5) &&
				MathUtils::IsAlmostEqualTo(uv.GetV(), 1.5));
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(uv.U(), 2.5) &&
				MathUtils::IsAlmostEqualTo(uv.V(), 1.5));
}

TEST(Test_UV, Common)
{
	UV uv = UV(2, 4);
	UV another = UV(5, 10);
	UV add = uv.Add(another);
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(add.GetU(), 2 + 5) &&
				MathUtils::IsAlmostEqualTo(add.GetV(), 4 + 10));
	UV sub = uv.Substract(another);
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(sub.GetU(), 2 - 5) &&
				MathUtils::IsAlmostEqualTo(sub.GetV(), 4 - 10));
	UV nor0 = UV(1, 0);
	UV nor1 = UV(1, 0);
	EXPECT_TRUE(nor0.IsUnit());
	EXPECT_TRUE(nor1.IsUnit());
	EXPECT_FALSE(uv.IsUnit());
	UV neg = uv.Negative();
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(neg.GetU(), -2) &&
				MathUtils::IsAlmostEqualTo(neg.GetV(), -4));
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(uv.Length(), 4.472135));
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(uv.SqrLength(), 20));
	double distance = uv.Distance(another);
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(distance, 6.708203));
	EXPECT_FALSE(uv.IsAlmostEqualTo(another));
	UV uv1 = UV(2 + Constants::DoubleEpsilon, 4 + Constants::DoubleEpsilon);
	EXPECT_TRUE(uv.IsAlmostEqualTo(uv1));
	double dot = uv.DotProduct(another);
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(dot, 50));
	double cross = uv.CrossProduct(another);
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(cross, 0));
	UV nor2 = uv.Normalize();
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(nor2.Length(), 1));
}

TEST(Test_UV, OperatorOverride)
{
	UV uv = UV(2, 4);
	UV another = UV(5, 10);
	UV add = uv + another;
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(add.GetU(), 2 + 5) &&
				MathUtils::IsAlmostEqualTo(add.GetV(), 4 + 10));
	UV sub = uv - another;
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(sub.GetU(), 2 - 5) &&
				MathUtils::IsAlmostEqualTo(sub.GetV(), 4 - 10));
	uv = another;
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(uv.GetU(), 5) &&
		        MathUtils::IsAlmostEqualTo(uv.GetV(), 10));
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(uv[0], 5) &&
				MathUtils::IsAlmostEqualTo(uv[1], 10));
	uv = -uv;
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(uv.GetU(), -5) &&
				MathUtils::IsAlmostEqualTo(uv.GetV(), -10));
	uv += uv;
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(uv.GetU(), -10) &&
				MathUtils::IsAlmostEqualTo(uv.GetV(), -20));
	uv -= uv;
	EXPECT_TRUE(uv.IsZero());
	uv = UV(2, 4);
	uv *= 2;
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(uv.GetU(), 4) &&
				MathUtils::IsAlmostEqualTo(uv.GetV(), 8));
	uv /= 2;
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(uv.GetU(), 2) &&
				MathUtils::IsAlmostEqualTo(uv.GetV(), 4));
	double dot = uv * another;
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(dot, 50));

	UV multi  = uv * 2;
	UV multi1 = 2 * uv;
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(multi.GetU(), 4) &&
				MathUtils::IsAlmostEqualTo(multi.GetV(), 8));
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(multi1.GetU(), 4) &&
				MathUtils::IsAlmostEqualTo(multi1.GetV(), 8));
	double cross = uv ^ another;
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(cross, 0));
	UV divide = uv / 2;
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(divide.GetU(), 1) &&
				MathUtils::IsAlmostEqualTo(divide.GetV(), 2));
	
}

