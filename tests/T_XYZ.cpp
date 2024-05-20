#include "gtest/gtest.h"
#include "XYZ.h"
#include "MathUtils.h"
using namespace LNLib;

TEST(Test_XYZ, Construct)
{
	XYZ xyz = XYZ();
	EXPECT_TRUE(xyz.IsZero());
	xyz = XYZ(1, 2, 3);
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(xyz.GetX(), 1) &&
				MathUtils::IsAlmostEqualTo(xyz.GetY(), 2) &&
				MathUtils::IsAlmostEqualTo(xyz.GetZ(), 3));
}

TEST(Test_XYZ, GetSet)
{
	XYZ xyz = XYZ(1.5, 2.5, 3.5);
	xyz.SetX(2.5);
	xyz.SetY(1.5);
	xyz.SetZ(4.5);
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(xyz.GetX(), 2.5) &&
				MathUtils::IsAlmostEqualTo(xyz.GetY(), 1.5) &&
				MathUtils::IsAlmostEqualTo(xyz.GetZ(), 4.5));
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(xyz.X(), 2.5) &&
				MathUtils::IsAlmostEqualTo(xyz.Y(), 1.5) &&
				MathUtils::IsAlmostEqualTo(xyz.Z(), 4.5));
}

TEST(Test_XYZ, Common)
{
	XYZ xyz = XYZ(1, 0, 0);
	EXPECT_TRUE(xyz.IsUnit());
	xyz = XYZ(0, 1, 0);
	EXPECT_TRUE(xyz.IsUnit());
	xyz = XYZ(0, 0, 1);
	EXPECT_TRUE(xyz.IsUnit());
	xyz = XYZ(2, 4, 6);
	XYZ another = XYZ(2 + Constants::DoubleEpsilon, 4 + Constants::DoubleEpsilon, 6 + Constants::DoubleEpsilon);
	EXPECT_TRUE(xyz.IsAlmostEqualTo(another));
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(xyz.Length(), 7.4833147));
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(xyz.SqrLength(), 56));
	xyz = XYZ(1, 0, 0);
	another = XYZ(0, 1, 0);
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(xyz.AngleTo(another), Constants::Pi / 2));
	xyz = XYZ(2, 4, 6);
	XYZ nor = xyz.Normalize();
	EXPECT_TRUE(nor.IsUnit());
	XYZ add = xyz.Add(another);
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(add.GetX(), 2) &&
				MathUtils::IsAlmostEqualTo(add.GetY(), 4 + 1) &&
				MathUtils::IsAlmostEqualTo(add.GetZ(), 6));
	XYZ sub = xyz.Substract(another);
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(sub.GetX(), 2) &&
				MathUtils::IsAlmostEqualTo(sub.GetY(), 4 - 1) &&
				MathUtils::IsAlmostEqualTo(sub.GetZ(), 6));
	XYZ neg = xyz.Negative();
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(neg.GetX(), -2) &&
				MathUtils::IsAlmostEqualTo(neg.GetY(), -4) &&
				MathUtils::IsAlmostEqualTo(neg.GetZ(), -6));
	double dot = xyz.DotProduct(another);
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(dot, 4));
	XYZ cross = xyz.CrossProduct(another);
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(cross.GetX(), -6) &&
				MathUtils::IsAlmostEqualTo(cross.GetY(), 0) &&
				MathUtils::IsAlmostEqualTo(cross.GetZ(), 2));
	double distance = xyz.Distance(another);
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(distance, 7));
}

TEST(Test_XYZ, OperatorOverride)
{
	XYZ xyz = XYZ(2, 4, 6);
	XYZ another = XYZ(5, 10, 15);
	XYZ add = xyz + another;
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(add.GetX(), 2 + 5) &&
				MathUtils::IsAlmostEqualTo(add.GetY(), 4 + 10) &&
				MathUtils::IsAlmostEqualTo(add.GetZ(), 6 + 15));
	XYZ sub = xyz - another;
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(sub.GetX(), 2 - 5) &&
				MathUtils::IsAlmostEqualTo(sub.GetY(), 4 - 10) &&
				MathUtils::IsAlmostEqualTo(sub.GetZ(), 6 - 15));
	xyz = another;
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(xyz.GetX(), 5)  &&
				MathUtils::IsAlmostEqualTo(xyz.GetY(), 10) &&
				MathUtils::IsAlmostEqualTo(xyz.GetZ(), 15));
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(xyz[0], 5)  &&
				MathUtils::IsAlmostEqualTo(xyz[1], 10) &&
				MathUtils::IsAlmostEqualTo(xyz[2], 15));
	xyz = -xyz;
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(xyz.GetX(), -5)  &&
				MathUtils::IsAlmostEqualTo(xyz.GetY(), -10) &&
				MathUtils::IsAlmostEqualTo(xyz.GetZ(), -15));
	xyz += xyz;
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(xyz.GetX(), -10) &&
				MathUtils::IsAlmostEqualTo(xyz.GetY(), -20) &&
				MathUtils::IsAlmostEqualTo(xyz.GetZ(), -30));
	xyz -= xyz;
	EXPECT_TRUE(xyz.IsZero());
	xyz = XYZ(2, 4, 6);
	xyz *= 2;
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(xyz.GetX(), 4) &&
				MathUtils::IsAlmostEqualTo(xyz.GetY(), 8) &&
				MathUtils::IsAlmostEqualTo(xyz.GetZ(), 12));
	xyz /= 2;
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(xyz.GetX(), 2) &&
				MathUtils::IsAlmostEqualTo(xyz.GetY(), 4) &&
				MathUtils::IsAlmostEqualTo(xyz.GetZ(), 6));
	double dot = xyz * another;
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(dot, 140));

	XYZ multi = xyz * 2;
	XYZ multi1 = 2 * xyz;
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(multi.GetX(), 4) &&
				MathUtils::IsAlmostEqualTo(multi.GetY(), 8) &&
				MathUtils::IsAlmostEqualTo(multi.GetZ(), 12));
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(multi1.GetX(), 4) &&
				MathUtils::IsAlmostEqualTo(multi1.GetY(), 8) &&
				MathUtils::IsAlmostEqualTo(multi1.GetZ(), 12));
	XYZ cross = xyz ^ another;
	EXPECT_TRUE(cross.IsZero());
	XYZ divide = xyz / 2;
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(divide.GetX(), 1) &&
				MathUtils::IsAlmostEqualTo(divide.GetY(), 2) &&
				MathUtils::IsAlmostEqualTo(divide.GetZ(), 3));
}

TEST(Test_XYZ, Orthogonal)
{
	XYZ current = XYZ(3, 1, 2);
	for (int i = 0; i < 100; i++)
	{
		XYZ orthogonal = XYZ::CreateRandomOrthogonal(current);
		double result = orthogonal.DotProduct(current);
		EXPECT_TRUE(MathUtils::IsAlmostEqualTo(result, 0.0));
	}
}