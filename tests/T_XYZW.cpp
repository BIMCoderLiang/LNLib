#include "gtest/gtest.h"
#include "XYZW.h"
#include "MathUtils.h"
using namespace LNLib;

TEST(Test_XYZW, Construct)
{
	XYZW xyzw = XYZW();
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(xyzw.GetWX(), 0) &&
				MathUtils::IsAlmostEqualTo(xyzw.GetWY(), 0) &&
				MathUtils::IsAlmostEqualTo(xyzw.GetWZ(), 0) &&
				MathUtils::IsAlmostEqualTo(xyzw.GetW(), 0));
	xyzw = XYZW(XYZ(1, 2, 3), 4);
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(xyzw.GetWX(), 4) &&
				MathUtils::IsAlmostEqualTo(xyzw.GetWY(), 8) &&
				MathUtils::IsAlmostEqualTo(xyzw.GetWZ(), 12) &&
				MathUtils::IsAlmostEqualTo(xyzw.GetW(), 4));
	xyzw = XYZW(4, 8, 12, 4);
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(xyzw.GetWX(), 4) &&
				MathUtils::IsAlmostEqualTo(xyzw.GetWY(), 8) &&
				MathUtils::IsAlmostEqualTo(xyzw.GetWZ(), 12) &&
				MathUtils::IsAlmostEqualTo(xyzw.GetW(), 4));
}

TEST(Test_XYZW, GetSet)
{
	XYZW xyzw = XYZW(4, 8, 12, 4);
	xyzw.SetW(5);
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(xyzw.GetWX(), 5) &&
				MathUtils::IsAlmostEqualTo(xyzw.GetWY(), 10) &&
				MathUtils::IsAlmostEqualTo(xyzw.GetWZ(), 15) &&
				MathUtils::IsAlmostEqualTo(xyzw.GetW(), 5));
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(xyzw.WX(), 5) &&
				MathUtils::IsAlmostEqualTo(xyzw.WY(), 10) &&
				MathUtils::IsAlmostEqualTo(xyzw.WZ(), 15) &&
				MathUtils::IsAlmostEqualTo(xyzw.W(), 5));
}

TEST(Test_XYZW, Common)
{
	XYZW xyzw = XYZW(4, 8, 12, 4);
	XYZ xyz = xyzw.ToXYZ(true);
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(xyz.GetX(), 1) &&
				MathUtils::IsAlmostEqualTo(xyz.GetY(), 2) &&
				MathUtils::IsAlmostEqualTo(xyz.GetZ(), 3));
	XYZ xyz1 = xyzw.ToXYZ(false);
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(xyz1.GetX(), 4) &&
				MathUtils::IsAlmostEqualTo(xyz1.GetY(), 8) &&
				MathUtils::IsAlmostEqualTo(xyz1.GetZ(), 12));
	XYZW xyzw1 = XYZW(5, 10, 15, 5);
	EXPECT_TRUE(xyzw.IsAlmostEqualTo(xyzw1));
	double distance = xyzw.Distance(xyzw1);
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(distance, 3.872983));
}

TEST(Test_XYZW, OperatorOverride)
{
	XYZW xyzw = XYZW(2, 4, 6, 2);
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(xyzw[0], 2) &&
				MathUtils::IsAlmostEqualTo(xyzw[1], 4) &&
				MathUtils::IsAlmostEqualTo(xyzw[2], 6) &&
				MathUtils::IsAlmostEqualTo(xyzw[3], 2));
	XYZW another = XYZW(5, 10, 15, 5);
	XYZW add = xyzw + another;
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(add.GetWX(), 2 + 5) &&
				MathUtils::IsAlmostEqualTo(add.GetWY(), 4 + 10) &&
				MathUtils::IsAlmostEqualTo(add.GetWZ(), 6 + 15) &&
				MathUtils::IsAlmostEqualTo(add.GetW(),  2 + 5));
	XYZW sub = xyzw - another;
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(sub.GetWX(), 2 - 5) &&
				MathUtils::IsAlmostEqualTo(sub.GetWY(), 4 - 10) &&
				MathUtils::IsAlmostEqualTo(sub.GetWZ(), 6 - 15) &&
				MathUtils::IsAlmostEqualTo(sub.GetW(),  2 - 5));
	xyzw += xyzw;
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(xyzw.GetWX(), 4) &&
				MathUtils::IsAlmostEqualTo(xyzw.GetWY(), 8) &&
				MathUtils::IsAlmostEqualTo(xyzw.GetWZ(), 12) &&
				MathUtils::IsAlmostEqualTo(xyzw.GetW(),  4));
	xyzw = XYZW(2, 4, 6, 2);
	XYZW multi = xyzw * 2;
	XYZW multi1 = 2 * xyzw;
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(multi.GetWX(), 4) &&
				MathUtils::IsAlmostEqualTo(multi.GetWY(), 8) &&
				MathUtils::IsAlmostEqualTo(multi.GetWZ(), 12) &&
				MathUtils::IsAlmostEqualTo(multi.GetW(),  4));
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(multi1.GetWX(), 4) &&
				MathUtils::IsAlmostEqualTo(multi1.GetWY(), 8) &&
				MathUtils::IsAlmostEqualTo(multi1.GetWZ(), 12) &&
				MathUtils::IsAlmostEqualTo(multi1.GetW(), 4));
	XYZW divide = xyzw / 2;
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(divide.GetWX(), 1) &&
				MathUtils::IsAlmostEqualTo(divide.GetWY(), 2) &&
				MathUtils::IsAlmostEqualTo(divide.GetWZ(), 3) &&
				MathUtils::IsAlmostEqualTo(divide.GetW(),  1));

}