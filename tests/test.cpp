#include "gtest/gtest.h"
#include "XYZ.h"
#include "XYZW.h"
#include "NurbsCurve.h"
using namespace LNLib;


int main(int argc, char** argv) {
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}

TEST(TestLNLib, GetPointOnCurve)
{
	constexpr auto degree = 3;
	std::vector<XYZW> ctrlPoints
	{
		{ -4, -4, 0, 1},
		{ -2,4,0,1 },
		{ 2,-4,0,1 },
		{ 4,4,0,1 }
	};
	std::vector<double>knots{0, 0, 0, 0, 1, 1, 1, 1};
	XYZ point = NurbsCurve::GetPointOnCurve(degree, knots, ctrlPoints, 0);
	EXPECT_TRUE(point.IsAlmostEqualTo({ -4,-4, 0 }));
}