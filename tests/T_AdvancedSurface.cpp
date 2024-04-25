#include "gtest/gtest.h"
#include "NurbsCurve.h"
#include "NurbsSurface.h"
#include "KnotVectorUtils.h"
#include "UV.h"
#include "XYZ.h"
#include "XYZW.h"
#include "MathUtils.h"

using namespace LNLib;

TEST(Test_AdvancedSurface, KnotVector)
{
	std::vector<double> kv0 = { 0,0,0,1.0 / 4,1.0 / 2,3.0 / 4,1,1,1 };
	std::vector<double> kv1 = { 0,0,0,0,3.0 / 10,7.0 / 10,1,1,1,1 };

	std::vector<std::vector<double>> kvs;
	kvs.emplace_back(kv0);
	kvs.emplace_back(kv1);
	auto result = KnotVectorUtils::GetInsertedKnotElements(kvs);
	EXPECT_TRUE(result[0].size() == 4);
	EXPECT_TRUE(result[1].size() == 3);
}

TEST(Test_AdvancedSurface, Bilinear)
{
	XYZ p00 = XYZ(0, 0, 0);
	XYZ p01 = XYZ(2, 0.5, 0);
	XYZ p10 = XYZ(0.2, 1.0, 0);
	XYZ p11 = XYZ(1.8, 1.5, 0);

	LN_NurbsSurface surface;
	NurbsSurface::CreateBilinearSurface(p00, p01, p10, p11, surface);
	XYZ point = NurbsSurface::GetPointOnSurface(surface, UV(0.5, 0.5));
	EXPECT_TRUE(point.IsAlmostEqualTo(XYZ(1.0, 0.75, 0.0)));
}