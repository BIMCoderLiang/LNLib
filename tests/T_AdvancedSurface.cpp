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