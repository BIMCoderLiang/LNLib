﻿#include "gtest/gtest.h"
#include "NurbsCurve.h"
#include "XYZ.h"
#include "XYZW.h"
#include "MathUtils.h"
#include "Intersection.h"

using namespace LNLib;

TEST(Test_Circles, All)
{
	{
		XYZ sp1 = XYZ(5.0, 181.34, 0);
		XYZ ep1 = XYZ(13.659999999999997, 176.34, 0);

		XYZ sp2 = XYZ(19.999779996773235, 189.9998729810778, 0);
		XYZ ep2 = XYZ(19.999652977851035, 180.00009298430456, 0);

		double param1, param2 = 0.0;
		XYZ interP;
		CurveCurveIntersectionType type = Intersection::ComputeRays(sp1, ep1 - sp1, sp2, ep2 - sp2, param1, param2, interP);
		EXPECT_TRUE(type != CurveCurveIntersectionType::Skew);
	}

	{
		XYZ center = XYZ(0, 0, 0);
		XYZ xAxis = XYZ(1, 0, 0);
		XYZ yAxis = XYZ(0, 1, 0);
		int degree;
		std::vector<double> kv;
		std::vector<XYZW> cps;
		bool result = NurbsCurve::CreateArc(center, xAxis, yAxis, 0, 2 * Constants::Pi, 10, 10, degree, kv, cps);
		EXPECT_TRUE(result);
		XYZ C0 = NurbsCurve::GetPointOnCurve(degree, kv, 0, cps);
		EXPECT_TRUE(C0.IsAlmostEqualTo(XYZ(10, 0, 0)));
		XYZ C1 = NurbsCurve::GetPointOnCurve(degree, kv, 0.25, cps);
		EXPECT_TRUE(C1.IsAlmostEqualTo(XYZ(0, 10, 0)));
		XYZ C2 = NurbsCurve::GetPointOnCurve(degree, kv, 0.5, cps);
		EXPECT_TRUE(C2.IsAlmostEqualTo(XYZ(-10, 0, 0)));
		XYZ C3 = NurbsCurve::GetPointOnCurve(degree, kv, 0.75, cps);
		EXPECT_TRUE(C3.IsAlmostEqualTo(XYZ(0, -10, 0)));
	}

	{
		XYZ start = XYZ(0, 10, 0);
		XYZ startTangent = XYZ(1, 0, 0);
		XYZ end = XYZ(0, -10, 0);
		XYZ endTangent = XYZ(-1, 0, 0);
		XYZ pointOnConic = XYZ(10, 0, 0);

		int degree;
		std::vector<double> kv;
		std::vector<XYZW> cps;
		bool result = NurbsCurve::CreateOpenConic(start, startTangent, end, endTangent, pointOnConic, degree, kv, cps);
		EXPECT_TRUE(result);
		XYZ C0 = NurbsCurve::GetPointOnCurve(degree, kv, 0, cps);
		EXPECT_TRUE(C0.IsAlmostEqualTo(XYZ(0, 10, 0)));
		XYZ C1 = NurbsCurve::GetPointOnCurve(degree, kv, 0.5, cps);
		EXPECT_TRUE(C1.IsAlmostEqualTo(XYZ(10, 0, 0)));
		XYZ C2 = NurbsCurve::GetPointOnCurve(degree, kv, 1, cps);
		EXPECT_TRUE(C2.IsAlmostEqualTo(XYZ(0, -10, 0)));
	}
}