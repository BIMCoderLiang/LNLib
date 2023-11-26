#include "gtest/gtest.h"
#include "NurbsCurve.h"
#include "NurbsSurface.h"
#include "UV.h"
#include "XYZ.h"
#include "XYZW.h"
#include "MathUtils.h"
#include "ValidationUtils.h"
#include "Projection.h"

using namespace LNLib;

TEST(Test_Advanced, All)
{
	{
		int degree = 2;
		std::vector<double> kv = { 0,0,0,1,2,3,3,3 };
		std::vector<XYZW> cps = { XYZW(XYZ(0,0,0),1), XYZW(XYZ(1,1,0),4), XYZW(XYZ(3,2,0),1), XYZW(XYZ(4,1,0),1), XYZW(XYZ(5,-1,0),1) };
		XYZ direction = XYZ(1, 1, 0);
		double distance = 5;
		int moveIndex = 3;

		LN_Curve curve;
		curve.Degree = degree;
		curve.KnotVector = kv;
		curve.ControlPoints = cps;

		LN_Curve newCurve;
		bool result = NurbsCurve::ControlPointReposition(curve, 2, moveIndex, direction, distance, newCurve);
		EXPECT_TRUE(result);
	}
	{
		int degree = 2;
		std::vector<double> kv = { 0,0,0,1,2,3,3,3 };
		std::vector<XYZW> cps = { XYZW(XYZ(0,0,0),1), XYZW(XYZ(1,1,0),4), XYZW(XYZ(3,2,0),1), XYZW(XYZ(4,1,0),1), XYZW(XYZ(5,-1,0),1) };
		XYZ direction = XYZ(1, 1, 0);
		double distance = 5;
		int moveIndex = 3;

		LN_Curve curve;
		curve.Degree = degree;
		curve.KnotVector = kv;
		curve.ControlPoints = cps;

		LN_Curve newCurve;
		NurbsCurve::WeightModification(curve, 2, moveIndex, distance, newCurve);

		std::vector<XYZW> updatedCps = newCurve.ControlPoints;
		EXPECT_TRUE(updatedCps[moveIndex].IsAlmostEqualTo(cps[moveIndex]));
		EXPECT_FALSE(MathUtils::IsAlmostEqualTo(updatedCps[moveIndex].GetW(), cps[moveIndex].GetW()));
	}
	{
		int degree = 2;
		std::vector<double> kv = { 0,0,0,1,2,3,3,3 };
		std::vector<XYZW> cps = { XYZW(XYZ(0,0,0),1), XYZW(XYZ(1,1,0),4), XYZW(XYZ(3,2,0),1), XYZW(XYZ(4,1,0),1), XYZW(XYZ(5,-1,0),1) };
		XYZ direction = XYZ(1, 1, 0);
		double distance = 5;
		int moveIndex = 3;

		LN_Curve curve;
		curve.Degree = degree;
		curve.KnotVector = kv;
		curve.ControlPoints = cps;

		LN_Curve newCurve;
		NurbsCurve::NeighborWeightsModification(curve, 2, moveIndex, distance, 2.0, newCurve);

		std::vector<XYZW> updatedCps = newCurve.ControlPoints;
		EXPECT_FALSE(MathUtils::IsAlmostEqualTo(updatedCps[moveIndex].GetW(), cps[moveIndex].GetW()));
		EXPECT_FALSE(MathUtils::IsAlmostEqualTo(updatedCps[moveIndex + 1].GetW(), cps[moveIndex + 1].GetW()));
	}
	{
		std::vector<double> kv = { 0,0,0,0,1.0/5,2.0/5,3.0/5,4.0/5,1,1,1,1 };
		std::vector<XYZW> cps = { XYZW(0,0,0,1), XYZW(10,0,0,1), XYZW(15,10,0,1), XYZW(12,20,0,1), XYZW(30,25,0,1), XYZW(40,18,0,1), XYZW(36,0,0,1), XYZW(28,0,0,1) };
		std::vector<double> constraintParams = { 0.5 };
		std::vector<XYZ> derivativeConstraint = { XYZ(1,-2,0)};
		std::vector<int> dr = {0};
		std::vector<int> dk = {0};
		std::vector<int> fixed;

		LN_Curve curve;
		curve.Degree = 3;
		curve.KnotVector = kv;
		curve.ControlPoints = cps;

		LN_Curve newCurve;
		NurbsCurve::ConstraintBasedModification(curve, constraintParams, derivativeConstraint, dr, dk, fixed, newCurve);
		EXPECT_TRUE(ValidationUtils::IsValidNurbs(newCurve.Degree, newCurve.KnotVector.size(), newCurve.ControlPoints.size()));
	}
	{
		XYZ origin = XYZ(0, 0, 0);
		XYZ dir = XYZ(1, 0, 0);
		XYZ result = Projection::PointToRay(origin, dir, XYZ(0, 10, 0));
		EXPECT_TRUE(result.IsAlmostEqualTo(origin));
		result = Projection::PointToRay(origin, dir, XYZ(10, 10, 0));
		EXPECT_TRUE(result.IsAlmostEqualTo(XYZ(10,0,0)));
	}
	{
		XYZ center = XYZ(0, 0, 0);
		XYZ xAxis = XYZ(1, 0, 0);
		XYZ yAxis = XYZ(0, 1, 0);
		int degree;

		LN_Curve curve;
		NurbsCurve::CreateArc(center, xAxis, yAxis, 0, Constants::Pi, 10, 10, curve);
		std::vector<double> shape = { 0,1,1.5,1,0 };

		LN_Curve newCurve;
		NurbsCurve::Warping(curve, shape, 2, XYZ(0, 0, 1), curve.KnotVector[2], curve.KnotVector[5], newCurve);
		EXPECT_TRUE(curve.ControlPoints.size() == newCurve.ControlPoints.size());
	}
	{
		int degree = 3;
		std::vector<double> kv = { 0,0,0,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1,1,1 };
		std::vector<XYZW> cps = { XYZW(XYZ(-10,2,0),1), XYZW(XYZ(-8,1,0),1), XYZW(XYZ(-6,0,0),1), XYZW(XYZ(-4,-1,0),1),XYZW(XYZ(-2,-2,0),1),XYZW(XYZ(-1,-1,0),1), XYZW(XYZ(0,0,0),1),XYZW(XYZ(1,2,0),1),XYZW(XYZ(2,3,0),1),XYZW(XYZ(4,4,0),1),XYZW(XYZ(6,8,0),1),XYZW(XYZ(7,9,0),1),XYZW(XYZ(8,10,0),1) };
		
		LN_Curve curve;
		curve.Degree = degree;
		curve.KnotVector = kv;
		curve.ControlPoints = cps;

		LN_Curve newCurve;
		auto result = NurbsCurve::Flattening(curve, XYZ(-10, 20, 0), XYZ(10, 20, 0), kv[0], kv[kv.size()-1], newCurve);
		EXPECT_TRUE(result);
	}
	{
		int degree = 3;
		std::vector<double> kv = { 0,0,0,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1,1,1 };
		std::vector<XYZW> cps = { XYZW(XYZ(-10,2,0),1), XYZW(XYZ(-8,1,0),1), XYZW(XYZ(-6,0,0),1), XYZW(XYZ(-4,-1,0),1),XYZW(XYZ(-2,-2,0),1),XYZW(XYZ(-1,-1,0),1), XYZW(XYZ(0,0,0),1),XYZW(XYZ(1,2,0),1),XYZW(XYZ(2,3,0),1),XYZW(XYZ(4,4,0),1),XYZW(XYZ(6,8,0),1),XYZW(XYZ(7,9,0),1),XYZW(XYZ(8,10,0),1) };
		
		LN_Curve curve;
		curve.Degree = degree;
		curve.KnotVector = kv;
		curve.ControlPoints = cps;

		LN_Curve newCurve;
		NurbsCurve::Bending(curve, kv[0], kv[kv.size() - 1], XYZ(0,50,0),10,1.5,newCurve);
		EXPECT_TRUE(curve.ControlPoints.size() == newCurve.ControlPoints.size());
	}
}