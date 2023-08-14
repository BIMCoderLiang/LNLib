#include "gtest/gtest.h"
#include "Matrix4d.h"
#include "XYZ.h"
#include "MathUtils.h"
using namespace LNLib;

TEST(Test_Matrix4d, Construct)
{
	Matrix4d m = Matrix4d();
	EXPECT_TRUE(m.IsIdentity());
	XYZ x = XYZ(1, 0, 0);
	XYZ y = XYZ(0, 1, 0);
	XYZ z = XYZ(0, 0, 1);
	XYZ w = XYZ(0, 0, 0);
	Matrix4d m1 = Matrix4d(x,y,z,w);
	EXPECT_TRUE(m1.IsIdentity());
	Matrix4d m2 = Matrix4d(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1);
	m = m2;
	bool result = true;
	for (int i = 0; i <= 3; i++)
	{
		if (!result) break;
		for (int j = 0; j <= 3; j++)
		{
			if (!MathUtils::IsAlmostEqualTo(m.GetElement(i, j), 1))
			{
				result = false;
				break;
			}
		}
	}
	EXPECT_TRUE(result);
}

TEST(Test_Matrix4d, Creation)
{
	XYZ origin = XYZ(0, 0, 0);
	XYZ t = XYZ(1, 1, 1);
	Matrix4d mt = Matrix4d::CreateTranslation(t);
	EXPECT_TRUE(mt.IsTranslation());
	XYZ co = mt.OfPoint(origin);
	EXPECT_TRUE(co.IsAlmostEqualTo(t));
	XYZ scale = XYZ(1, 2, 3);
	Matrix4d ms = Matrix4d::CreateScale(scale);
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(ms.GetElement(0, 0), 1) &&
				MathUtils::IsAlmostEqualTo(ms.GetElement(1, 1), 2) &&
				MathUtils::IsAlmostEqualTo(ms.GetElement(2, 2), 3));
	XYZ _scale = ms.GetScale();
	EXPECT_TRUE(scale.IsAlmostEqualTo(_scale));
	Matrix4d rotation = Matrix4d::CreateRotation(XYZ(0,0,1),Constants::Pi);
	XYZ rt = XYZ(-1, 0, 1);
	XYZ rted = rotation.OfPoint(rt);
	EXPECT_TRUE(rted.IsAlmostEqualTo(XYZ(1, 0, 1)));
	rotation = Matrix4d::CreateRotationAtPoint(XYZ(1,0,0), XYZ(0, 0, 1), Constants::Pi);
	rted = rotation.OfPoint(rt);
	EXPECT_TRUE(rted.IsAlmostEqualTo(XYZ(3, 0, 1)));
	Matrix4d reflection = Matrix4d::CreateReflection(XYZ(-1, 0, 0));
	XYZ re = XYZ(-1, 0, 1);
	XYZ red = reflection.OfPoint(re);
	EXPECT_TRUE(red.IsAlmostEqualTo(XYZ(1, 0, 1)));
}

TEST(Test_Matrix4d, View)
{

}

TEST(Test_Matrix4d, GetSet)
{
	Matrix4d m = Matrix4d();
	XYZ xyz = XYZ(5, 5, 5);
	m.SetBasisX(xyz);
	XYZ x = m.GetBasisX();
	EXPECT_TRUE(xyz.IsAlmostEqualTo(x));
	m.SetBasisY(xyz);
	XYZ y = m.GetBasisY();
	EXPECT_TRUE(xyz.IsAlmostEqualTo(y));
	m.SetBasisZ(xyz);
	XYZ z = m.GetBasisZ();
	EXPECT_TRUE(xyz.IsAlmostEqualTo(z));
	m.SetBasisW(xyz);
	XYZ w = m.GetBasisW();
	EXPECT_TRUE(xyz.IsAlmostEqualTo(w));
	m.SetElement(0, 0, 6);
	double d = m.GetElement(0, 0);
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(d, 6));
}

TEST(Test_Matrix4d, Matrix)
{
	Matrix4d m = Matrix4d();
	Matrix4d inverse;
	bool result = m.GetInverse(inverse);
	EXPECT_TRUE(result);
	EXPECT_TRUE(inverse.IsIdentity());
	Matrix4d m1 = Matrix4d(0, 0, 5, 2, 0, 0, 2, 1, 1, -2, 0, 0, 1, 1, 0, 0);
	result = m1.GetInverse(inverse);
	EXPECT_TRUE(result);
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(inverse.GetElement(0, 0), 0)   &&
				MathUtils::IsAlmostEqualTo(inverse.GetElement(0, 1), 0)   &&
				MathUtils::IsAlmostEqualTo(inverse.GetElement(0, 2), 1/3) &&
				MathUtils::IsAlmostEqualTo(inverse.GetElement(0, 3), 2/3) &&
				MathUtils::IsAlmostEqualTo(inverse.GetElement(1, 0), 0)   &&
				MathUtils::IsAlmostEqualTo(inverse.GetElement(1, 1), 0)   &&
				MathUtils::IsAlmostEqualTo(inverse.GetElement(1, 2), -1/3)&&
				MathUtils::IsAlmostEqualTo(inverse.GetElement(1, 3), 1/3) &&
				MathUtils::IsAlmostEqualTo(inverse.GetElement(2, 0), 1)   &&
				MathUtils::IsAlmostEqualTo(inverse.GetElement(2, 1), -2)  &&
				MathUtils::IsAlmostEqualTo(inverse.GetElement(2, 2), 0)   &&
				MathUtils::IsAlmostEqualTo(inverse.GetElement(2, 3), 0)   &&
				MathUtils::IsAlmostEqualTo(inverse.GetElement(3, 0), -2)  &&
				MathUtils::IsAlmostEqualTo(inverse.GetElement(3, 1), 5)   &&
				MathUtils::IsAlmostEqualTo(inverse.GetElement(3, 2), 0)   &&
				MathUtils::IsAlmostEqualTo(inverse.GetElement(3, 3), 0));
	Matrix4d transpose = m.GetTranspose();
	EXPECT_TRUE(transpose.IsIdentity());
	m = Matrix4d(4, 3, 2, 2, 0, 1, -3, 3, 0, -1, 3, 3, 0, 3, 1, 1);
	double d = m.GetDeterminant();
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(d, -240));
}

TEST(Test_Matrix4d, OperatorOverride)
{

}



