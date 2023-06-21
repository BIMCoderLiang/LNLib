/*
 * Author:
 * 2023/06/08 - Yuqing Liang (BIMCoder Liang)
 * bim.frankliang@foxmail.com
 * 微信公众号：BIMCoder梁老师
 *
 * Use of this source code is governed by a GPL-3.0 license license that can be found in
 * the LICENSE file.
 */

#include "Matrix4d.h"
#include "XYZ.h"
#include "XYZW.h"
#include "MathUtils.h"
#include <cstring>

using namespace LNLib;

LNLib::Matrix4d::Matrix4d()
{
	memset(m_matrix4d, 0, sizeof(double) * 16);
	m_matrix4d[0][0] = 1;
	m_matrix4d[1][1] = 1;
	m_matrix4d[2][2] = 1;
	m_matrix4d[3][3] = 1;
}

LNLib::Matrix4d::Matrix4d(XYZ basisX, XYZ basisY, XYZ basisZ, XYZ origin)
{
	m_matrix4d[0][0] = basisX[0];
	m_matrix4d[1][0] = basisX[1];
	m_matrix4d[2][0] = basisX[2];
	m_matrix4d[3][0] = 0;

	m_matrix4d[0][1] = basisY[0];
	m_matrix4d[1][1] = basisY[1];
	m_matrix4d[2][1] = basisY[2];
	m_matrix4d[3][1] = 0;

	m_matrix4d[0][2] = basisZ[0];
	m_matrix4d[1][2] = basisZ[1];
	m_matrix4d[2][2] = basisZ[2];
	m_matrix4d[3][2] = 0;

	m_matrix4d[0][3] = origin[0];
	m_matrix4d[1][3] = origin[1];
	m_matrix4d[2][3] = origin[2];
	m_matrix4d[3][3] = 0;
}

LNLib::Matrix4d::Matrix4d(double a00, double a01, double a02, double a03, double a10, double a11, double a12, double a13, double a20, double a21, double a22, double a23, double a30, double a31, double a32, double a33)
{
	m_matrix4d[0][0] = a00;
	m_matrix4d[0][1] = a01;
	m_matrix4d[0][2] = a02;
	m_matrix4d[0][3] = a03;

	m_matrix4d[1][0] = a10;
	m_matrix4d[1][1] = a11;
	m_matrix4d[1][2] = a12;
	m_matrix4d[1][3] = a13;

	m_matrix4d[2][0] = a20;
	m_matrix4d[2][1] = a21;
	m_matrix4d[2][2] = a22;
	m_matrix4d[2][3] = a23;

	m_matrix4d[3][0] = a30;
	m_matrix4d[3][1] = a31;
	m_matrix4d[3][2] = a32;
	m_matrix4d[3][3] = a33;
}

Matrix4d LNLib::Matrix4d::CreateTranslation(XYZ vector)
{
	Matrix4d result = Matrix4d();
	result.m_matrix4d[0][3] = vector[0];
	result.m_matrix4d[1][3] = vector[1];
	result.m_matrix4d[2][3] = vector[2];
	return result;
}

Matrix4d LNLib::Matrix4d::Multiply(const Matrix4d& right)
{
	Matrix4d result = Matrix4d();

	result.m_matrix4d[0][0] = m_matrix4d[0][0] * right.m_matrix4d[0][0] + m_matrix4d[0][1] * right.m_matrix4d[1][0] + m_matrix4d[0][2] * right.m_matrix4d[2][0] + m_matrix4d[0][3] * right.m_matrix4d[3][0];
	result.m_matrix4d[0][1] = m_matrix4d[0][0] * right.m_matrix4d[0][1] + m_matrix4d[0][1] * right.m_matrix4d[1][1] + m_matrix4d[0][2] * right.m_matrix4d[2][1] + m_matrix4d[0][3] * right.m_matrix4d[3][1];
	result.m_matrix4d[0][2] = m_matrix4d[0][0] * right.m_matrix4d[0][2] + m_matrix4d[0][1] * right.m_matrix4d[1][2] + m_matrix4d[0][2] * right.m_matrix4d[2][2] + m_matrix4d[0][3] * right.m_matrix4d[3][2];
	result.m_matrix4d[0][3] = m_matrix4d[0][0] * right.m_matrix4d[0][3] + m_matrix4d[0][1] * right.m_matrix4d[1][3] + m_matrix4d[0][2] * right.m_matrix4d[2][3] + m_matrix4d[0][3] * right.m_matrix4d[3][3];

	result.m_matrix4d[1][0] = m_matrix4d[1][0] * right.m_matrix4d[0][0] + m_matrix4d[1][1] * right.m_matrix4d[1][0] + m_matrix4d[1][2] * right.m_matrix4d[2][0] + m_matrix4d[1][3] * right.m_matrix4d[3][0];
	result.m_matrix4d[1][1] = m_matrix4d[1][0] * right.m_matrix4d[0][1] + m_matrix4d[1][1] * right.m_matrix4d[1][1] + m_matrix4d[1][2] * right.m_matrix4d[2][1] + m_matrix4d[1][3] * right.m_matrix4d[3][1];
	result.m_matrix4d[1][2] = m_matrix4d[1][0] * right.m_matrix4d[0][2] + m_matrix4d[1][1] * right.m_matrix4d[1][2] + m_matrix4d[1][2] * right.m_matrix4d[2][2] + m_matrix4d[1][3] * right.m_matrix4d[3][2];
	result.m_matrix4d[1][3] = m_matrix4d[1][0] * right.m_matrix4d[0][3] + m_matrix4d[1][1] * right.m_matrix4d[1][3] + m_matrix4d[1][2] * right.m_matrix4d[2][3] + m_matrix4d[1][3] * right.m_matrix4d[3][3];

	result.m_matrix4d[2][0] = m_matrix4d[2][0] * right.m_matrix4d[0][0] + m_matrix4d[2][1] * right.m_matrix4d[1][0] + m_matrix4d[2][2] * right.m_matrix4d[2][0] + m_matrix4d[2][3] * right.m_matrix4d[3][0];
	result.m_matrix4d[2][1] = m_matrix4d[2][0] * right.m_matrix4d[0][1] + m_matrix4d[2][1] * right.m_matrix4d[1][1] + m_matrix4d[2][2] * right.m_matrix4d[2][1] + m_matrix4d[2][3] * right.m_matrix4d[3][1];
	result.m_matrix4d[2][2] = m_matrix4d[2][0] * right.m_matrix4d[0][2] + m_matrix4d[2][1] * right.m_matrix4d[1][2] + m_matrix4d[2][2] * right.m_matrix4d[2][2] + m_matrix4d[2][3] * right.m_matrix4d[3][2];
	result.m_matrix4d[2][3] = m_matrix4d[2][0] * right.m_matrix4d[0][3] + m_matrix4d[2][1] * right.m_matrix4d[1][3] + m_matrix4d[2][2] * right.m_matrix4d[2][3] + m_matrix4d[2][3] * right.m_matrix4d[3][3];

	result.m_matrix4d[3][0] = m_matrix4d[3][0] * right.m_matrix4d[0][0] + m_matrix4d[3][1] * right.m_matrix4d[1][0] + m_matrix4d[3][2] * right.m_matrix4d[2][0] + m_matrix4d[3][3] * right.m_matrix4d[3][0];
	result.m_matrix4d[3][1] = m_matrix4d[3][0] * right.m_matrix4d[0][1] + m_matrix4d[3][1] * right.m_matrix4d[1][1] + m_matrix4d[3][2] * right.m_matrix4d[2][1] + m_matrix4d[3][3] * right.m_matrix4d[3][1];
	result.m_matrix4d[3][2] = m_matrix4d[3][0] * right.m_matrix4d[0][2] + m_matrix4d[3][1] * right.m_matrix4d[1][2] + m_matrix4d[3][2] * right.m_matrix4d[2][2] + m_matrix4d[3][3] * right.m_matrix4d[3][2];
	result.m_matrix4d[3][3] = m_matrix4d[3][0] * right.m_matrix4d[0][3] + m_matrix4d[3][1] * right.m_matrix4d[1][3] + m_matrix4d[3][2] * right.m_matrix4d[2][3] + m_matrix4d[3][3] * right.m_matrix4d[3][3];
	
	return result;
}

XYZ LNLib::Matrix4d::OfPoint(const XYZ& point)
{
	double x = m_matrix4d[0][0] * point[0] + m_matrix4d[0][1] * point[1] + m_matrix4d[0][2] * point[2] + m_matrix4d[0][3] * 1;
	double y = m_matrix4d[1][0] * point[0] + m_matrix4d[1][1] * point[1] + m_matrix4d[1][2] * point[2] + m_matrix4d[1][3] * 1;
	double z = m_matrix4d[2][0] * point[0] + m_matrix4d[2][1] * point[1] + m_matrix4d[2][2] * point[2] + m_matrix4d[2][3] * 1;
	
	double w = m_matrix4d[3][0] * point[0] + m_matrix4d[3][1] * point[1] + m_matrix4d[3][2] * point[2] + m_matrix4d[3][3] * 1;
	
	return XYZ(x,y,z)/w;
}

XYZ LNLib::Matrix4d::OfVector(const XYZ& vector)
{
	double x = m_matrix4d[0][0] * vector[0] + m_matrix4d[0][1] * vector[1] + m_matrix4d[0][2] * vector[2];
	double y = m_matrix4d[1][0] * vector[0] + m_matrix4d[1][1] * vector[1] + m_matrix4d[1][2] * vector[2];
	double z = m_matrix4d[2][0] * vector[0] + m_matrix4d[2][1] * vector[1] + m_matrix4d[2][2] * vector[2];

	return XYZ(x,y,z);
}

Matrix4d LNLib::Matrix4d::GetTranspose()
{
	return Matrix4d(m_matrix4d[0][0], m_matrix4d[1][0], m_matrix4d[2][0], m_matrix4d[3][0],
		            m_matrix4d[0][1], m_matrix4d[1][1], m_matrix4d[2][1], m_matrix4d[3][1],
		            m_matrix4d[0][2], m_matrix4d[1][2], m_matrix4d[2][2], m_matrix4d[3][2],
		            m_matrix4d[0][3], m_matrix4d[1][3], m_matrix4d[2][3], m_matrix4d[3][3]);
}

bool LNLib::Matrix4d::IsIdentity()
{
	bool c1 = MathUtils::IsAlmostEqualTo(m_matrix4d[0][0], 1);
	if (!c1)
	{
		return false;
	}
	bool c2 = MathUtils::IsAlmostEqualTo(m_matrix4d[1][1], 1);
	if (!c2)
	{
		return false;
	}
	bool c3 = MathUtils::IsAlmostEqualTo(m_matrix4d[2][2], 1);
	if (!c3)
	{
		return false;
	}

	double sum = m_matrix4d[0][1] + m_matrix4d[0][2] + m_matrix4d[0][3] +
		         m_matrix4d[1][0] + m_matrix4d[1][2] + m_matrix4d[1][3] +
		         m_matrix4d[2][0] + m_matrix4d[2][1] + m_matrix4d[2][3] +
		         m_matrix4d[3][0] + m_matrix4d[3][1] + m_matrix4d[3][2];
	bool c4 = MathUtils::IsAlmostEqualTo(sum, 0.0);
	if (!c4)
	{
		return false;
	}
	return true;
}

bool LNLib::Matrix4d::IsTranslation()
{
	bool c1 = MathUtils::IsAlmostEqualTo(m_matrix4d[0][0], 1);
	bool c2 = MathUtils::IsAlmostEqualTo(m_matrix4d[1][1], 1);
	bool c3 = MathUtils::IsAlmostEqualTo(m_matrix4d[2][2], 1);

	double sum = m_matrix4d[0][1] + m_matrix4d[0][2] +
				 m_matrix4d[1][0] + m_matrix4d[1][2] +
				 m_matrix4d[2][0] + m_matrix4d[2][1] +
				 m_matrix4d[3][0] + m_matrix4d[3][1] + m_matrix4d[3][2];

	bool c4 = MathUtils::IsAlmostEqualTo(sum, 0.0);
	return c1 && c2 && c3 && c4;
}


Matrix4d& LNLib::Matrix4d::operator=(const Matrix4d& another)
{
	std::memcpy(m_matrix4d, another.m_matrix4d, sizeof(double) * 16);
	return *this;
}
