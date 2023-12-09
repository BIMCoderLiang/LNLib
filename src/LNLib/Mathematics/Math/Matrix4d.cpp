/*
 * Author:
 * 2023/06/22 - Yuqing Liang (BIMCoder Liang)
 * bim.frankliang@foxmail.com
 * 微信公众号：BIMCoder梁老师
 *
 * Use of this source code is governed by a GPL-3.0 license that can be found in
 * the LICENSE file.
 */

#include "Matrix4d.h"
#include "XYZ.h"
#include "XYZW.h"
#include "MathUtils.h"
#include <algorithm>

using namespace LNLib;

LNLib::Matrix4d::Matrix4d()
{
	memset(m_matrix4d, 0, sizeof(double) * 16);
	for (int i = 0; i <= 3; i++)
	{
		m_matrix4d[i][i] = 1;
	}
}

LNLib::Matrix4d::Matrix4d(XYZ basisX, XYZ basisY, XYZ basisZ, XYZ basisW)
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

	m_matrix4d[0][3] = basisW[0];
	m_matrix4d[1][3] = basisW[1];
	m_matrix4d[2][3] = basisW[2];
	m_matrix4d[3][3] = 1;
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

Matrix4d LNLib::Matrix4d::CreateReflection(const XYZ& normal)
{
	Matrix4d result = Matrix4d();

	XYZ t = const_cast<XYZ&>(normal).Normalize();
	double x = t[0];
	double y = t[1];
	double z = t[2];

	result.m_matrix4d[0][0] = 1 - 2 * x * x;
	result.m_matrix4d[0][1] = -2 * x * y;
	result.m_matrix4d[0][2] = -2 * x * z;
	result.m_matrix4d[0][3] = 0;

	result.m_matrix4d[1][0] = - 2 * x * y;
	result.m_matrix4d[1][1] = 1 - 2 * y * y;
	result.m_matrix4d[1][2] = - 2 * y * z;
	result.m_matrix4d[1][3] = 0;

	result.m_matrix4d[2][0] = - 2 * x * z;
	result.m_matrix4d[2][1] = - 2 * y * z;
	result.m_matrix4d[2][2] = 1 - 2 * z * z;
	result.m_matrix4d[2][3] = 0;

	result.m_matrix4d[3][0] = 0.0;
	result.m_matrix4d[3][1] = 0.0;
	result.m_matrix4d[3][2] = 0.0;
	result.m_matrix4d[3][3] = 1.0;
	return result;
}

Matrix4d LNLib::Matrix4d::CreateRotation(const XYZ& axis, double rad)
{
	double c = cos(rad);
	double s = sin(rad);
	double t = 1.0 - c;

	XYZ nAxis = const_cast<XYZ&>(axis).Normalize();
	double x = nAxis.GetX();
	double y = nAxis.GetY();
	double z = nAxis.GetZ();

	Matrix4d result = Matrix4d();
	result.m_matrix4d[0][0] = c + x * x * t;
	result.m_matrix4d[0][1] = x * y * t - z * s;
	result.m_matrix4d[0][2] = x * z * t + y * s;
	result.m_matrix4d[0][3] = 0.0;

	result.m_matrix4d[1][0] = y * x * t + z * s;
	result.m_matrix4d[1][1] = c * y * y * t;
	result.m_matrix4d[1][2] = y * z * t - x * s;
	result.m_matrix4d[1][3] = 0.0;

	result.m_matrix4d[2][0] = z * x * t - y * s;
	result.m_matrix4d[2][1] = z * y * t + x * s;
	result.m_matrix4d[2][2] = c + z * z * t;
	result.m_matrix4d[2][3] = 0.0;

	result.m_matrix4d[3][0] = 0.0;
	result.m_matrix4d[3][1] = 0.0;
	result.m_matrix4d[3][2] = 0.0;
	result.m_matrix4d[3][3] = 1.0;
	return result;
}

Matrix4d LNLib::Matrix4d::CreateRotationAtPoint(const XYZ& origin, const XYZ& axis, double rad)
{
	Matrix4d rodrigues = CreateRotation(axis, rad);
	double x = rodrigues.GetElement(0, 0) * origin[0] + rodrigues.GetElement(0, 1) * origin[1] + rodrigues.GetElement(0, 2) * origin[2];
	double y = rodrigues.GetElement(1, 0) * origin[0] + rodrigues.GetElement(1, 1) * origin[1] + rodrigues.GetElement(1, 2) * origin[2];
	double z = rodrigues.GetElement(2, 0) * origin[0] + rodrigues.GetElement(2, 1) * origin[1] + rodrigues.GetElement(2, 2) * origin[2];
	XYZ Rv = XYZ(x,y,z);
	XYZ offset = origin - Rv;
	rodrigues.SetBasisW(offset);
	return rodrigues;
}

Matrix4d LNLib::Matrix4d::CreateTranslation(const XYZ& vector)
{
	Matrix4d result = Matrix4d();
	result.m_matrix4d[0][3] = vector[0];
	result.m_matrix4d[1][3] = vector[1];
	result.m_matrix4d[2][3] = vector[2];
	return result;
}

Matrix4d LNLib::Matrix4d::CreateScale(const XYZ& scale)
{
	Matrix4d result = Matrix4d();
	result.m_matrix4d[0][0] *= scale[0];
	result.m_matrix4d[1][1] *= scale[1];
	result.m_matrix4d[2][2] *= scale[2];
	return result;
}

void LNLib::Matrix4d::SetBasisX(const XYZ& basisX)
{
	m_matrix4d[0][0] = basisX[0];
	m_matrix4d[1][0] = basisX[1];
	m_matrix4d[2][0] = basisX[2];
}

XYZ LNLib::Matrix4d::GetBasisX() const
{
	return XYZ(m_matrix4d[0][0], m_matrix4d[1][0], m_matrix4d[2][0]);
}

void LNLib::Matrix4d::SetBasisY(const XYZ& basisY)
{
	m_matrix4d[0][1] = basisY[0];
	m_matrix4d[1][1] = basisY[1];
	m_matrix4d[2][1] = basisY[2];
}

XYZ LNLib::Matrix4d::GetBasisY() const
{
	return XYZ(m_matrix4d[0][1], m_matrix4d[1][1], m_matrix4d[2][1]);
}

void LNLib::Matrix4d::SetBasisZ(const XYZ& basisZ)
{
	m_matrix4d[0][2] = basisZ[0];
	m_matrix4d[1][2] = basisZ[1];
	m_matrix4d[2][2] = basisZ[2];
}

XYZ LNLib::Matrix4d::GetBasisZ() const
{
	return XYZ(m_matrix4d[0][2], m_matrix4d[1][2], m_matrix4d[2][2]);
}

void LNLib::Matrix4d::SetBasisW(const XYZ& basisW)
{
	m_matrix4d[0][3] = basisW[0];
	m_matrix4d[1][3] = basisW[1];
	m_matrix4d[2][3] = basisW[2];
}

XYZ LNLib::Matrix4d::GetBasisW() const
{
	return XYZ(m_matrix4d[0][3], m_matrix4d[1][3], m_matrix4d[2][3]);
}

double LNLib::Matrix4d::GetElement(int row, int column) const
{
	return m_matrix4d[row][column];
}

void LNLib::Matrix4d::SetElement(int row, int column, double value)
{
	m_matrix4d[row][column] = value;
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

XYZW LNLib::Matrix4d::OfWeightedPoint(const XYZW& point)
{
	double x = m_matrix4d[0][0] * point[0] + m_matrix4d[0][1] * point[1] + m_matrix4d[0][2] * point[2] + m_matrix4d[0][3] * point[3];
	double y = m_matrix4d[1][0] * point[0] + m_matrix4d[1][1] * point[1] + m_matrix4d[1][2] * point[2] + m_matrix4d[1][3] * point[3];
	double z = m_matrix4d[2][0] * point[0] + m_matrix4d[2][1] * point[1] + m_matrix4d[2][2] * point[2] + m_matrix4d[2][3] * point[3];

	double w = m_matrix4d[3][0] * point[0] + m_matrix4d[3][1] * point[1] + m_matrix4d[3][2] * point[2] + m_matrix4d[3][3] * point[3];

	return XYZW(XYZ(x,y,z),w);
}

XYZ LNLib::Matrix4d::OfVector(const XYZ& vector)
{
	double x = m_matrix4d[0][0] * vector[0] + m_matrix4d[0][1] * vector[1] + m_matrix4d[0][2] * vector[2];
	double y = m_matrix4d[1][0] * vector[0] + m_matrix4d[1][1] * vector[1] + m_matrix4d[1][2] * vector[2];
	double z = m_matrix4d[2][0] * vector[0] + m_matrix4d[2][1] * vector[1] + m_matrix4d[2][2] * vector[2];

	return XYZ(x,y,z);
}

bool LNLib::Matrix4d::GetInverse(Matrix4d& inverse)
{
	int n = 4;
	std::vector<std::vector<double>> currentMatrix(n);
	for (int i = 0; i < n; i++)
	{
		currentMatrix[i].resize(n);
		for (int j = 0; j < n; j++)
		{
			currentMatrix[i][j] = m_matrix4d[i][j];
		}
	}
	std::vector<std::vector<double>> inverseMatrix;
	if (!MathUtils::MakeInverse(currentMatrix, inverseMatrix))
	{
		return false;
	}
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			inverse.m_matrix4d[i][j] = inverseMatrix[i][j];
		}
	}
	return true;
}

Matrix4d LNLib::Matrix4d::GetTranspose()
{
	return Matrix4d(m_matrix4d[0][0], m_matrix4d[1][0], m_matrix4d[2][0], m_matrix4d[3][0],
		            m_matrix4d[0][1], m_matrix4d[1][1], m_matrix4d[2][1], m_matrix4d[3][1],
		            m_matrix4d[0][2], m_matrix4d[1][2], m_matrix4d[2][2], m_matrix4d[3][2],
		            m_matrix4d[0][3], m_matrix4d[1][3], m_matrix4d[2][3], m_matrix4d[3][3]);
}

XYZ LNLib::Matrix4d::GetScale()
{
	return XYZ(GetBasisX().Length(),GetBasisY().Length(),GetBasisZ().Length());
}

double LNLib::Matrix4d::GetDeterminant()
{
	std::vector<std::vector<double>> matrix(4, (std::vector<double>(4,0.0)));
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			matrix[i][j] = m_matrix4d[i][j];
		}
	}
	return MathUtils::GetDeterminant(matrix, 4);
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

bool LNLib::Matrix4d::HasReflection()
{
	XYZ cross = GetBasisX().CrossProduct(GetBasisY());
	double dot = cross.DotProduct(GetBasisZ());
	return MathUtils::IsLessThan(dot, 0.0);
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

LNLIB_EXPORT Matrix4d LNLib::operator*(const Matrix4d& left, const Matrix4d& right)
{
	Matrix4d l = left;
	Matrix4d r = right;
	return l.Multiply(r);
}

LNLIB_EXPORT Matrix4d LNLib::operator+(const Matrix4d& left, const Matrix4d& right)
{
	Matrix4d result;
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			result.SetElement(i, j, left.GetElement(i, j) + right.GetElement(i, j));
		}
	}
	return result;
}

LNLIB_EXPORT Matrix4d LNLib::operator-(const Matrix4d& left, const Matrix4d& right)
{
	Matrix4d result;
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			result.SetElement(i,j,left.GetElement(i,j)-right.GetElement(i,j));
		}
	}
	return result;
}
