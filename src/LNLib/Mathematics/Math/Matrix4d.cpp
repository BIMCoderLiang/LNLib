/*
 * Author:
 * 2023/06/22 - Yuqing Liang (BIMCoder Liang)
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

Matrix4d LNLib::Matrix4d::CreateReflection(const XYZ& origin, const XYZ& normal)
{
	Matrix4d result = Matrix4d();

	XYZ temp = normal;
	XYZ nTemp = temp.Normalize();
	XYZW wNormal = XYZW(nTemp.GetX(), nTemp.GetY(), nTemp.GetZ(), -origin.DotProduct(nTemp));

	result.m_matrix4d[0][0] = 1 - 2 * wNormal[0] * wNormal[0];
	result.m_matrix4d[0][1] = - 2 * wNormal[1] * wNormal[0];
	result.m_matrix4d[0][2] = - 2 * wNormal[2] * wNormal[0];
	result.m_matrix4d[0][3] = - 2 * wNormal[3] * wNormal[0];

	result.m_matrix4d[0][0] = - 2 * wNormal[0] * wNormal[1];
	result.m_matrix4d[0][0] = 1 - 2 * wNormal[1] * wNormal[1];
	result.m_matrix4d[0][0] = - 2 * wNormal[2] * wNormal[1];
	result.m_matrix4d[0][0] = - 2 * wNormal[3] * wNormal[1];

	result.m_matrix4d[0][0] = - 2 * wNormal[0] * wNormal[2];
	result.m_matrix4d[0][0] = - 2 * wNormal[1] * wNormal[2];
	result.m_matrix4d[0][0] = 1 - 2 * wNormal[2] * wNormal[2];
	result.m_matrix4d[0][0] = - 2 * wNormal[3] * wNormal[0];

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

	XYZ tempAxis = axis;
	XYZ nAxis = tempAxis.Normalize();
	double x = nAxis.GetX();
	double y = nAxis.GetY();
	double z = nAxis.GetZ();

	Matrix4d result = Matrix4d();
	result.m_matrix4d[0][0] = 1 + t * (x * x - 1);
	result.m_matrix4d[0][1] = z * s + t * x * y;
	result.m_matrix4d[0][2] = -y * s + t * x * z;
	result.m_matrix4d[0][3] = 0.0;

	result.m_matrix4d[1][0] = -z * s + t * x * y;
	result.m_matrix4d[1][1] = 1 + t * (y * y - 1);
	result.m_matrix4d[1][2] = x * s + t * y * z;
	result.m_matrix4d[1][3] = 0.0;

	result.m_matrix4d[2][0] = y * s + t * x * z;
	result.m_matrix4d[2][1] = -x * s + t * y * z;
	result.m_matrix4d[2][2] = 1 + t * (z * z - 1);
	result.m_matrix4d[2][3] = 0.0;

	result.m_matrix4d[3][0] = 0.0;
	result.m_matrix4d[3][1] = 0.0;
	result.m_matrix4d[3][2] = 0.0;
	result.m_matrix4d[3][3] = 1.0;
	return result;
}

Matrix4d LNLib::Matrix4d::CreateTranslation(const XYZ& vector)
{
	Matrix4d result = Matrix4d();
	result.m_matrix4d[0][3] = vector[0];
	result.m_matrix4d[1][3] = vector[1];
	result.m_matrix4d[2][3] = vector[2];
	return result;
}

Matrix4d LNLib::Matrix4d::CreateScale(double scale, bool isScaleOrigin)
{
	Matrix4d result = Matrix4d();
	result.m_matrix4d[0][0] *= scale;
	result.m_matrix4d[1][1] *= scale;
	result.m_matrix4d[2][2] *= scale;

	if (isScaleOrigin)
	{
		result.m_matrix4d[0][3] *= scale;
		result.m_matrix4d[1][3] *= scale;
		result.m_matrix4d[2][3] *= scale;
	}
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

Matrix4d LNLib::Matrix4d::CreateCamera(const XYZ& eyePoint, const XYZ& lookDirection, const XYZ& upDirection, const XYZ& rightDirection)
{
	XYZ lTemp = lookDirection;
	XYZ look = lTemp.Normalize();

	XYZ uTemp = upDirection;
	XYZ up = uTemp.Normalize();

	XYZ rTemp = upDirection;
	XYZ right = rTemp.Normalize();

	Matrix4d result = Matrix4d();
	result.m_matrix4d[0][0] = right.GetX();
	result.m_matrix4d[1][0] = right.GetY();
	result.m_matrix4d[2][0] = right.GetZ();
	result.m_matrix4d[3][0] = -right.DotProduct(eyePoint);

	result.m_matrix4d[0][1] = up.GetX();
	result.m_matrix4d[1][1] = up.GetY();
	result.m_matrix4d[2][1] = up.GetZ();
	result.m_matrix4d[3][1] = -up.DotProduct(eyePoint);

	result.m_matrix4d[0][2] = look.GetX();
	result.m_matrix4d[1][2] = look.GetY();
	result.m_matrix4d[2][2] = look.GetZ();
	result.m_matrix4d[3][2] = -look.DotProduct(eyePoint);

	result.m_matrix4d[0][3] = 0.0;
	result.m_matrix4d[1][3] = 0.0;
	result.m_matrix4d[2][3] = 0.0;
	result.m_matrix4d[3][3] = 1.0;
	return result;
}

Matrix4d LNLib::Matrix4d::CreateLookAt(const XYZ& eyePoint, const XYZ& targetPoint, const XYZ& upDirection)
{
	XYZ xAxis, yAxis, zAxis;
	zAxis = (eyePoint - targetPoint).Normalize();
	xAxis = (upDirection.CrossProduct(zAxis)).Normalize();
	yAxis = (zAxis.CrossProduct(xAxis)).Normalize();

	Matrix4d result = Matrix4d();
	result.m_matrix4d[0][0] = xAxis.GetX();
	result.m_matrix4d[1][0] = xAxis.GetY();
	result.m_matrix4d[2][0] = xAxis.GetZ();
	result.m_matrix4d[3][0] = -xAxis.DotProduct(eyePoint);

	result.m_matrix4d[0][1] = yAxis.GetX();
	result.m_matrix4d[1][1] = yAxis.GetY();
	result.m_matrix4d[2][1] = yAxis.GetZ();
	result.m_matrix4d[3][1] = -yAxis.DotProduct(eyePoint);

	result.m_matrix4d[0][2] = zAxis.GetX();
	result.m_matrix4d[1][2] = zAxis.GetY();
	result.m_matrix4d[2][2] = zAxis.GetZ();
	result.m_matrix4d[3][2] = -zAxis.DotProduct(eyePoint);

	result.m_matrix4d[0][3] = 0.0;
	result.m_matrix4d[1][3] = 0.0;
	result.m_matrix4d[2][3] = 0.0;
	result.m_matrix4d[3][3] = 1.0;
	return result;
}

Matrix4d LNLib::Matrix4d::Orthogonal(double width, double height, double zNear, double zFar)
{
	Matrix4d result = Matrix4d();
	result.m_matrix4d[0][0] = 2.0 / width;
	result.m_matrix4d[1][0] = 0.0;
	result.m_matrix4d[2][0] = 0.0;
	result.m_matrix4d[3][0] = 0.0;

	result.m_matrix4d[0][1] = 0.0;
	result.m_matrix4d[1][1] = 2.0 / height;
	result.m_matrix4d[2][1] = 0.0;
	result.m_matrix4d[3][1] = 0.0;

	result.m_matrix4d[0][2] = 0.0;
	result.m_matrix4d[1][2] = 0.0;
	result.m_matrix4d[2][2] = 1.0 / (zNear - zFar);
	result.m_matrix4d[3][2] = zNear / (zNear - zFar);

	result.m_matrix4d[0][3] = 0.0;
	result.m_matrix4d[1][3] = 0.0;
	result.m_matrix4d[2][3] = 0.0;
	result.m_matrix4d[3][3] = 1.0;
	return result;
}

Matrix4d LNLib::Matrix4d::Perspective(double width, double height, double zNear, double zFar)
{
	Matrix4d result = Matrix4d();
	result.m_matrix4d[0][0] = 2.0 * zNear / width;
	result.m_matrix4d[1][0] = 0.0;
	result.m_matrix4d[2][0] = 0.0;
	result.m_matrix4d[3][0] = 0.0;

	result.m_matrix4d[0][1] = 0.0;
	result.m_matrix4d[1][1] = 2.0 * zNear / height;
	result.m_matrix4d[2][1] = 0.0;
	result.m_matrix4d[3][1] = 0.0;

	result.m_matrix4d[0][2] = 0.0;
	result.m_matrix4d[1][2] = 0.0;
	result.m_matrix4d[2][2] = zFar / (zNear - zFar);
	result.m_matrix4d[3][2] = zFar * zNear / (zNear - zFar);

	result.m_matrix4d[0][3] = 0.0;
	result.m_matrix4d[1][3] = 0.0;
	result.m_matrix4d[2][3] = -1.0;
	result.m_matrix4d[3][3] = 0.0;
	return result;
}

Matrix4d LNLib::Matrix4d::PerspectiveFov(double fov, double aspect, double zNear, double zFar)
{
	double width = 1.0f / tan(fov / 2.0);
    double height = aspect / tan(fov / 2.0);

	Matrix4d result = Matrix4d();
	result.m_matrix4d[0][0] = width;
	result.m_matrix4d[1][0] = 0.0;
	result.m_matrix4d[2][0] = 0.0;
	result.m_matrix4d[3][0] = 0.0;

	result.m_matrix4d[0][1] = 0.0;
	result.m_matrix4d[1][1] = height;
	result.m_matrix4d[2][1] = 0.0;
	result.m_matrix4d[3][1] = 0.0;

	result.m_matrix4d[0][2] = 0.0;
	result.m_matrix4d[1][2] = 0.0;
	result.m_matrix4d[2][2] = zFar / (zNear - zFar);
	result.m_matrix4d[3][2] = zFar * zNear / (zNear - zFar);

	result.m_matrix4d[0][3] = 0.0;
	result.m_matrix4d[1][3] = 0.0;
	result.m_matrix4d[2][3] = -1.0;
	result.m_matrix4d[3][3] = 0.0;
	return result;
}

Matrix4d LNLib::Matrix4d::PerspectiveMultiFovs(double fovX, double fovY, double zNear, double zFar)
{
	double width = 1.0 / tan(fovX / 2.0);
	double height = 1.0f / tan(fovY / 2.0);

	Matrix4d result = Matrix4d();
	result.m_matrix4d[0][0] = width;
	result.m_matrix4d[1][0] = 0.0;
	result.m_matrix4d[2][0] = 0.0;
	result.m_matrix4d[3][0] = 0.0;

	result.m_matrix4d[0][1] = 0.0;
	result.m_matrix4d[1][1] = height;
	result.m_matrix4d[2][1] = 0.0;
	result.m_matrix4d[3][1] = 0.0;

	result.m_matrix4d[0][2] = 0.0;
	result.m_matrix4d[1][2] = 0.0;
	result.m_matrix4d[2][2] = zFar / (zNear - zFar);
	result.m_matrix4d[3][2] = zFar * zNear / (zNear - zFar);

	result.m_matrix4d[0][3] = 0.0;
	result.m_matrix4d[1][3] = 0.0;
	result.m_matrix4d[2][3] = -1.0;
	result.m_matrix4d[3][3] = 0.0;
	return result;
}

void LNLib::Matrix4d::SetBasisX(const XYZ basisX)
{
	m_matrix4d[0][0] = basisX[0];
	m_matrix4d[1][0] = basisX[1];
	m_matrix4d[2][0] = basisX[2];
}

XYZ LNLib::Matrix4d::GetBasisX() const
{
	return XYZ(m_matrix4d[0][0], m_matrix4d[1][0], m_matrix4d[2][0]);
}

void LNLib::Matrix4d::SetBasisY(const XYZ basisY)
{
	m_matrix4d[0][1] = basisY[0];
	m_matrix4d[1][1] = basisY[1];
	m_matrix4d[2][1] = basisY[2];
}

XYZ LNLib::Matrix4d::GetBasisY() const
{
	return XYZ(m_matrix4d[0][1], m_matrix4d[1][1], m_matrix4d[2][1]);
}

void LNLib::Matrix4d::SetBasisZ(const XYZ basisZ)
{
	m_matrix4d[0][2] = basisZ[0];
	m_matrix4d[1][2] = basisZ[1];
	m_matrix4d[2][2] = basisZ[2];
}

XYZ LNLib::Matrix4d::GetBasisZ() const
{
	return XYZ(m_matrix4d[0][2], m_matrix4d[1][2], m_matrix4d[2][2]);
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

bool LNLib::Matrix4d::GetInverse(Matrix4d& inverse)
{
	double tempResult[4][4];
	double tmp[12];
	double src[16]; 
	double det;
	for (int i = 0; i < 4; i++)
	{
		src[i + 0] = (*this).m_matrix4d[i][0];
		src[i + 4] = (*this).m_matrix4d[i][1];
		src[i + 8] = (*this).m_matrix4d[i][2];
		src[i + 12] = (*this).m_matrix4d[i][3];
	}
	tmp[0] = src[10] * src[15];
	tmp[1] = src[11] * src[14];
	tmp[2] = src[9] * src[15];
	tmp[3] = src[11] * src[13];
	tmp[4] = src[9] * src[14];
	tmp[5] = src[10] * src[13];
	tmp[6] = src[8] * src[15];
	tmp[7] = src[11] * src[12];
	tmp[8] = src[8] * src[14];
	tmp[9] = src[10] * src[12];
	tmp[10] = src[8] * src[13];
	tmp[11] = src[9] * src[12];

	tempResult[0][0] = tmp[0] * src[5] + tmp[3] * src[6] + tmp[4] * src[7];
	tempResult[0][0] -= tmp[1] * src[5] + tmp[2] * src[6] + tmp[5] * src[7];
	tempResult[0][1] = tmp[1] * src[4] + tmp[6] * src[6] + tmp[9] * src[7];
	tempResult[0][1] -= tmp[0] * src[4] + tmp[7] * src[6] + tmp[8] * src[7];
	tempResult[0][2] = tmp[2] * src[4] + tmp[7] * src[5] + tmp[10] * src[7];
	tempResult[0][2] -= tmp[3] * src[4] + tmp[6] * src[5] + tmp[11] * src[7];
	tempResult[0][3] = tmp[5] * src[4] + tmp[8] * src[5] + tmp[11] * src[6];
	tempResult[0][3] -= tmp[4] * src[4] + tmp[9] * src[5] + tmp[10] * src[6];
	tempResult[1][0] = tmp[1] * src[1] + tmp[2] * src[2] + tmp[5] * src[3];
	tempResult[1][0] -= tmp[0] * src[1] + tmp[3] * src[2] + tmp[4] * src[3];
	tempResult[1][1] = tmp[0] * src[0] + tmp[7] * src[2] + tmp[8] * src[3];
	tempResult[1][1] -= tmp[1] * src[0] + tmp[6] * src[2] + tmp[9] * src[3];
	tempResult[1][2] = tmp[3] * src[0] + tmp[6] * src[1] + tmp[11] * src[3];
	tempResult[1][2] -= tmp[2] * src[0] + tmp[7] * src[1] + tmp[10] * src[3];
	tempResult[1][3] = tmp[4] * src[0] + tmp[9] * src[1] + tmp[10] * src[2];
	tempResult[1][3] -= tmp[5] * src[0] + tmp[8] * src[1] + tmp[11] * src[2];

	tmp[0] = src[2] * src[7];
	tmp[1] = src[3] * src[6];
	tmp[2] = src[1] * src[7];
	tmp[3] = src[3] * src[5];
	tmp[4] = src[1] * src[6];
	tmp[5] = src[2] * src[5];

	tmp[6] = src[0] * src[7];
	tmp[7] = src[3] * src[4];
	tmp[8] = src[0] * src[6];
	tmp[9] = src[2] * src[4];
	tmp[10] = src[0] * src[5];
	tmp[11] = src[1] * src[4];

	tempResult[2][0] = tmp[0] * src[13] + tmp[3] * src[14] + tmp[4] * src[15];
	tempResult[2][0] -= tmp[1] * src[13] + tmp[2] * src[14] + tmp[5] * src[15];
	tempResult[2][1] = tmp[1] * src[12] + tmp[6] * src[14] + tmp[9] * src[15];
	tempResult[2][1] -= tmp[0] * src[12] + tmp[7] * src[14] + tmp[8] * src[15];
	tempResult[2][2] = tmp[2] * src[12] + tmp[7] * src[13] + tmp[10] * src[15];
	tempResult[2][2] -= tmp[3] * src[12] + tmp[6] * src[13] + tmp[11] * src[15];
	tempResult[2][3] = tmp[5] * src[12] + tmp[8] * src[13] + tmp[11] * src[14];
	tempResult[2][3] -= tmp[4] * src[12] + tmp[9] * src[13] + tmp[10] * src[14];
	tempResult[3][0] = tmp[2] * src[10] + tmp[5] * src[11] + tmp[1] * src[9];
	tempResult[3][0] -= tmp[4] * src[11] + tmp[0] * src[9] + tmp[3] * src[10];
	tempResult[3][1] = tmp[8] * src[11] + tmp[0] * src[8] + tmp[7] * src[10];
	tempResult[3][1] -= tmp[6] * src[10] + tmp[9] * src[11] + tmp[1] * src[8];
	tempResult[3][2] = tmp[6] * src[9] + tmp[11] * src[11] + tmp[3] * src[8];
	tempResult[3][2] -= tmp[10] * src[11] + tmp[2] * src[8] + tmp[7] * src[9];
	tempResult[3][3] = tmp[10] * src[10] + tmp[4] * src[8] + tmp[9] * src[9];
	tempResult[3][3] -= tmp[8] * src[9] + tmp[11] * src[10] + tmp[5] * src[8];

	det = src[0] * tempResult[0][0] + src[1] * tempResult[0][1] + src[2] * tempResult[0][2] + src[3] * tempResult[0][3];
	if (MathUtils::IsAlmostEqualTo(det, 0.0))
	{
		return false;
	}
	det = 1.0 / det;

	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			inverse.m_matrix4d[i][j] = tempResult[i][j] * det;
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
	return m_matrix4d[0][0] * m_matrix4d[1][1] * m_matrix4d[2][2] * m_matrix4d[3][3] +
		   m_matrix4d[0][1] * m_matrix4d[1][2] * m_matrix4d[2][3] * m_matrix4d[3][0] +
		   m_matrix4d[0][2] * m_matrix4d[1][3] * m_matrix4d[2][0] * m_matrix4d[3][1] +
		   m_matrix4d[0][3] * m_matrix4d[1][0] * m_matrix4d[2][1] * m_matrix4d[3][2] -
		   m_matrix4d[0][3] * m_matrix4d[1][1] * m_matrix4d[2][2] * m_matrix4d[3][0] -
		   m_matrix4d[0][2] * m_matrix4d[1][3] * m_matrix4d[2][0] * m_matrix4d[3][1] -
		   m_matrix4d[0][1] * m_matrix4d[1][0] * m_matrix4d[2][3] * m_matrix4d[3][2] -
		   m_matrix4d[0][0] * m_matrix4d[1][2] * m_matrix4d[2][1] * m_matrix4d[3][3];
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

	if (MathUtils::IsLessThan(dot, 0.0))
	{
		return true;
	}
	return false;
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
