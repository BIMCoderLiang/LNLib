/*
 * Author:
 * 2023/06/08 - Yuqing Liang (BIMCoder Liang)
 * bim.frankliang@foxmail.com
 * 微信公众号：BIMCoder梁老师
 *
 * Use of this source code is governed by a GPL-3.0 license license that can be found in
 * the LICENSE file.
 */

#include "XYZ.h"
#include <math.h>
#include "MathUtils.h"

using namespace LNLib;

XYZ::XYZ()
{
	m_xyz[0] = 0;
	m_xyz[1] = 0;
	m_xyz[2] = 0;
}

XYZ::XYZ(double x, double y, double z){

	m_xyz[0] = x;
	m_xyz[1] = y;
	m_xyz[2] = z;
}


void XYZ::SetX(const double x)
{
	m_xyz[0] = x;
}

double XYZ::GetX() const { return m_xyz[0]; };

void XYZ::SetY(const double y)
{
	m_xyz[1] = y;
}


double XYZ::GetY() const { return m_xyz[1]; };

void XYZ::SetZ(const double z) 
{
	m_xyz[2] = z;
}

double XYZ::GetZ() const { return m_xyz[2]; }

double XYZ::X() const
{
	return m_xyz[0];
}

double& XYZ::X()
{
	return m_xyz[0];
}

double XYZ::Y() const
{
	return m_xyz[1];
}

double& XYZ::Y()
{
	return m_xyz[1];
}

double XYZ::Z() const
{
	return m_xyz[2];
}

double& XYZ::Z()
{
	return m_xyz[2];
}

bool LNLib::XYZ::IsZero(const double epsilon) const
{
	return SqrLength() <= epsilon * epsilon;
}

bool LNLib::XYZ::IsUnit(const double epsilon) const
{
	return fabs(SqrLength() - 1) < epsilon * epsilon;
}

bool LNLib::XYZ::IsAlmostEqualTo(const XYZ& another) const
{
	return (*this - another).SqrLength() <= Constants::DoubleEpsilon * Constants::DoubleEpsilon;
}

double XYZ::Length() const
{
	return sqrt(m_xyz[0]* m_xyz[0] + m_xyz[1] * m_xyz[1] + m_xyz[2] * m_xyz[2]);
}

double XYZ::SqrLength() const
{
	return m_xyz[0] * m_xyz[0] + m_xyz[1] * m_xyz[1] + m_xyz[2] * m_xyz[2];
}

double LNLib::XYZ::AngleTo(const XYZ& another) const
{
	XYZ tTemp = (*this);
	XYZ ntTemp = tTemp.Normalize();

	XYZ aTemp = another;
	XYZ naTemp = aTemp.Normalize();

	if (ntTemp.IsZero() ||
		naTemp.IsZero())
	{
		return 0.0;
	}
	double dot = ntTemp.DotProduct(naTemp);
	dot = dot < -1.0 ? -1.0 : dot > 1.0 ? 1.0 : dot;
	return acos(dot);
}

XYZ XYZ::Normalize()
{
	double length = Length();
	if (length > 0)
	{
		double invLength = (double)(1.0 / length);
		m_xyz[0] *= invLength;
		m_xyz[1] *= invLength;
		m_xyz[2] *= invLength;
	}
	return *this;
}

XYZ XYZ::Add(const XYZ& another) const
{
	return *this + another;
}

XYZ XYZ::Substract(const XYZ& another) const
{
	return *this - another;
}

XYZ XYZ::Negative() const
{
	return XYZ(-m_xyz[0], -m_xyz[1], -m_xyz[2]);
}

double XYZ::DotProduct(const XYZ& another) const
{
	return m_xyz[0] * another.m_xyz[0] + m_xyz[1] * another.m_xyz[1] + m_xyz[2] * another.m_xyz[2];
}

XYZ XYZ::CrossProduct(const XYZ& another) const
{
	return XYZ(m_xyz[1] * another.m_xyz[2] - another.m_xyz[1] * m_xyz[2],
				m_xyz[2] * another.m_xyz[0] - another.m_xyz[2] * m_xyz[0],
				m_xyz[0] * another[1] - another.m_xyz[0] * m_xyz[1]);
}

double LNLib::XYZ::Distance(const XYZ& another) const
{
	double squareValue = pow((another.GetX() - m_xyz[0]),2) + pow((another.GetY() - m_xyz[1]),2) + pow((another.GetZ() - m_xyz[2]),2);
	return sqrt(squareValue);
}

XYZ& XYZ::operator=(const XYZ& xyz)
{
	m_xyz[0] = xyz.m_xyz[0];
	m_xyz[1] = xyz.m_xyz[1];
	m_xyz[2] = xyz.m_xyz[2];
	return *this;
}

double& XYZ::operator[](int index)
{
	return m_xyz[index];
}

const double& XYZ::operator[](int index) const
{
	return m_xyz[index];
}

XYZ XYZ::operator+(const XYZ& xyz) const
{
	return XYZ(m_xyz[0] + xyz.m_xyz[0], m_xyz[1] + xyz.m_xyz[1], m_xyz[2] + xyz.m_xyz[2]);
}

XYZ XYZ::operator-(const XYZ& xyz) const
{
	return XYZ(m_xyz[0] - xyz.m_xyz[0], m_xyz[1] - xyz.m_xyz[1], m_xyz[2] - xyz.m_xyz[2]);
}

double XYZ::operator*(const XYZ& xyz) const
{
	return DotProduct(xyz);
}

XYZ& XYZ::operator*=(const double& d)
{
	m_xyz[0] *= d;
	m_xyz[1] *= d;
	m_xyz[2] *= d;
	return *this;
}

XYZ& XYZ::operator/=(const double& d)
{
	m_xyz[0] /= d;
	m_xyz[1] /= d;
	m_xyz[2] /= d;
	return *this;
}

XYZ& XYZ::operator+=(const XYZ& xyz)
{
	m_xyz[0] += xyz.m_xyz[0];
	m_xyz[1] += xyz.m_xyz[1];
	m_xyz[2] += xyz.m_xyz[2];
	return *this;
}

XYZ& XYZ::operator-=(const XYZ& xyz)
{
	m_xyz[0] -= xyz.m_xyz[0];
	m_xyz[1] -= xyz.m_xyz[1];
	m_xyz[2] -= xyz.m_xyz[2];
	return *this;
}

XYZ XYZ::operator-() const
{
	return XYZ(-m_xyz[0],-m_xyz[1],-m_xyz[2]);
}
XYZ LNLib::operator*(const XYZ& source, const double d)
{
	return XYZ(source.GetX()*d, source.GetY()*d, source.GetZ()*d);
}
XYZ LNLib::operator*(const double& d, const XYZ& source)
{
	return XYZ(source.GetX() * d, source.GetY() * d, source.GetZ() * d);
}
XYZ LNLib::operator^(const XYZ& xyz1, const XYZ& xyz2)
{
	return xyz1.CrossProduct(xyz2);
}
XYZ LNLib::operator/(const XYZ& source, double d)
{
	return XYZ(source.GetX() / d, source.GetY() / d, source.GetZ() / d);
}
