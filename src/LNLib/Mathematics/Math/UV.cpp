/*
 * Author:
 * 2023/06/08 - Yuqing Liang (BIMCoder Liang)
 * bim.frankliang@foxmail.com
 * 
 *
 * Use of this source code is governed by a GPL-3.0 license that can be found in
 * the LICENSE file.
 */

#include "UV.h"
#include "MathUtils.h"

#include <cmath>

using namespace LNLib;

LNLib::UV::UV()
{
	m_uv[0] = 0;
	m_uv[1] = 0;
}


UV::UV(double u, double v) {

	m_uv[0] = u;
	m_uv[1] = v;
}


void UV::SetU(const double x)
{
	m_uv[0] = x;
}

double UV::GetU() const { return m_uv[0]; };

void UV::SetV(const double y)
{
	m_uv[1] = y;
}


double UV::GetV() const { return m_uv[1]; };


double UV::U() const
{
	return m_uv[0];
}

double& UV::U()
{
	return m_uv[0];
}

double UV::V() const
{
	return m_uv[1];
}

double& UV::V()
{
	return m_uv[1];
}


bool LNLib::UV::IsZero(const double epsilon) const
{
	return SqrLength() <= epsilon * epsilon;
}

bool LNLib::UV::IsUnit(const double epsilon) const
{
	return std::fabs(SqrLength() - 1) < epsilon * epsilon;
}

bool LNLib::UV::IsAlmostEqualTo(const UV& another) const
{
	return MathUtils::IsAlmostEqualTo(m_uv[0], another.m_uv[0]) &&
		   MathUtils::IsAlmostEqualTo(m_uv[1], another.m_uv[1]);
}

double UV::Length() const
{
	return std::sqrt(m_uv[0] * m_uv[0] + m_uv[1] * m_uv[1]);
}

double UV::SqrLength() const
{
	return m_uv[0] * m_uv[0] + m_uv[1] * m_uv[1];
}

double LNLib::UV::AngleTo(const UV& another) const
{
	return 0.0;
}

UV UV::Normalize()
{
	double length = Length();
	UV newUV = *this;
	if (length > 0)
	{
		double invLength = (double)(1.0 / length);
		newUV.m_uv[0] *= invLength;
		newUV.m_uv[1] *= invLength;
	}
	return newUV;
}

UV UV::Add(const UV& another) const
{
	return *this + another;
}

UV UV::Substract(const UV& another) const
{
	return *this - another;
}

UV UV::Negative() const
{
	return UV(-m_uv[0], -m_uv[1]);
}

double UV::DotProduct(const UV& another) const
{
	return m_uv[0] * another.m_uv[0] + m_uv[1] * another.m_uv[1];
}

double UV::CrossProduct(const UV& another) const
{
	return m_uv[0] * another[1] - another.m_uv[0] * m_uv[1];
}

double LNLib::UV::Distance(const UV& another) const
{
	double squareValue = pow((another.GetU() - m_uv[0]), 2) + pow((another.GetV() - m_uv[1]), 2);
	return std::sqrt(squareValue);
}

UV& UV::operator=(const UV& uv)
{
	m_uv[0] = uv.m_uv[0];
	m_uv[1] = uv.m_uv[1];
	return *this;
}

double& UV::operator[](int index)
{
	return m_uv[index];
}

const double& UV::operator[](int index) const
{
	return m_uv[index];
}

UV UV::operator+(const UV& uv) const
{
	return LNLib::UV(m_uv[0] + uv.m_uv[0], m_uv[1] + uv.m_uv[1]);
}

UV UV::operator-(const UV& uv) const
{
	return LNLib::UV(m_uv[0] - uv.m_uv[0], m_uv[1] - uv.m_uv[1]);
}

double UV::operator*(const UV& uv) const
{
	return DotProduct(uv);
}

UV& UV::operator*=(const double& d)
{
	m_uv[0] *= d;
	m_uv[1] *= d;
	return *this;
}

UV& UV::operator/=(const double& d)
{
	m_uv[0] /= d;
	m_uv[1] /= d;
	return *this;
}

UV& UV::operator+=(const UV& uv)
{
	m_uv[0] += uv.m_uv[0];
	m_uv[1] += uv.m_uv[1];
	return *this;
}

UV& UV::operator-=(const UV& uv)
{
	m_uv[0] -= uv.m_uv[0];
	m_uv[1] -= uv.m_uv[1];
	return *this;
}

UV UV::operator-() const
{
	return UV(-m_uv[0], -m_uv[1]);
}

LNLib::UV LNLib::operator*(const UV& source, const double d)
{
	return UV(source.GetU() * d, source.GetV() * d);
}
LNLib::UV LNLib::operator*(const double& d, const UV& source)
{
	return UV(source.GetU() * d, source.GetV() * d);
}
double LNLib::operator^(const UV& uv1, const UV& uv2)
{
	return uv1.CrossProduct(uv2);
}
LNLib::UV LNLib::operator/(const UV& source, double d)
{
	return UV(source.GetU() / d, source.GetV() / d);
}
