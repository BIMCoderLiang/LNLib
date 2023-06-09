/*
 * Author:
 * 2023/06/08 - Yuqing Liang (BIMCoder Liang)
 * bim.frankliang@foxmail.com
 * 微信公众号：BIMCoder梁老师
 *
 * Use of this source code is governed by a GPL-3.0 license license that can be found in
 * the LICENSE file.
 */

#include "UV.h"
#include <math.h>

using namespace LNLib;

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
	return fabs(SqrLength() - 1) < epsilon * epsilon;
}

bool LNLib::UV::IsAlmostEqualTo(const UV& another) const
{
	return (*this - another).SqrLength() <= Constants::DoubleEpsilon * Constants::DoubleEpsilon;
}

double UV::Length() const
{
	return sqrt(m_uv[0] * m_uv[0] + m_uv[1] * m_uv[1]);
}

double UV::SqrLength() const
{
	return m_uv[0] * m_uv[0] + m_uv[1] * m_uv[1];
}

double UV::Normalize()
{
	double length = Length();
	if (length > 0)
	{
		double invLength = (double)(1.0 / length);
		m_uv[0] *= invLength;
		m_uv[1] *= invLength;
	}
	return length;
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

UV& UV::operator=(const UV& UV)
{
	m_uv[0] = UV.m_uv[0];
	m_uv[1] = UV.m_uv[1];
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

UV UV::operator+(const UV& UV) const
{
	return LNLib::UV(m_uv[0] + UV.m_uv[0], m_uv[1] + UV.m_uv[1]);
}

UV UV::operator-(const UV& UV) const
{
	return LNLib::UV(m_uv[0] - UV.m_uv[0], m_uv[1] - UV.m_uv[1]);
}

double UV::operator*(const UV& UV) const
{
	return DotProduct(UV);
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

UV& UV::operator+=(const UV& UV)
{
	m_uv[0] += UV.m_uv[0];
	m_uv[1] += UV.m_uv[1];
	return *this;
}

UV& UV::operator-=(const UV& UV)
{
	m_uv[0] -= UV.m_uv[0];
	m_uv[1] -= UV.m_uv[1];
	return *this;
}

UV UV::operator-() const
{
	return UV(-m_uv[0], -m_uv[1]);
}

UV LNLib::operator*(const UV& source, const double d)
{
	return UV(source.GetU() * d, source.GetV() * d);
}
UV LNLib::operator*(const double& d, const UV& source)
{
	return UV(source.GetU() * d, source.GetV() * d);
}
double LNLib::operator^(const UV& UV1, const UV& UV2)
{
	return UV1.CrossProduct(UV2);
}
UV LNLib::operator/(const UV& source, double d)
{
	return UV(source.GetU() / d, source.GetV() / d);
}
