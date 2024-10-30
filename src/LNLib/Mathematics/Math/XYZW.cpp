/*
 * Author:
 * 2023/06/08 - Yuqing Liang (BIMCoder Liang)
 * bim.frankliang@foxmail.com
 * 
 *
 * Use of this source code is governed by a GPL-3.0 license that can be found in
 * the LICENSE file.
 */

#include "XYZW.h"
#include "XYZ.h"
#include "MathUtils.h"

#include <cmath>

using namespace LNLib;

LNLib::XYZW::XYZW()
{
	m_xyzw[0] = 0;
	m_xyzw[1] = 0;
	m_xyzw[2] = 0;
	m_xyzw[3] = 0;
}

LNLib::XYZW::XYZW(XYZ xyz, double w)
{
	m_xyzw[0] = xyz.GetX() * w;
	m_xyzw[1] = xyz.GetY() * w;
	m_xyzw[2] = xyz.GetZ() * w;
	m_xyzw[3] = w;
}

LNLib::XYZW::XYZW(double wx, double wy, double wz, double w)
{
	m_xyzw[0] = wx;
	m_xyzw[1] = wy;
	m_xyzw[2] = wz;
	m_xyzw[3] = w;
}

double XYZW::GetWX() const { return m_xyzw[0]; };

double XYZW::GetWY() const { return m_xyzw[1]; };

double XYZW::GetWZ() const { return m_xyzw[2]; }

void XYZW::SetW(const double w)
{
	XYZ origin = this->ToXYZ(true);
	*this = XYZW(origin, w);
}

double XYZW::GetW() const { return m_xyzw[3]; }

double XYZW::WX() const
{
	return m_xyzw[0];
}

double& XYZW::WX()
{
	return m_xyzw[0];
}

double XYZW::WY() const
{
	return m_xyzw[1];
}

double& XYZW::WY()
{
	return m_xyzw[1];
}

double XYZW::WZ() const
{
	return m_xyzw[2];
}

double& XYZW::WZ()
{
	return m_xyzw[2];
}

double XYZW::W() const
{
	return m_xyzw[3];
}

double& XYZW::W()
{
	return m_xyzw[3];
}

XYZ LNLib::XYZW::ToXYZ(bool divideWeight)const
{
	if (divideWeight)
	{
		double w = m_xyzw[3];
		if (MathUtils::IsAlmostEqualTo(w, 0.0))
		{
			return XYZ(m_xyzw[0], m_xyzw[1], m_xyzw[2]);
		}
		else
		{
			return XYZ(m_xyzw[0] / w, m_xyzw[1] / w, m_xyzw[2] / w);
		}
	}
	else
	{
		return XYZ(m_xyzw[0], m_xyzw[1], m_xyzw[2]);
	}
}

bool LNLib::XYZW::IsAlmostEqualTo(const XYZW& another) const
{
	XYZW self = *this;
	XYZW temp = another;
	return  self.ToXYZ(true).IsAlmostEqualTo(temp.ToXYZ(true));
}

double LNLib::XYZW::Distance(const XYZW& another) const
{
	double squareValue = pow((another.GetWX() - m_xyzw[0]), 2) + pow((another.GetWY() - m_xyzw[1]), 2) + pow((another.GetWZ() - m_xyzw[2]), 2) + pow((another.GetW() - m_xyzw[3]), 2);
	return std::sqrt(squareValue);
}

double& LNLib::XYZW::operator[](int index)
{
	return m_xyzw[index];
}

const double& LNLib::XYZW::operator[](int index) const
{
	return m_xyzw[index];
}

XYZW LNLib::XYZW::operator+(const XYZW& xyzw) const
{
	return XYZW(m_xyzw[0] + xyzw.m_xyzw[0], m_xyzw[1] + xyzw.m_xyzw[1], m_xyzw[2] + xyzw.m_xyzw[2], m_xyzw[3] + xyzw.m_xyzw[3]);
}

XYZW LNLib::XYZW::operator-(const XYZW& xyzw) const
{
	return XYZW(m_xyzw[0] - xyzw.m_xyzw[0], m_xyzw[1] - xyzw.m_xyzw[1], m_xyzw[2] - xyzw.m_xyzw[2], m_xyzw[3] - xyzw.m_xyzw[3]);
}

XYZW& LNLib::XYZW::operator+=(const XYZW& xyzw)
{
	m_xyzw[0] += xyzw.m_xyzw[0];
	m_xyzw[1] += xyzw.m_xyzw[1];
	m_xyzw[2] += xyzw.m_xyzw[2];
	m_xyzw[3] += xyzw.m_xyzw[3];
	return *this;
}


XYZW LNLib::operator*(const XYZW& source, const double d)
{
	return XYZW(source.GetWX() * d, source.GetWY() * d, source.GetWZ() * d, source.GetW() *d);
}

XYZW LNLib::operator*(const double& d, const XYZW& source)
{
	return XYZW(source.GetWX() * d, source.GetWY() * d, source.GetWZ() * d, source.GetW() * d);
}

XYZW LNLib::operator/(const XYZW& source, double d)
{
	return XYZW(source.GetWX() / d, source.GetWY() / d, source.GetWZ() / d, source.GetW()/d);
}
