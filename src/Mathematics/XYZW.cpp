#include "XYZW.h"
#include "XYZ.h"
#include "MathUtils.h"

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
	m_xyzw[0] = xyz.GetX();
	m_xyzw[1] = xyz.GetY();
	m_xyzw[2] = xyz.GetZ();
	m_xyzw[3] = w;
}

LNLib::XYZW::XYZW(double wx, double wy, double wz, double w)
{
	m_xyzw[0] = wx;
	m_xyzw[1] = wy;
	m_xyzw[2] = wz;
	m_xyzw[3] = w;
}

void XYZW::SetWX(const double wx)
{
	m_xyzw[0] = wx;
}

double XYZW::GetWX() const { return m_xyzw[0]; };

void XYZW::SetWY(const double wy)
{
	m_xyzw[1] = wy;
}


double XYZW::GetWY() const { return m_xyzw[1]; };

void XYZW::SetWZ(const double wz)
{
	m_xyzw[2] = wz;
}

double XYZW::GetWZ() const { return m_xyzw[2]; }

void XYZW::SetW(const double w)
{
	m_xyzw[3] = w;
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

XYZ LNLib::XYZW::ToXYZ(bool divideWeight)
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
