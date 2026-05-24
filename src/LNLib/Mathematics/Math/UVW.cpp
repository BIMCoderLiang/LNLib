/*
 * Author:
 * 2026/05/24 - Yuqing Liang (BIMCoder Liang)
 * bim.frankliang@foxmail.com
 * 
 *
 * Use of this source code is governed by a LGPL-2.1 license that can be found in
 * the LICENSE file.
 */

#include "UVW.h"

#include <cmath>

using namespace LNLib;

LNLib::UVW::UVW()
{
	m_uvw[0] = 0;
	m_uvw[1] = 0;
	m_uvw[2] = 0;
}


UVW::UVW(double u, double v, double w) {

	m_uvw[0] = u;
	m_uvw[1] = v;
	m_uvw[2] = w;
}


void UVW::SetU(const double x)
{
	m_uvw[0] = x;
}

double UVW::GetU() const { return m_uvw[0]; };

void UVW::SetV(const double y)
{
	m_uvw[1] = y;
}

double UVW::GetV() const { return m_uvw[1]; };

void UVW::SetW(const double z)
{
	m_uvw[2] = z;
}

double UVW::GetW() const { return m_uvw[2]; };


double UVW::U() const
{
	return m_uvw[0];
}

double& UVW::U()
{
	return m_uvw[0];
}

double UVW::V() const
{
	return m_uvw[1];
}

double& UVW::V()
{
	return m_uvw[1];
}
double UVW::W() const
{
	return m_uvw[2];
}

double& UVW::W()
{
	return m_uvw[2];
}