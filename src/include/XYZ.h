/*
 * Author:
 * 2023/06/08 - Yuqing Liang (BIMCoder Liang)
 * bim.frankliang@foxmail.com
 * 微信公众号：BIMCoder梁老师
 *
 * Use of this source code is governed by a GPL-3.0 license license that can be found in
 * the LICENSE file.
 */

#pragma once
#include "Constants.h"
#include "LNLibDefinitions.h"

namespace LNLib
{
	class LNLIB_EXPORT XYZ
	{

	public:
		XYZ();
		XYZ(double x, double y, double z);

	public:

		void SetX(const double x);
		double GetX() const;
		void SetY(const double y);
		double GetY() const;
		void SetZ(const double z);
		double GetZ() const;

		double X() const;
		double& X();
		double Y() const;
		double& Y();
		double Z() const;
		double& Z();

	public:

		bool IsZero(const double epsilon = Constants::DoubleEpsilon) const;
		bool IsUnit(const double epsilon = Constants::DoubleEpsilon) const;
		bool IsAlmostEqualTo(const XYZ& another) const;
		double Length() const;
		double SqrLength() const;
		XYZ Normalize();
		XYZ Add(const XYZ& another) const;
		XYZ Substract(const XYZ& another) const;
		XYZ Negative() const;
		double DotProduct(const XYZ& another) const;
		XYZ CrossProduct(const XYZ& another) const;

	public:
		double Distance(const XYZ& another) const;

	public:

		XYZ& operator =(const XYZ& xyz);
		double& operator[](int index);
		const double& operator[](int index) const;
		XYZ operator +(const XYZ& xyz) const;
		XYZ operator -(const XYZ& xyz) const;
		double operator *(const XYZ& xyz) const;
		XYZ& operator *=(const double& d);
		XYZ& operator /=(const double& d);
		XYZ& operator +=(const XYZ& xyz);
		XYZ& operator -=(const XYZ& xyz);
		XYZ  operator-() const;

	private:

		double m_xyz[3];

	};


	XYZ operator *(const XYZ& source, const double d);
	XYZ operator *(const double& d, const XYZ& source);
	XYZ operator ^(const XYZ& xyz1, const XYZ& xyz2);
	XYZ operator /(const XYZ& source, double d);
}


