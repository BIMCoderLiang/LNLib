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
	class LNLIB_EXPORT UV
	{

	public:

		UV(double u, double v);

	public:

		void SetU(const double x);
		double GetU() const;
		void SetV(const double y);
		double GetV() const;

		double U() const;
		double& U();
		double V() const;
		double& V();

	public:

		bool IsZero(const double epsilon = Constants::DoubleEpsilon) const;
		bool IsUnit(const double epsilon = Constants::DoubleEpsilon) const;
		bool IsAlmostEqualTo(const UV& another) const;
		double Length() const;
		double SqrLength() const;
		double Normalize();
		UV Add(const UV& another) const;
		UV Substract(const UV& another) const;
		UV Negative() const;
		double DotProduct(const UV& another) const;
		double CrossProduct(const UV& another) const;


	public:

		UV& operator =(const UV& xyz);
		double& operator[](int index);
		const double& operator[](int index) const;
		UV operator +(const UV& xyz) const;
		UV operator -(const UV& xyz) const;
		double operator *(const UV& xyz) const;
		UV& operator *=(const double& d);
		UV& operator /=(const double& d);
		UV& operator +=(const UV& xyz);
		UV& operator -=(const UV& xyz);
		UV  operator-() const;

	private:

		double m_uv[2];
	};

	UV operator *(const UV& source, const double d);
	UV operator *(const double& d, const UV& source);
	double operator ^(const UV& xyz1, const UV& xyz2);
	UV operator /(const UV& source, double d);
}



