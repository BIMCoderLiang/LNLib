/*
 * Author:
 * 2026/04/19 - Yuqing Liang (BIMCoder Liang)
 * bim.frankliang@foxmail.com
 * 
 *
 * Use of this source code is governed by a LGPL-2.1 license that can be found in
 * the LICENSE file.
 */

#pragma once

#include "Constants.h"
#include "LNLibDefinitions.h"

namespace LNLib
{
	/// <summary>
	/// Represents three-dimension location/vector/offset
	/// </summary>
	class LNLIB_EXPORT UVW
	{

	public:

		UVW();
		UVW(double u, double v ,double w);

	public:

		void SetU(const double x);
		double GetU() const;
		void SetV(const double y);
		double GetV() const;
		void SetW(const double z);
		double GetW() const;

		double U() const;
		double& U();
		double V() const;
		double& V();
		double W() const;
		double& W();

	private:

		double m_uvw[3];
	};
}



