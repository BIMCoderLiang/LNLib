#pragma once
#include "LNLibDefinitions.h"
#include "XYZ.h"

namespace LNLib
{

	class LNLIB_EXPORT XYZW
	{
	public:

		XYZW();
		XYZW(XYZ xyz, double w);
		XYZW(double wx, double wy, double wz, double w);

	public:

		void SetWX(const double wx);
		double GetWX() const;
		void SetWY(const double wy);
		double GetWY() const;
		void SetWZ(const double wz);
		double GetWZ() const;
		void SetW(const double w);
		double GetW() const;

		double WX() const;
		double& WX();
		double WY() const;
		double& WY();
		double WZ() const;
		double& WZ();
		double W() const;
		double& W();

	public:
		XYZ ToXYZ(bool divideWeight);

	public:
		XYZW& operator +=(const XYZW& xyzw);

	private:

		double m_xyzw[4];
	};

	XYZW operator *(const XYZW& source, const double d);
	XYZW operator *(const double& d, const XYZW& source);
}

