#pragma once

#include "LNLibDefinitions.h"

namespace LNLib
{

	class LNLIB_EXPORT ValidationUtils
	{
	public:

		static bool IsInRange(double input, double min, double max);

		static bool IsValidBezier(unsigned int degree, unsigned int controlPointsCount);

		static bool IsValidBspline(unsigned int degree, unsigned int knotVectorCount, unsigned int controlPointsCount);

		static bool IsValidNurbs(unsigned int degree, unsigned int knotVectorCount, unsigned int controlPointsCount, unsigned int weightsCount);
	};
}

