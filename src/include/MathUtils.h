#pragma once

#include "LNLibDefinitions.h"
#include "Constants.h"

namespace LNLib
{
	class LNLIB_EXPORT MathUtils
	{

	public:

		static bool IsAlmostEqualTo(double value1, double value2, double tolerance = Constants::DoubleEpsilon);

		static bool IsGreaterThan(double value1, double value2, double tolerance = Constants::DoubleEpsilon);

		static bool IsGreaterThanOrEqual(double value1, double value2, double tolerance = Constants::DoubleEpsilon);

		static bool IsLessThan(double value1, double value2, double tolerance = Constants::DoubleEpsilon);

		static bool IsLessThanOrEqual(double value1, double value2, double tolerance = Constants::DoubleEpsilon);

		static bool IsInfinite(double value);

		static bool IsNaN(double value);

		static double RadiansToAngle(double radians);

		static double AngleToRadians(double angle);

		static int Factorial(unsigned int number);

		static double Binomial(unsigned int number, unsigned int i);
	};
}


