/**
 * Author:
 * 2023/06/08 - Yuqing Liang (BIMCoder Liang)
 * bim.frankliang@foxmail.com
 * 微信公众号：BIMCoder梁老师
 *
 * Use of this source code is governed by a GPL-3.0 license license that can be found in
 * the LICENSE file.
 */

#include "MathUtils.h"
#include <limits>


bool LNLib::MathUtils::IsAlmostEqualTo(double value1, double value2, double tolerance)
{
    if (IsNaN(value1) || IsNaN(value2))
        return false;
    double eps = (abs(value1) + abs(value2) + 10) * tolerance;
    double delta = value1 - value2;
    return (-eps < delta) && (eps > delta);
}

bool LNLib::MathUtils::IsGreaterThan(double value1, double value2, double tolerance)
{
    if (IsNaN(value1) || IsNaN(value2))
        return false;
    return value1 > value2 && !IsAlmostEqualTo(value1, value2, tolerance);
}

bool LNLib::MathUtils::IsGreaterThanOrEqual(double value1, double value2, double tolerance)
{
    if (IsNaN(value1) || IsNaN(value2))
        return false;
    return (value1 - value2 > tolerance) || IsAlmostEqualTo(value1, value2, tolerance);
}

bool LNLib::MathUtils::IsLessThan(double value1, double value2, double tolerance)
{
    if (IsNaN(value1) || IsNaN(value2))
        return false;
    return value1 < value2 && !IsAlmostEqualTo(value1, value2, tolerance);
}

bool LNLib::MathUtils::IsLessThanOrEqual(double value1, double value2, double tolerance)
{
    if (IsNaN(value1) || IsNaN(value2))
        return false;
    return value1 < value2 || IsAlmostEqualTo(value1, value2, tolerance);
}

bool LNLib::MathUtils::IsInfinite(double value)
{
    constexpr double maxValue = (std::numeric_limits<double>::max)();
    double minValue = -maxValue;
    return !(minValue <= value && value <= maxValue);
}

bool LNLib::MathUtils::IsNaN(double value)
{
    return value != value;
}

double LNLib::MathUtils::RadiansToAngle(double radians)
{
    return radians * 180.0 / Constants::Pi;
}

double LNLib::MathUtils::AngleToRadians(double angle)
{
    return angle * Constants::Pi / 180.0;
}

int LNLib::MathUtils::Factorial(unsigned int number)
{
    if (number == 0)
        return 1;
    else
        return number * Factorial(number - 1);
}

double LNLib::MathUtils::Binomial(unsigned int number, unsigned int i)
{
    return Factorial(number) / (Factorial(i) * Factorial(number - 1));
}
