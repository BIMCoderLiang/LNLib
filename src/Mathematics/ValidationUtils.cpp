#include "ValidationUtils.h"
#include "MathUtils.h"

bool LNLib::ValidationUtils::IsInRange(double input, double min, double max)
{
	return (MathUtils::IsGreaterThanOrEqual(input, min) && MathUtils::IsLessThanOrEqual(input, max));
}

bool LNLib::ValidationUtils::IsValidBezier(unsigned int degree, unsigned int controlPointsCount)
{
	return controlPointsCount == degree + 1;
}

bool LNLib::ValidationUtils::IsValidBspline(unsigned int degree, unsigned int knotVectorCount, unsigned int controlPointsCount)
{
	return knotVectorCount - 1 == degree + controlPointsCount;
}

bool LNLib::ValidationUtils::IsValidNurbs(unsigned int degree, unsigned int knotVectorCount, unsigned int controlPointsCount, unsigned int weightsCount)
{
	return knotVectorCount - 1 == degree + controlPointsCount && 
			controlPointsCount == weightsCount;
}
