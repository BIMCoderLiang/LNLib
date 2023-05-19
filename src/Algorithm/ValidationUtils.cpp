#include "ValidationUtils.h"
#include "XYZW.h"
#include "Constants.h"
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

double LNLib::ValidationUtils::ComputeCurveModifyTolerance(const std::vector<XYZW>& controlPoints)
{
	double minWeight = 1.0;
	double maxDistance = 0.0;

	int size = static_cast<int>(controlPoints.size());
	for (int i = 0; i < size; i++)
	{
		XYZW temp = controlPoints[i];
		minWeight = std::min(minWeight, temp.GetW());
		maxDistance = std::max(maxDistance, temp.ToXYZ(true).Length());
	}

	return Constants::DoubleEpsilon * minWeight / (1 + std::abs(maxDistance));
}



double LNLib::ValidationUtils::ComputeMaxErrorOfBezierReduction(unsigned int degree, const std::vector<XYZW>& current, const std::vector<XYZW>& reduction)
{
	int r = (degree - 1)/2;

	if (degree % 2 == 0)
	{
		XYZW Pr = (current[r + 1] - (1 - (r + 1) / degree) * reduction[r + 1]) / ((r + 1) / degree);
		XYZW Pr1 = (current[r + 1 + 1] - (1 - (r + 1 + 1) / degree) * reduction[r + 1 + 1]) / ((r + 1 + 1) / degree);

		XYZW temp = 0.5 * (Pr + Pr1);
		return abs(current[r + 1].Distance(temp));
	}
	else
	{
		XYZW Plr = (current[r] - (r / degree) * reduction[r + 1]) / (1- (r / degree));
		XYZW PRr = (current[r + 1] - (1 - ((r + 1) / degree)) * reduction[r + 1]) / ((r + 1) / degree);

		return abs(Plr.Distance(PRr));
	}
}
