#include "ValidationUtils.h"
#include "XYZW.h"
#include "Constants.h"
#include "MathUtils.h"
#include <algorithm>


namespace LNLib
{
	static double GetCoefficient(int index, unsigned int degree)
	{
		return index / degree;
	}
}

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

bool LNLib::ValidationUtils::IsValidDegreeReduction(unsigned int degree)
{
	return degree > 1;
}

bool LNLib::ValidationUtils::IsValidKnotVector(const std::vector<double>& knotVector)
{
	return std::is_sorted(knotVector.begin(), knotVector.end());
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

	return Constants::DistanceEpsilon * minWeight / (1 + std::abs(maxDistance));
}

double LNLib::ValidationUtils::ComputeMaxErrorOfBezierReduction(unsigned int degree, const std::vector<XYZW>& current, const std::vector<XYZW>& reduction)
{
	int r = (degree - 1)/2;

	if (degree % 2 == 0)
	{
		XYZW cr = current[r];
		XYZW rmo = reduction[r-1];
		XYZ Pr = (cr.ToXYZ(true) - GetCoefficient(r, degree) * rmo.ToXYZ(true)) / (1 - GetCoefficient(r, degree));

		XYZW cpt = current[r + 2];
		XYZW rpt = reduction[r + 2];
		XYZ Pr1 = (cpt.ToXYZ(true) - (1 - GetCoefficient(r + 2, degree)) * rpt.ToXYZ(true)) / GetCoefficient(r + 2,degree);

		XYZW cpo = current[r + 1];
		XYZ temp = 0.5 * (Pr + Pr1);
		return std::abs(cpo.ToXYZ(true).Distance(temp));
	}
	else
	{
		XYZW cr = current[r];
		XYZW rmo = reduction[r - 1];
		XYZ PLr = (cr.ToXYZ(true) - GetCoefficient(r, degree) * rmo.ToXYZ(true)) / (1 - GetCoefficient(r, degree));

		XYZW cpo = current[r + 1];
		XYZW rpo = reduction[r + 1];
		XYZ PRr = (cpo.ToXYZ(true) - (1 - GetCoefficient(r + 1, degree)) * rpo.ToXYZ(true)) / GetCoefficient(r + 1, degree);

		return std::abs(PLr.Distance(PRr));
	}
}
