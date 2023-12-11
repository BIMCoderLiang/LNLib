/*
 * Author:
 * 2023/06/08 - Yuqing Liang (BIMCoder Liang)
 * bim.frankliang@foxmail.com
 * 微信公众号：BIMCoder梁老师
 *
 * Use of this source code is governed by a GPL-3.0 license that can be found in
 * the LICENSE file.
 */

#include "NurbsSurface.h"
#include "Polynomials.h"
#include "UV.h"
#include "XYZ.h"
#include "XYZW.h"
#include "Matrix4d.h"
#include "MathUtils.h"
#include "NurbsCurve.h"
#include "BsplineSurface.h"
#include "Projection.h"
#include "Intersection.h"
#include "Interpolation.h"
#include "ValidationUtils.h"
#include "KnotVectorUtils.h"
#include "ControlPointsUtils.h"
#include "LNLibExceptions.h"
#include <algorithm>

namespace LNLib
{
	std::vector<int> GetIndex(int size)
	{
		std::vector<int> ind(2 * (size - 1) + 2);
		ind[0] = 0;
		ind[ind.size() - 1] = 3 * size - 3;
		int ii = 1;
		int jj = 1;
		for (int i = 0; i < size - 1; i++)
		{
			ind[ii] = jj;
			ind[ii + 1] = jj + 1;
			ii = ii + 2;
			jj = jj + 3;
		}
		return ind;
	}
}

LNLib::XYZ LNLib::NurbsSurface::GetPointOnSurface(const LN_Surface& surface, UV uv)
{
	int degreeU = surface.DegreeU;
	int degreeV = surface.DegreeV;
	std::vector<double> knotVectorU = surface.KnotVectorU;
	std::vector<double> knotVectorV = surface.KnotVectorV;
	std::vector<std::vector<XYZW>> controlPoints = surface.ControlPoints;

	VALIDATE_ARGUMENT(degreeU > 0, "degreeU", "Degree must greater than zero.");
	VALIDATE_ARGUMENT(degreeV > 0, "degreeV", "Degree must greater than zero.");
	VALIDATE_ARGUMENT(knotVectorU.size() > 0, "knotVectorU", "KnotVector size must greater than zero.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidKnotVector(knotVectorU), "knotVectorU", "KnotVector must be a nondecreasing sequence of real numbers.");
	VALIDATE_ARGUMENT_RANGE(uv.GetU(), knotVectorU[0], knotVectorU[knotVectorU.size() - 1]);
	VALIDATE_ARGUMENT(knotVectorV.size() > 0, "knotVectorV", "KnotVector size must greater than zero.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidKnotVector(knotVectorV), "knotVectorV", "KnotVector must be a nondecreasing sequence of real numbers.");
	VALIDATE_ARGUMENT_RANGE(uv.GetV(), knotVectorV[0], knotVectorV[knotVectorV.size() - 1]);
	VALIDATE_ARGUMENT(controlPoints.size() > 0, "controlPoints", "ControlPoints must contains one point at least.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidNurbs(degreeU, knotVectorU.size(), controlPoints.size()), "controlPoints", "Arguments must fit: m = n + p + 1");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidNurbs(degreeV, knotVectorV.size(), controlPoints[0].size()), "controlPoints", "Arguments must fit: m = n + p + 1");

	XYZW result = BsplineSurface::GetPointOnSurface(degreeU, degreeV, knotVectorU, knotVectorV, uv, controlPoints);
	return result.ToXYZ(true);
}


std::vector<std::vector<LNLib::XYZ>> LNLib::NurbsSurface::ComputeRationalSurfaceDerivatives(const LN_Surface& surface, int derivative, UV uv)
{
	int degreeU = surface.DegreeU;
	int degreeV = surface.DegreeV;
	std::vector<double> knotVectorU = surface.KnotVectorU;
	std::vector<double> knotVectorV = surface.KnotVectorV;
	std::vector<std::vector<XYZW>> controlPoints = surface.ControlPoints;

	VALIDATE_ARGUMENT(degreeU > 0, "degreeU", "Degree must greater than zero.");
	VALIDATE_ARGUMENT(degreeV > 0, "degreeU", "Degree must greater than zero.");
	VALIDATE_ARGUMENT(derivative > 0, "derivative", "derivative must greater than zero.");
	VALIDATE_ARGUMENT(knotVectorU.size() > 0, "knotVectorU", "KnotVector size must greater than zero.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidKnotVector(knotVectorU), "knotVectorU", "KnotVector must be a nondecreasing sequence of real numbers.");
	VALIDATE_ARGUMENT_RANGE(uv.GetU(), knotVectorU[0], knotVectorU[knotVectorU.size() - 1]);
	VALIDATE_ARGUMENT(knotVectorV.size() > 0, "knotVectorV", "KnotVector size must greater than zero.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidKnotVector(knotVectorV), "knotVectorV", "KnotVector must be a nondecreasing sequence of real numbers.");
	VALIDATE_ARGUMENT_RANGE(uv.GetV(), knotVectorV[0], knotVectorV[knotVectorV.size() - 1]);
	VALIDATE_ARGUMENT(controlPoints.size() > 0, "controlPoints", "ControlPoints must contains one point at least.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidNurbs(degreeU, knotVectorU.size(), controlPoints.size()), "controlPoints", "Arguments must fit: m = n + p + 1");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidNurbs(degreeV, knotVectorV.size(), controlPoints[0].size()), "controlPoints", "Arguments must fit: m = n + p + 1");

	std::vector<std::vector<XYZ>> derivatives(derivative + 1, std::vector<XYZ>(derivative + 1));

	std::vector<std::vector<XYZW>> ders = BsplineSurface::ComputeDerivatives(degreeU, degreeV, derivative, knotVectorU, knotVectorV, uv, controlPoints);
	std::vector<std::vector<XYZ>> Aders(derivative + 1, std::vector<XYZ>(derivative + 1));
	std::vector<std::vector<double>> wders(derivative + 1, std::vector<double>(derivative + 1));
	for (int i = 0; i < ders.size(); i++)
	{
		for (int j = 0; j < ders[0].size(); j++)
		{
			Aders[i][j] = ders[i][j].ToXYZ(false);
			wders[i][j] = ders[i][j].GetW();
		}
	}

	for (int k = 0; k <= derivative; k++)
	{
		for (int l = 0; l <= derivative - k; l++)
		{
			XYZ v = Aders[k][l];
			for ( int j = 1; j <= l; j++)
			{
				v = v - MathUtils::Binomial(l, j) * wders[0][j] * derivatives[k][l - j];
			}

			for (int i = 1; i <= k; i++)
			{
				v = v - MathUtils::Binomial(k, i) * wders[i][0] * derivatives[k - i][l];

				XYZ v2 = XYZ(0,0,0);
				for (int j = 1; j <= l; j++)
				{
					v2 = v2 + MathUtils::Binomial(l, j) * wders[i][j] * derivatives[k - i][l - j];
				}
				v = v - MathUtils::Binomial(k, i) * v2;
			}
			derivatives[k][l] = v / wders[0][0];
		}
	}
	return derivatives;
}

double LNLib::NurbsSurface::Curvature(const LN_Surface& surface, SurfaceCurvature curvature, UV uv)
{
	int degreeU = surface.DegreeU;
	int degreeV = surface.DegreeV;
	std::vector<double> knotVectorU = surface.KnotVectorU;
	std::vector<double> knotVectorV = surface.KnotVectorV;
	std::vector<std::vector<XYZW>> controlPoints = surface.ControlPoints;

	VALIDATE_ARGUMENT(degreeU > 0, "degreeU", "Degree must greater than zero.");
	VALIDATE_ARGUMENT(degreeV > 0, "degreeU", "Degree must greater than zero.");
	VALIDATE_ARGUMENT(knotVectorU.size() > 0, "knotVectorU", "KnotVector size must greater than zero.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidKnotVector(knotVectorU), "knotVectorU", "KnotVector must be a nondecreasing sequence of real numbers.");
	VALIDATE_ARGUMENT_RANGE(uv.GetU(), knotVectorU[0], knotVectorU[knotVectorU.size() - 1]);
	VALIDATE_ARGUMENT(knotVectorV.size() > 0, "knotVectorV", "KnotVector size must greater than zero.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidKnotVector(knotVectorV), "knotVectorV", "KnotVector must be a nondecreasing sequence of real numbers.");
	VALIDATE_ARGUMENT_RANGE(uv.GetV(), knotVectorV[0], knotVectorV[knotVectorV.size() - 1]);
	VALIDATE_ARGUMENT(controlPoints.size() > 0, "controlPoints", "ControlPoints must contains one point at least.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidNurbs(degreeU, knotVectorU.size(), controlPoints.size()), "controlPoints", "Arguments must fit: m = n + p + 1");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidNurbs(degreeV, knotVectorV.size(), controlPoints[0].size()), "controlPoints", "Arguments must fit: m = n + p + 1");

	std::vector<std::vector<XYZ>> ders = ComputeRationalSurfaceDerivatives(surface, 2, uv);
	XYZ Suu = ders[2][0];
	XYZ Svv = ders[0][2];
	XYZ Suv = ders[1][1];

	XYZ Su = ders[1][0];
	XYZ Sv = ders[0][1];
	XYZ normal = Normal(surface, uv);
	
	double L = Suu.DotProduct(normal);
	double M = Suv.DotProduct(normal);
	double N = Svv.DotProduct(normal);

	double E = Su.DotProduct(Su);
	double F = Su.DotProduct(Sv);
	double G = Sv.DotProduct(Sv);

	double denominator = E * G - F * F;
	if (MathUtils::IsAlmostEqualTo(denominator, 0.0))
	{
		return 0.0;
	}

	double K = (L * N - M * M) / denominator;
	double H = (E * N + G * L - 2 * F * M) / (2 * denominator);
	double k1 = H + sqrt(abs(H * H - K));
	double k2 = H - sqrt(abs(H * H - K));

	if (curvature == SurfaceCurvature::Gauss)
	{
		return K;
	}
	else if (curvature == SurfaceCurvature::Mean)
	{
		return H;
	}
	else if (curvature == SurfaceCurvature::Maximum)
	{
		return k1;
	}
	else if (curvature == SurfaceCurvature::Minimum)
	{
		return k2;
	}
	else if (curvature == SurfaceCurvature::Abs)
	{
		return abs(k1) + abs(k2);
	}
	else if (curvature == SurfaceCurvature::Rms)
	{
		return sqrt(k1* k1 + k2 * k2);
	}
	return 0.0;
}

LNLib::XYZ LNLib::NurbsSurface::Normal(const LN_Surface& surface, UV uv)
{
	int degreeU = surface.DegreeU;
	int degreeV = surface.DegreeV;
	std::vector<double> knotVectorU = surface.KnotVectorU;
	std::vector<double> knotVectorV = surface.KnotVectorV;
	std::vector<std::vector<XYZW>> controlPoints = surface.ControlPoints;

	VALIDATE_ARGUMENT(degreeU > 0, "degreeU", "Degree must greater than zero.");
	VALIDATE_ARGUMENT(degreeV > 0, "degreeU", "Degree must greater than zero.");
	VALIDATE_ARGUMENT(knotVectorU.size() > 0, "knotVectorU", "KnotVector size must greater than zero.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidKnotVector(knotVectorU), "knotVectorU", "KnotVector must be a nondecreasing sequence of real numbers.");
	VALIDATE_ARGUMENT_RANGE(uv.GetU(), knotVectorU[0], knotVectorU[knotVectorU.size() - 1]);
	VALIDATE_ARGUMENT(knotVectorV.size() > 0, "knotVectorV", "KnotVector size must greater than zero.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidKnotVector(knotVectorV), "knotVectorV", "KnotVector must be a nondecreasing sequence of real numbers.");
	VALIDATE_ARGUMENT_RANGE(uv.GetV(), knotVectorV[0], knotVectorV[knotVectorV.size() - 1]);
	VALIDATE_ARGUMENT(controlPoints.size() > 0, "controlPoints", "ControlPoints must contains one point at least.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidNurbs(degreeU, knotVectorU.size(), controlPoints.size()), "controlPoints", "Arguments must fit: m = n + p + 1");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidNurbs(degreeV, knotVectorV.size(), controlPoints[0].size()), "controlPoints", "Arguments must fit: m = n + p + 1");

	std::vector<std::vector<XYZ>> derivatives = ComputeRationalSurfaceDerivatives(surface, 1, uv);
	return derivatives[1][0].Normalize().CrossProduct(derivatives[0][1]).Normalize();
}

void LNLib::NurbsSurface::Swap(const LN_Surface& surface, LN_Surface& result)
{
	int degreeU = surface.DegreeU;
	int degreeV = surface.DegreeV;
	std::vector<double> knotVectorU = surface.KnotVectorU;
	std::vector<double> knotVectorV = surface.KnotVectorV;
	std::vector<std::vector<XYZW>> controlPoints = surface.ControlPoints;

	VALIDATE_ARGUMENT(degreeU > 0, "degreeU", "Degree must greater than zero.");
	VALIDATE_ARGUMENT(degreeV > 0, "degreeU", "Degree must greater than zero.");
	VALIDATE_ARGUMENT(knotVectorU.size() > 0, "knotVectorU", "KnotVector size must greater than zero.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidKnotVector(knotVectorU), "knotVectorU", "KnotVector must be a nondecreasing sequence of real numbers.");
	VALIDATE_ARGUMENT(knotVectorV.size() > 0, "knotVectorV", "KnotVector size must greater than zero.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidKnotVector(knotVectorV), "knotVectorV", "KnotVector must be a nondecreasing sequence of real numbers.");
	VALIDATE_ARGUMENT(controlPoints.size() > 0, "controlPoints", "ControlPoints must contains one point at least.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidNurbs(degreeU, knotVectorU.size(), controlPoints.size()), "controlPoints", "Arguments must fit: m = n + p + 1");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidNurbs(degreeV, knotVectorV.size(), controlPoints[0].size()), "controlPoints", "Arguments must fit: m = n + p + 1");

	std::vector<std::vector<XYZW>> transposedControlPoints;
	MathUtils::Transpose(surface.ControlPoints, transposedControlPoints);
	result.DegreeU = surface.DegreeV;
	result.DegreeV = surface.DegreeU;
	result.KnotVectorU = surface.KnotVectorV;
	result.KnotVectorV = surface.KnotVectorU;
	result.ControlPoints = transposedControlPoints;
}

void LNLib::NurbsSurface::Reverse(const LN_Surface& surface, SurfaceDirection direction, LN_Surface& result)
{
	int degreeU = surface.DegreeU;
	int degreeV = surface.DegreeV;
	std::vector<double> knotVectorU = surface.KnotVectorU;
	std::vector<double> knotVectorV = surface.KnotVectorV;
	std::vector<std::vector<XYZW>> controlPoints = surface.ControlPoints;

	VALIDATE_ARGUMENT(degreeU > 0, "degreeU", "Degree must greater than zero.");
	VALIDATE_ARGUMENT(degreeV > 0, "degreeU", "Degree must greater than zero.");
	VALIDATE_ARGUMENT(knotVectorU.size() > 0, "knotVectorU", "KnotVector size must greater than zero.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidKnotVector(knotVectorU), "knotVectorU", "KnotVector must be a nondecreasing sequence of real numbers.");
	VALIDATE_ARGUMENT(knotVectorV.size() > 0, "knotVectorV", "KnotVector size must greater than zero.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidKnotVector(knotVectorV), "knotVectorV", "KnotVector must be a nondecreasing sequence of real numbers.");
	VALIDATE_ARGUMENT(controlPoints.size() > 0, "controlPoints", "ControlPoints must contains one point at least.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidNurbs(degreeU, knotVectorU.size(), controlPoints.size()), "controlPoints", "Arguments must fit: m = n + p + 1");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidNurbs(degreeV, knotVectorV.size(), controlPoints[0].size()), "controlPoints", "Arguments must fit: m = n + p + 1");

	result.DegreeU = degreeU;
	result.DegreeV = degreeV;

	if (direction == SurfaceDirection::All || direction == SurfaceDirection::UDirection)
	{
		int size = knotVectorU.size();
		std::vector<double> reversedKnotVectorU(size);
		double min = knotVectorU[0];
		double max = knotVectorU[size - 1];
		reversedKnotVectorU[0] = min;
		for (int i = 1; i < size; i++)
		{
			reversedKnotVectorU[i] = reversedKnotVectorU[i - 1] + (knotVectorU[size - i] - knotVectorU[size - i - 1]);
		}
		result.KnotVectorU = reversedKnotVectorU;
	}
	
	if (direction == SurfaceDirection::All || direction == SurfaceDirection::VDirection)
	{
		int size = knotVectorV.size();
		std::vector<double> reversedKnotVectorV(size);
		double min = knotVectorV[0];
		double max = knotVectorV[size - 1];
		reversedKnotVectorV[0] = min;
		for (int i = 1; i < size; i++)
		{
			reversedKnotVectorV[i] = reversedKnotVectorV[i - 1] + (knotVectorV[size - i] - knotVectorV[size - i - 1]);
		}
		result.KnotVectorV = reversedKnotVectorV;
	}

	if (direction == SurfaceDirection::UDirection)
	{
		int row = controlPoints.size();
		int column = controlPoints[0].size();

		for (int i = 0; i < column; i++) {
			int start = 0;
			int end = row - 1;
			while (start < end) {
				std::swap(controlPoints[start][i], controlPoints[end][i]);
				start++;
				end--;
			}
		}
		result.ControlPoints = controlPoints;
	}
	else if (direction == SurfaceDirection::VDirection)
	{
		int size = controlPoints.size();
		std::vector<std::vector<XYZW>> newControlPoints(size);
		for (int i = 0; i < size; i++)
		{
			std::reverse(controlPoints[i].begin(), controlPoints[i].end());
			newControlPoints.emplace_back(controlPoints[i]);
		}
		result.ControlPoints = newControlPoints;
	}
	else
	{
		std::reverse(controlPoints.begin(), controlPoints.end());
		result.ControlPoints = controlPoints;
	}
	
}

void LNLib::NurbsSurface::InsertKnot(const LN_Surface& surface, double insertKnot, int times, bool isUDirection, LN_Surface& result)
{
	int degreeU = surface.DegreeU;
	int degreeV = surface.DegreeV;
	std::vector<double> knotVectorU = surface.KnotVectorU;
	std::vector<double> knotVectorV = surface.KnotVectorV;
	std::vector<std::vector<XYZW>> controlPoints = surface.ControlPoints;

	VALIDATE_ARGUMENT(degreeU > 0, "degreeU", "Degree must greater than zero.");
	VALIDATE_ARGUMENT(degreeV > 0, "degreeV", "Degree must greater than zero.");
	VALIDATE_ARGUMENT(knotVectorU.size() > 0, "knotVectorU", "KnotVector size must greater than zero.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidKnotVector(knotVectorU), "knotVectorU", "KnotVector must be a nondecreasing sequence of real numbers.");
	VALIDATE_ARGUMENT(knotVectorV.size() > 0, "knotVectorV", "KnotVector size must greater than zero.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidKnotVector(knotVectorV), "knotVectorV", "KnotVector must be a nondecreasing sequence of real numbers.");
	VALIDATE_ARGUMENT(controlPoints.size() > 0, "controlPoints", "ControlPoints must contains one point at least.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidNurbs(degreeU, knotVectorU.size(), controlPoints.size()), "controlPoints", "Arguments must fit: m = n + p + 1");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidNurbs(degreeV, knotVectorV.size(), controlPoints[0].size()), "controlPoints", "Arguments must fit: m = n + p + 1");

	if (isUDirection)
	{
		VALIDATE_ARGUMENT(ValidationUtils::IsValidNurbs(degreeU, knotVectorU.size(), controlPoints.size()), "controlPoints", "Arguments must fit: m = n + p + 1");
	}
	else
	{
		VALIDATE_ARGUMENT(ValidationUtils::IsValidNurbs(degreeV, knotVectorV.size(), controlPoints[0].size()), "controlPoints", "Arguments must fit: m = n + p + 1");
	}
	
	int degree = isUDirection ? surface.DegreeU : surface.DegreeV;
	std::vector<double> knotVector = isUDirection ? surface.KnotVectorU : surface.KnotVectorV;
	int knotSpanIndex = Polynomials::GetKnotSpanIndex(degree, knotVector, insertKnot);
	int multiplicity = Polynomials::GetKnotMultiplicity(knotVector, insertKnot);

	if (multiplicity == degree)
	{
		result = surface;
		return;
	}

	if ((times + multiplicity) > degree)
	{
		times = degree - multiplicity;
	}

	std::vector<double> insertedKnotVector(knotVector.size() + times);
	for (int i = 0; i <= knotSpanIndex; i++)
	{
		insertedKnotVector[i] = knotVector[i];
	}

	for (int i = 1; i <= times; i++)
	{
		insertedKnotVector[knotSpanIndex + i] = insertKnot;
	}

	for (int i = knotSpanIndex + 1; i < knotVector.size(); i++)
	{
		insertedKnotVector[i + times] = knotVector[i];
	}

	std::vector<std::vector<double>> alpha(degree - multiplicity,std::vector<double>(times + 1));
	for (int j = 1; j <= times; j++)
	{
		int L = knotSpanIndex - degree + j;
		for (int i = 0; i <= degree - j - multiplicity; i++)
		{
			alpha[i][j] = (insertKnot - knotVector[L + i]) / (knotVector[i + knotSpanIndex + 1] - knotVector[L + i]);
		}
	}

	std::vector<XYZW> temp(degree + 1);

	int rows = controlPoints.size();
	int columns = controlPoints[0].size();

	result = surface;

	std::vector<std::vector<XYZW>> updatedControlPoints;
	if (isUDirection)
	{
		updatedControlPoints.resize(rows + times,std::vector<XYZW>(columns));
		for (int col = 0; col < columns; col++)
		{
			for (int i = 0; i <= knotSpanIndex - degree; i++)
			{
				updatedControlPoints[i][col] = controlPoints[i][col];
			}

			for (int i = knotSpanIndex - multiplicity; i < rows; i++)
			{
				updatedControlPoints[i + times][col] = controlPoints[i][col];
			}

			for (int i = 0; i < degree - multiplicity + 1; i++)
			{
				temp[i] = controlPoints[knotSpanIndex - degree + i][col];
			}

			int L = 0;
			for (int j = 1; j <= times; j++)
			{
				L = knotSpanIndex - degree + j;
				for (int i = 0; i <= degree - j - multiplicity; i++)
				{
					double a = alpha[i][j];
					temp[i] = a * temp[i + 1] + (1.0 - a) * temp[i];
				}
				updatedControlPoints[L][col] = temp[0];
				updatedControlPoints[knotSpanIndex + times - j - multiplicity][col] = temp[degree - j - multiplicity];
			}

			for (int i = L + 1; i < knotSpanIndex - multiplicity; i++)
			{
				updatedControlPoints[i][col] = temp[i - L];
			}
		}
		result.KnotVectorU = insertedKnotVector;
		result.ControlPoints = updatedControlPoints;
	}
	else
	{
		updatedControlPoints.resize(rows,std::vector<XYZW>(columns+times));
		for (int row = 0; row < rows; row++)
		{
			for (int i = 0; i <= knotSpanIndex - degree; i++)
			{
				updatedControlPoints[row][i] = controlPoints[row][i];
			}

			for (int i = knotSpanIndex - multiplicity; i < columns; i++)
			{
				updatedControlPoints[row][i+times] = controlPoints[row][i];
			}

			for (int i = 0; i < degree - multiplicity + 1; i++)
			{
				temp[i] = controlPoints[row][knotSpanIndex - degree + i];
			}

			int L = 0;
			for (int j = 1; j <= times; j++)
			{
				L = knotSpanIndex - degree + j;
				for (int i = 0; i <= degree - j - multiplicity; i++)
				{
					double a = alpha[i][j];
					temp[i] = a * temp[i + 1] + (1.0 - a) * temp[i];
				}
				updatedControlPoints[row][L] = temp[0];
				updatedControlPoints[row][knotSpanIndex + times - j - multiplicity] = temp[degree - j - multiplicity];
			}

			for (int i = L + 1; i < knotSpanIndex - multiplicity; i++)
			{
				updatedControlPoints[row][i] = temp[i - L];
			}
		}
		result.KnotVectorV = insertedKnotVector;
		result.ControlPoints = updatedControlPoints;
	}
}

void LNLib::NurbsSurface::RefineKnotVector(const LN_Surface& surface, std::vector<double>& insertKnotElements, bool isUDirection, LN_Surface& result)
{
	int degreeU = surface.DegreeU;
	int degreeV = surface.DegreeV;
	std::vector<double> knotVectorU = surface.KnotVectorU;
	std::vector<double> knotVectorV = surface.KnotVectorV;
	std::vector<std::vector<XYZW>> controlPoints = surface.ControlPoints;

	VALIDATE_ARGUMENT(degreeU > 0, "degreeU", "Degree must greater than zero.");
	VALIDATE_ARGUMENT(degreeV > 0, "degreeV", "Degree must greater than zero.");
	VALIDATE_ARGUMENT(knotVectorU.size() > 0, "knotVectorU", "KnotVector size must greater than zero.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidKnotVector(knotVectorU), "knotVectorU", "KnotVector must be a nondecreasing sequence of real numbers.");
	VALIDATE_ARGUMENT(knotVectorV.size() > 0, "knotVectorV", "KnotVector size must greater than zero.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidKnotVector(knotVectorV), "knotVectorV", "KnotVector must be a nondecreasing sequence of real numbers.");
	VALIDATE_ARGUMENT(controlPoints.size() > 0, "controlPoints", "ControlPoints must contains one point at least.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidNurbs(degreeU, knotVectorU.size(), controlPoints.size()), "controlPoints", "Arguments must fit: m = n + p + 1");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidNurbs(degreeV, knotVectorV.size(), controlPoints[0].size()), "controlPoints", "Arguments must fit: m = n + p + 1");
	VALIDATE_ARGUMENT(insertKnotElements.size() > 0, "insertKnotElements", "insertKnotElements size must greater than zero.");

	std::vector<double> tempKnotVector;
	std::vector<std::vector<XYZW>> tempControlPoints;

	result = surface;

	if (isUDirection)
	{
		std::vector<std::vector<XYZW>> transposedControlPoints;
		MathUtils::Transpose(controlPoints, transposedControlPoints);

		for (int i = 0; i < transposedControlPoints.size(); i++)
		{
			LN_Curve tc;
			tc.Degree = degreeU;
			tc.KnotVector = knotVectorU;
			tc.ControlPoints = transposedControlPoints[i];

			LN_Curve newtc;
			NurbsCurve::RefineKnotVector(tc, insertKnotElements, newtc);
			tempControlPoints.emplace_back(newtc.ControlPoints);
			tempKnotVector = newtc.KnotVector;
		}
		result.KnotVectorU = tempKnotVector;
		result.KnotVectorV = knotVectorV;
		std::vector<std::vector<XYZW>> updatedControlPoints;
		MathUtils::Transpose(tempControlPoints, updatedControlPoints);
		result.ControlPoints = updatedControlPoints;
	}
	else
	{
		for (int i = 0; i < controlPoints.size(); i++)
		{
			LN_Curve tc;
			tc.Degree = degreeU;
			tc.KnotVector = knotVectorU;
			tc.ControlPoints = controlPoints[i];

			LN_Curve newtc;
			NurbsCurve::RefineKnotVector(tc, insertKnotElements, newtc);
			tempControlPoints.emplace_back(newtc.ControlPoints);
			tempKnotVector = newtc.KnotVector;
		}
		result.KnotVectorU = knotVectorU;
		result.KnotVectorV = tempKnotVector;
		result.ControlPoints = tempControlPoints;
	}
}

std::vector<std::vector<std::vector<LNLib::XYZW>>> LNLib::NurbsSurface::DecomposeToBeziers(const LN_Surface& surface)
{
	int degreeU = surface.DegreeU;
	int degreeV = surface.DegreeV;
	std::vector<double> knotVectorU = surface.KnotVectorU;
	std::vector<double> knotVectorV = surface.KnotVectorV;
	std::vector<std::vector<XYZW>> controlPoints = surface.ControlPoints;

	VALIDATE_ARGUMENT(degreeU > 0, "degreeU", "Degree must greater than zero.");
	VALIDATE_ARGUMENT(degreeV > 0, "degreeV", "Degree must greater than zero.");
	VALIDATE_ARGUMENT(knotVectorU.size() > 0, "knotVectorU", "KnotVector size must greater than zero.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidKnotVector(knotVectorU), "knotVectorU", "KnotVector must be a nondecreasing sequence of real numbers.");
	VALIDATE_ARGUMENT(knotVectorV.size() > 0, "knotVectorV", "KnotVector size must greater than zero.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidKnotVector(knotVectorV), "knotVectorV", "KnotVector must be a nondecreasing sequence of real numbers.");
	VALIDATE_ARGUMENT(controlPoints.size() > 0, "controlPoints", "ControlPoints must contains one point at least.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidNurbs(degreeU, knotVectorU.size(), controlPoints.size()), "controlPoints", "Arguments must fit: m = n + p + 1");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidNurbs(degreeV, knotVectorV.size(), controlPoints[0].size()), "controlPoints", "Arguments must fit: m = n + p + 1");

	std::vector<std::vector<std::vector<LNLib::XYZW>>> decomposedControlPoints;
	
	int rows = static_cast<int>(controlPoints.size());
	int columns = static_cast<int>(controlPoints[0].size());

	std::vector<std::vector<std::vector<LNLib::XYZW>>> temp(rows - degreeU, std::vector<std::vector<XYZW>>(degreeU + 1, std::vector<XYZW>(columns)));

	int m = rows - 1 + degreeU + 1;
	int a = degreeU;
	int b = degreeU + 1;

	int nb = 0;

	for (int i = 0; i <= degreeU; i++)
	{
		for (int col = 0; col < columns; col++)
		{
			temp[nb][i][col] = controlPoints[i][col];
		}
	}

	std::vector<double> alphaVector(std::max(degreeU,degreeV) + 1);
	while (b < m)
	{
		int i = b;
		while (b < m && MathUtils::IsLessThanOrEqual(knotVectorU[b + 1], knotVectorU[b]))
		{
			b++;
		}
		int multi = b - i + 1;
		if (multi < degreeU)
		{
			double numerator = knotVectorU[b] - knotVectorU[a];
			
			for (int j = degreeU; j > multi; j--)
			{
				alphaVector[j - multi - 1] = numerator / (knotVectorU[a + j]-knotVectorU[a]);
			}
			int r = degreeU - multi;
			for (int j = 1; j <= r; j++)
			{
				int save = r - j;
				int s = multi + j;
				for (int k = degreeU; k >= s; k--)
				{
					double alpha = alphaVector[k - s];
					for (int col = 0; col < columns; col++)
					{
						temp[nb][k][col] = alpha * decomposedControlPoints[nb][k][col] + (1.0 - alpha) * decomposedControlPoints[nb][k - 1][col];
					}
				}
				if (b < m)
				{
					for (int col = 0; col < columns; col++)
					{
						temp[nb + 1][save][col] = decomposedControlPoints[nb][degreeU][col];
					}
				}
			}
		}
		nb++;
		if (b < m)
		{
			for (int i = degreeU - multi; i <= degreeU; i++)
			{
				for (int col = 0; col < columns; col++)
				{
					temp[nb][i][col] = controlPoints[b - degreeU + i][col];
				}
			}
			a = b;
			b++;
		}
	}

	decomposedControlPoints.resize(nb * (columns - degreeV), std::vector<std::vector<XYZW>>(degreeU + 1, std::vector<XYZW>(degreeV + 1)));

	nb = 0;
	for (int np = 0; np < nb; np++)
	{
		for (int i = 0; i <= degreeU; i++)
		{
			for (int j = 0; j <= degreeV; j++)
			{
				decomposedControlPoints[nb][i][j] = temp[np][i][j];
			}
		}

		m = columns + degreeV;
		a = degreeV;
		b = degreeV + 1;

		while (b < m)
		{
			int i = b;
			while (b < m && MathUtils::IsLessThanOrEqual(knotVectorV[b + 1], knotVectorV[b]))
			{
				b++;
			}
			int multi = b - i + 1;
			if (multi < degreeV)
			{
				double numrator = knotVectorV[b] - knotVectorV[a];
				for (int j = degreeV; j > multi; j--)
				{
					alphaVector[j - multi - 1] = numrator / (knotVectorV[a + j] - knotVectorV[a]);
				}
				int r = degreeV - multi;
				for (int j = 1; j <= r; j++)
				{
					int save = r - j;
					int s = multi + j;
					for (int k = degreeV; k >= s; k--)
					{
						double alpha = alphaVector[k - s];
						for (int row = 0; row <= degreeU; row++)
						{
							decomposedControlPoints[nb][row][k] = alpha * decomposedControlPoints[nb][row][k] + (1.0 - alpha) * decomposedControlPoints[nb][row][k - 1];
						}
					}
				}
			}
			nb++;
			if (b < m)
			{
				for (int i = degreeV - multi; i <= degreeV; i++)
				{
					for (int row = 0; row <= degreeU; row++)
					{
						decomposedControlPoints[nb][row][i] = temp[np][row][b - degreeV + i];
					}
				}
				a = b;
				b++;
			}
		}
	}
	return decomposedControlPoints;
}

void LNLib::NurbsSurface::RemoveKnot(const LN_Surface& surface, double removeKnot, int times, bool isUDirection, LN_Surface& result)
{
	int degreeU = surface.DegreeU;
	int degreeV = surface.DegreeV;
	std::vector<double> knotVectorU = surface.KnotVectorU;
	std::vector<double> knotVectorV = surface.KnotVectorV;
	std::vector<std::vector<XYZW>> controlPoints = surface.ControlPoints;

	VALIDATE_ARGUMENT(degreeU > 0, "degreeU", "Degree must greater than zero.");
	VALIDATE_ARGUMENT(degreeV > 0, "degreeV", "Degree must greater than zero.");
	VALIDATE_ARGUMENT(knotVectorU.size() > 0, "knotVectorU", "KnotVector size must greater than zero.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidKnotVector(knotVectorU), "knotVectorU", "KnotVector must be a nondecreasing sequence of real numbers.");
	VALIDATE_ARGUMENT(knotVectorV.size() > 0, "knotVectorV", "KnotVector size must greater than zero.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidKnotVector(knotVectorV), "knotVectorV", "KnotVector must be a nondecreasing sequence of real numbers.");
	VALIDATE_ARGUMENT(controlPoints.size() > 0, "controlPoints", "ControlPoints must contains one point at least.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidNurbs(degreeU, knotVectorU.size(), controlPoints.size()), "controlPoints", "Arguments must fit: m = n + p + 1");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidNurbs(degreeV, knotVectorV.size(), controlPoints[0].size()), "controlPoints", "Arguments must fit: m = n + p + 1");
	if (isUDirection)
	{
		VALIDATE_ARGUMENT_RANGE(removeKnot, knotVectorU[0], knotVectorU[knotVectorU.size() - 1]);
	}
	else
	{
		VALIDATE_ARGUMENT_RANGE(removeKnot, knotVectorV[0], knotVectorV[knotVectorV.size() - 1]);
	}
	VALIDATE_ARGUMENT(times > 0, "times", "Times must greater than zero.");

	std::vector<double> tempKnotVector;
	std::vector<std::vector<XYZW>> tempControlPoints;

	result = surface;

	if (isUDirection)
	{
		std::vector<std::vector<XYZW>> transposedControlPoints;
		MathUtils::Transpose(controlPoints, transposedControlPoints);

		for (int i = 0; i < transposedControlPoints.size(); i++)
		{
			LN_Curve tc;
			tc.Degree = degreeU;
			tc.KnotVector = knotVectorU;
			tc.ControlPoints = transposedControlPoints[i];

			LN_Curve newtc;

			NurbsCurve::RemoveKnot(tc, removeKnot, times, newtc);
			tempControlPoints.emplace_back(newtc.ControlPoints);
			tempKnotVector = newtc.KnotVector;
		}
		result.KnotVectorU = tempKnotVector;
		result.KnotVectorV = knotVectorV;
		std::vector<std::vector<XYZW>> updatedControlPoints;
		MathUtils::Transpose(tempControlPoints, updatedControlPoints);
		result.ControlPoints = updatedControlPoints;
	}
	else
	{
		for (int i = 0; i < controlPoints.size(); i++)
		{
			LN_Curve tc;
			tc.Degree = degreeU;
			tc.KnotVector = knotVectorU;
			tc.ControlPoints = controlPoints[i];

			LN_Curve newtc;
			NurbsCurve::RemoveKnot(tc, removeKnot, times, newtc);
			tempControlPoints.emplace_back(newtc.ControlPoints);
			tempKnotVector = newtc.KnotVector;
		}
		result.KnotVectorU = knotVectorU;
		result.KnotVectorV = tempKnotVector;
		result.ControlPoints = tempControlPoints;
	}
}

void LNLib::NurbsSurface::ElevateDegree(const LN_Surface& surface, int times, bool isUDirection, LN_Surface& result)
{
	int degreeU = surface.DegreeU;
	int degreeV = surface.DegreeV;
	std::vector<double> knotVectorU = surface.KnotVectorU;
	std::vector<double> knotVectorV = surface.KnotVectorV;
	std::vector<std::vector<XYZW>> controlPoints = surface.ControlPoints;

	VALIDATE_ARGUMENT(degreeU > 0, "degreeU", "Degree must greater than zero.");
	VALIDATE_ARGUMENT(degreeV > 0, "degreeV", "Degree must greater than zero.");
	VALIDATE_ARGUMENT(knotVectorU.size() > 0, "knotVectorU", "KnotVector size must greater than zero.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidKnotVector(knotVectorU), "knotVectorU", "KnotVector must be a nondecreasing sequence of real numbers.");
	VALIDATE_ARGUMENT(knotVectorV.size() > 0, "knotVectorV", "KnotVector size must greater than zero.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidKnotVector(knotVectorV), "knotVectorV", "KnotVector must be a nondecreasing sequence of real numbers.");
	VALIDATE_ARGUMENT(controlPoints.size() > 0, "controlPoints", "ControlPoints must contains one point at least.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidNurbs(degreeU, knotVectorU.size(), controlPoints.size()), "controlPoints", "Arguments must fit: m = n + p + 1");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidNurbs(degreeV, knotVectorV.size(), controlPoints[0].size()), "controlPoints", "Arguments must fit: m = n + p + 1");
	VALIDATE_ARGUMENT(times > 0, "times", "Times must greater than zero.");

	std::vector<double> tempKnotVector;
	std::vector<std::vector<XYZW>> tempControlPoints;

	if (isUDirection)
	{
		std::vector<std::vector<XYZW>> transposedControlPoints;
		MathUtils::Transpose(controlPoints, transposedControlPoints);

		for (int i = 0; i < transposedControlPoints.size(); i++)
		{
			LN_Curve tc;
			tc.Degree = degreeU;
			tc.KnotVector = knotVectorU;
			tc.ControlPoints = transposedControlPoints[i];

			LN_Curve newtc;

			NurbsCurve::ElevateDegree(tc, times, newtc);
			tempControlPoints.emplace_back(newtc.ControlPoints);
			tempKnotVector = newtc.KnotVector;
		}
		result.DegreeU = degreeU + times;
		result.DegreeV = degreeV;
		result.KnotVectorU = tempKnotVector;
		result.KnotVectorV = knotVectorV;
		std::vector<std::vector<XYZW>> updatedControlPoints;
		MathUtils::Transpose(tempControlPoints, updatedControlPoints);
		result.ControlPoints = updatedControlPoints;
	}
	else
	{
		for (int i = 0; i < controlPoints.size(); i++)
		{
			LN_Curve tc;
			tc.Degree = degreeU;
			tc.KnotVector = knotVectorU;
			tc.ControlPoints = controlPoints[i];

			LN_Curve newtc;

			NurbsCurve::ElevateDegree(tc, times, newtc);
			tempControlPoints.emplace_back(newtc.ControlPoints);
			tempKnotVector = newtc.KnotVector;
		}
		result.DegreeU = degreeU;
		result.DegreeV = degreeV + times;
		result.KnotVectorU = knotVectorU;
		result.KnotVectorV = tempKnotVector;
		result.ControlPoints = tempControlPoints;
	}
}

bool LNLib::NurbsSurface::ReduceDegree(const LN_Surface& surface, bool isUDirection, LN_Surface& result)
{
	int degreeU = surface.DegreeU;
	int degreeV = surface.DegreeV;
	std::vector<double> knotVectorU = surface.KnotVectorU;
	std::vector<double> knotVectorV = surface.KnotVectorV;
	std::vector<std::vector<XYZW>> controlPoints = surface.ControlPoints;

	VALIDATE_ARGUMENT(degreeU > 0, "degreeU", "Degree must greater than zero.");
	VALIDATE_ARGUMENT(degreeV > 0, "degreeV", "Degree must greater than zero.");
	VALIDATE_ARGUMENT(knotVectorU.size() > 0, "knotVectorU", "KnotVector size must greater than zero.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidKnotVector(knotVectorU), "knotVectorU", "KnotVector must be a nondecreasing sequence of real numbers.");
	VALIDATE_ARGUMENT(knotVectorV.size() > 0, "knotVectorV", "KnotVector size must greater than zero.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidKnotVector(knotVectorV), "knotVectorV", "KnotVector must be a nondecreasing sequence of real numbers.");
	VALIDATE_ARGUMENT(controlPoints.size() > 0, "controlPoints", "ControlPoints must contains one point at least.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidNurbs(degreeU, knotVectorU.size(), controlPoints.size()), "controlPoints", "Arguments must fit: m = n + p + 1");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidNurbs(degreeV, knotVectorV.size(), controlPoints[0].size()), "controlPoints", "Arguments must fit: m = n + p + 1");

	std::vector<double> tempKnotVector;
	std::vector<std::vector<XYZW>> tempControlPoints;

	if (isUDirection)
	{
		std::vector<std::vector<XYZW>> transposedControlPoints;
		MathUtils::Transpose(controlPoints, transposedControlPoints);

		for (int i = 0; i < transposedControlPoints.size(); i++)
		{
			LN_Curve tc;
			tc.Degree = degreeU;
			tc.KnotVector = knotVectorU;
			tc.ControlPoints = transposedControlPoints[i];

			LN_Curve newtc;

			bool result = NurbsCurve::ReduceDegree(tc, newtc);
			if (!result)
			{
				return false;
			}	
			tempControlPoints.emplace_back(newtc.ControlPoints);
			tempKnotVector = newtc.KnotVector;
		}
		result.DegreeU = degreeU - 1;
		result.DegreeV = degreeV;
		result.KnotVectorU = tempKnotVector;
		result.KnotVectorV = knotVectorV;
		std::vector<std::vector<XYZW>> updatedControlPoints;
		MathUtils::Transpose(tempControlPoints, updatedControlPoints);
		result.ControlPoints = updatedControlPoints;
	}
	else
	{
		std::vector<XYZW> temp;
		for (int i = 0; i < controlPoints.size(); i++)
		{
			LN_Curve tc;
			tc.Degree = degreeU;
			tc.KnotVector = knotVectorU;
			tc.ControlPoints = controlPoints[i];

			LN_Curve newtc;

			bool result = NurbsCurve::ReduceDegree(tc, newtc);
			if (!result)
			{
				return false;
			}	
			tempControlPoints.emplace_back(newtc.ControlPoints);
			tempKnotVector = newtc.KnotVector;
		}
		result.DegreeU = degreeU;
		result.DegreeV = degreeV - 1;
		result.KnotVectorU = knotVectorU;
		result.KnotVectorV = tempKnotVector;
		result.ControlPoints = tempControlPoints;
	}
	return true;
}

void LNLib::NurbsSurface::EquallyTessellate(const LN_Surface& surface, std::vector<XYZ>& tessellatedPoints, std::vector<UV>& correspondingKnots)
{
	int degreeU = surface.DegreeU;
	int degreeV = surface.DegreeV;
	std::vector<double> knotVectorU = surface.KnotVectorU;
	std::vector<double> knotVectorV = surface.KnotVectorV;
	std::vector<std::vector<XYZW>> controlPoints = surface.ControlPoints;

	VALIDATE_ARGUMENT(degreeU > 0, "degreeU", "Degree must greater than zero.");
	VALIDATE_ARGUMENT(degreeV > 0, "degreeV", "Degree must greater than zero.");
	VALIDATE_ARGUMENT(knotVectorU.size() > 0, "knotVectorU", "KnotVector size must greater than zero.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidKnotVector(knotVectorU), "knotVectorU", "KnotVector must be a nondecreasing sequence of real numbers.");
	VALIDATE_ARGUMENT(knotVectorV.size() > 0, "knotVectorV", "KnotVector size must greater than zero.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidKnotVector(knotVectorV), "knotVectorV", "KnotVector must be a nondecreasing sequence of real numbers.");
	VALIDATE_ARGUMENT(controlPoints.size() > 0, "controlPoints", "ControlPoints must contains one point at least.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidNurbs(degreeU, knotVectorU.size(), controlPoints.size()), "controlPoints", "Arguments must fit: m = n + p + 1");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidNurbs(degreeV, knotVectorV.size(), controlPoints[0].size()), "controlPoints", "Arguments must fit: m = n + p + 1");

	std::vector<double> uniqueKvU = knotVectorU;
	uniqueKvU.erase(unique(uniqueKvU.begin(), uniqueKvU.end()), uniqueKvU.end());
	int sizeU = uniqueKvU.size();

	std::vector<double> uniqueKvV = knotVectorV;
	uniqueKvV.erase(unique(uniqueKvV.begin(), uniqueKvV.end()), uniqueKvV.end());
	int sizeV = uniqueKvV.size();

	std::vector<double> tessellatedU;
	int intervals = 100;
	for (int i = 0; i < sizeU - 1; i++)
	{
		double currentU = uniqueKvU[i];
		double nextU = uniqueKvU[i + 1];
		double stepU = (nextU - currentU) / intervals;
		for (int j = 0; j < intervals; j++)
		{
			double u = currentU + stepU * j;
			tessellatedU.emplace_back(u);
		}
	}
	std::vector<double> tessellatedV;
	for (int i = 0; i < sizeV - 1; i++)
	{
		double currentV = uniqueKvV[i];
		double nextV = uniqueKvV[i + 1];
		double stepV = (nextV - currentV) / intervals;
		for (int j = 0; j < intervals; j++)
		{
			double v = currentV + stepV * j;
			tessellatedV.emplace_back(v);
		}
	}

	for (int i = 0; i < tessellatedU.size(); i++)
	{
		for (int j = 0; j < tessellatedV.size(); j++)
		{
			UV uv = UV(tessellatedU[i],tessellatedV[j]);
			correspondingKnots.emplace_back(uv);
			tessellatedPoints.emplace_back(GetPointOnSurface(surface,uv));
		}
	}

	correspondingKnots.emplace_back(UV(knotVectorU[knotVectorU.size() - 1], knotVectorV[knotVectorV.size() - 1]));
	tessellatedPoints.emplace_back(const_cast<XYZW&>(controlPoints[controlPoints.size() - 1][controlPoints[0].size()-1]).ToXYZ(true));
}

LNLib::UV LNLib::NurbsSurface::GetParamOnSurface(const LN_Surface& surface, const XYZ& givenPoint)
{
	int degreeU = surface.DegreeU;
	int degreeV = surface.DegreeV;
	std::vector<double> knotVectorU = surface.KnotVectorU;
	std::vector<double> knotVectorV = surface.KnotVectorV;
	std::vector<std::vector<XYZW>> controlPoints = surface.ControlPoints;

	VALIDATE_ARGUMENT(degreeU > 0, "degreeU", "Degree must greater than zero.");
	VALIDATE_ARGUMENT(degreeV > 0, "degreeV", "Degree must greater than zero.");
	VALIDATE_ARGUMENT(knotVectorU.size() > 0, "knotVectorU", "KnotVector size must greater than zero.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidKnotVector(knotVectorU), "knotVectorU", "KnotVector must be a nondecreasing sequence of real numbers.");
	VALIDATE_ARGUMENT(knotVectorV.size() > 0, "knotVectorV", "KnotVector size must greater than zero.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidKnotVector(knotVectorV), "knotVectorV", "KnotVector must be a nondecreasing sequence of real numbers.");
	VALIDATE_ARGUMENT(controlPoints.size() > 0, "controlPoints", "ControlPoints must contains one point at least.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidNurbs(degreeU, knotVectorU.size(), controlPoints.size()), "controlPoints", "Arguments must fit: m = n + p + 1");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidNurbs(degreeV, knotVectorV.size(), controlPoints[0].size()), "controlPoints", "Arguments must fit: m = n + p + 1");

	double minValue = Constants::MaxDistance;

	int maxIterations = 10;
	UV param = UV(Constants::DoubleEpsilon, Constants::DoubleEpsilon);

	double minUParam = knotVectorU[0];
	double maxUParam = knotVectorU[knotVectorU.size() - 1];
	double minVParam = knotVectorV[0];
	double maxVParam = knotVectorV[knotVectorV.size() - 1];

	double a = minUParam;
	double b = maxUParam;
	double c = minVParam;
	double d = maxVParam;

	bool isClosedU = ValidationUtils::IsClosedU(controlPoints);
	bool isClosedV = ValidationUtils::IsClosedV(controlPoints);

	std::vector<XYZ> tessellatedPoints;
	std::vector<UV> correspondingKnots;
	EquallyTessellate(surface, tessellatedPoints, correspondingKnots);
	for (int i = 0; i < tessellatedPoints.size() - 1; i++)
	{
		UV currentUV = correspondingKnots[i];
		UV nextUV = correspondingKnots[i + 1];

		XYZ currentPoint = tessellatedPoints[i];
		XYZ nextPoint = tessellatedPoints[i + 1];

		XYZ vector1 = currentPoint - givenPoint;
		XYZ vector2 = nextPoint - currentPoint;
		double dot = vector1.DotProduct(vector2);

		XYZ projectPoint;
		UV project = UV(Constants::DoubleEpsilon, Constants::DoubleEpsilon);

		if (dot < 0)
		{
			projectPoint = currentPoint;
			project = currentUV;
		}
		else if (dot > 1)
		{
			projectPoint = nextPoint;
			project = nextUV;
		}
		else
		{
			projectPoint = currentPoint + dot * vector1.Normalize();
			project = currentUV + (nextUV - currentUV) * dot;
		}

		double distance = (givenPoint - projectPoint).Length();
		if (distance < minValue)
		{
			minValue = distance;
			param = project;
		}
	}

	int counters = 0;
	while (counters < maxIterations)
	{
		std::vector<std::vector<XYZ>> derivatives = ComputeRationalSurfaceDerivatives(surface, 2, param);
		XYZ difference = derivatives[0][0] - givenPoint;
		double fa = derivatives[1][0].DotProduct(difference);
		double fb = derivatives[0][1].DotProduct(difference);

		double condition1 = difference.Length();
		double condition2a = std::abs(fa / (derivatives[1][0].Length() * condition1));
		double condition2b = std::abs(fb / (derivatives[0][1].Length() * condition1));

		if (condition1 < Constants::DistanceEpsilon &&
			condition2a < Constants::DistanceEpsilon &&
			condition2b < Constants::DistanceEpsilon)
		{
			return param;
		}

		XYZ Su = derivatives[1][0];
		XYZ Sv = derivatives[0][1];

		XYZ Suu = derivatives[2][0];
		XYZ Svv = derivatives[0][2];

		XYZ Suv,Svu = derivatives[1][1];

		double fuv = -Su.DotProduct(difference);
		double guv = -Sv.DotProduct(difference);

		double fu = Su.DotProduct(Su) + difference.DotProduct(Suu);
		double fv = Su.DotProduct(Sv) + difference.DotProduct(Suv);
		double gu = Su.DotProduct(Sv) + difference.DotProduct(Svu);
		double gv = Sv.DotProduct(Sv) + difference.DotProduct(Svv);
		
		if (MathUtils::IsAlmostEqualTo(fu * gv, fv * gu))
		{
			counters++;
			continue;
		}

		double deltaU = ((-fuv * gv) - fv * (-guv)) / (fu * gv - fv * gu);
		double deltaV = (fu * (-guv) - (-fuv) * gu)/ (fu * gv - fv * gu);

		UV delta = UV(deltaU,deltaV);
		UV temp = param + delta;

		if (!isClosedU)
		{
			if (param[0] < a)
			{
				param = UV(a, param[1]);
			}
			if (param[0] > b)
			{
				param = UV(b, param[1]);
			}
		}
		if (!isClosedV)
		{
			if (param[1] < c)
			{
				param = UV(param[0], c);
			}
			if (param[1] > d)
			{
				param = UV(param[0], d);
			}
		}
		if (isClosedU)
		{
			if (param[0] < a)
			{
				param = UV(b-(a-param[0]), param[1]);
			}
			if (param[0] > b)
			{
				param = UV(a+(param[0]-b), param[1]);
			}
		}
		if (isClosedV)
		{
			if (param[1] < c)
			{
				param = UV(param[0], d-(c-param[1]));
			}
			if (param[1] > d)
			{
				param = UV(param[0], c+(param[1]-d));
			}
		}

		double condition4a = ((temp[0] - param[0]) * derivatives[1][0]).Length();
		double condition4b = ((temp[1] - param[1]) * derivatives[0][1]).Length();
		if (condition4a + condition4b < Constants::DistanceEpsilon) {
			return param;
		}

		param = temp;
		counters++;
	}
	return param;
}

bool LNLib::NurbsSurface::GetUVTangent(const LN_Surface& surface, const UV param, const XYZ& tangent, UV& uvTangent)
{
	int degreeU = surface.DegreeU;
	int degreeV = surface.DegreeV;
	std::vector<double> knotVectorU = surface.KnotVectorU;
	std::vector<double> knotVectorV = surface.KnotVectorV;
	std::vector<std::vector<XYZW>> controlPoints = surface.ControlPoints;

	VALIDATE_ARGUMENT(degreeU > 0, "degreeU", "Degree must greater than zero.");
	VALIDATE_ARGUMENT(degreeV > 0, "degreeV", "Degree must greater than zero.");
	VALIDATE_ARGUMENT(knotVectorU.size() > 0, "knotVectorU", "KnotVector size must greater than zero.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidKnotVector(knotVectorU), "knotVectorU", "KnotVector must be a nondecreasing sequence of real numbers.");
	VALIDATE_ARGUMENT(knotVectorV.size() > 0, "knotVectorV", "KnotVector size must greater than zero.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidKnotVector(knotVectorV), "knotVectorV", "KnotVector must be a nondecreasing sequence of real numbers.");
	VALIDATE_ARGUMENT(controlPoints.size() > 0, "controlPoints", "ControlPoints must contains one point at least.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidNurbs(degreeU, knotVectorU.size(), controlPoints.size()), "controlPoints", "Arguments must fit: m = n + p + 1");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidNurbs(degreeV, knotVectorV.size(), controlPoints[0].size()), "controlPoints", "Arguments must fit: m = n + p + 1");

	std::vector<std::vector<XYZ>> derivatives = ComputeRationalSurfaceDerivatives(surface, 1, param);
	XYZ Su = derivatives[1][0];
	XYZ Sv = derivatives[0][1];

	double a = Su.DotProduct(Su);
	double b = Su.DotProduct(Sv);

	double c = Su.DotProduct(Sv);
	double d = Sv.DotProduct(Sv);

	double e = Su.DotProduct(tangent);
	double f = Sv.DotProduct(tangent);

	if (MathUtils::IsAlmostEqualTo(a * d, b * c))
	{
		return false;
	}
	
	double u = (e * d - b * f) / (a * d - b * c);
	double v = (a * f - e * c) / (a * d - b * c);

	uvTangent = UV(u, v);
	return true;
}

void LNLib::NurbsSurface::CreateBilinearSurface(const XYZ& point1, const XYZ& point2, const XYZ& point3, const XYZ& point4, LN_Surface& surface)
{
	int degree = 3;
	surface.DegreeU = surface.DegreeV = degree;

	std::vector<double> knotVectorU;
	std::vector<double> knotVectorV;
	std::vector<std::vector<XYZW>> controlPoints;

	for (int i = 0; i <= degree; i++) 
	{
		std::vector<XYZW> row;
		for (int j = 0; j <= degree; j++) 
		{
			double l = 1.0 - i / (double)degree;
			XYZ inter12 = l * point1 + (1 - l) * point2;
			XYZ inter43 = l * point4 + (1 - l) * point3;

			XYZ res = inter12 * (1- j / (double)degree) + (j / (double)degree) * inter43;
			row.emplace_back(XYZW(res, 1.0));
		}
		controlPoints.emplace_back(row);

		knotVectorU.insert(knotVectorU.begin(),0.0);
		knotVectorU.emplace_back(1.0);

		knotVectorV.insert(knotVectorV.begin(), 0.0);
		knotVectorV.emplace_back(1.0);
	}
	surface.KnotVectorU = knotVectorU;
	surface.KnotVectorV = knotVectorV;
	surface.ControlPoints = controlPoints;
}

bool LNLib::NurbsSurface::CreateCylindricalSurface(const XYZ& origin, const XYZ& xAxis, const XYZ& yAxis, double startRad, double endRad, double radius, double height, LN_Surface& surface)
{
	VALIDATE_ARGUMENT(!xAxis.IsZero(), "xAxis", "xAxis must not be zero vector.");
	VALIDATE_ARGUMENT(!yAxis.IsZero(), "yAxis", "yAxis must not be zero vector.")
	VALIDATE_ARGUMENT(MathUtils::IsGreaterThan(endRad, startRad), "endRad", "endRad must greater than startRad.");
	VALIDATE_ARGUMENT(MathUtils::IsGreaterThan(radius, 0.0), "radius", "Radius must greater than zero.");
	VALIDATE_ARGUMENT(MathUtils::IsGreaterThan(height, 0.0), "height", "Height must greater than zero.");

	XYZ nX = const_cast<XYZ&>(xAxis).Normalize();
	XYZ nY = const_cast<XYZ&>(yAxis).Normalize();

	LN_Curve arc;
	bool isCreated = NurbsCurve::CreateArc(origin, nX, nY, startRad, endRad, radius, radius, arc);
	if (!isCreated) return false;

	XYZ axis = nX.CrossProduct(nY);
	XYZ translation = height * axis;
	XYZ halfTranslation = 0.5 * height * axis;

	int size = arc.ControlPoints.size();
	std::vector<std::vector<XYZW>> controlPoints(3,std::vector<XYZW>(size));

	for (int i = 0; i < size; i++)
	{
		double w = arc.ControlPoints[i].GetW();

		controlPoints[2][i] = XYZW(arc.ControlPoints[i].ToXYZ(true), w);
		controlPoints[1][i] = XYZW(halfTranslation + arc.ControlPoints[i].ToXYZ(true), w);
		controlPoints[0][i] = XYZW(translation + arc.ControlPoints[i].ToXYZ(true), w);
	}

	surface.DegreeU = 2;
	surface.DegreeV = arc.Degree;
	surface.KnotVectorU = { 0,0,0,1,1,1 };
	surface.KnotVectorV = arc.KnotVector;
	surface.ControlPoints = controlPoints;
	return true;
}

void LNLib::NurbsSurface::CreateRuledSurface(const LN_Curve& curve0, const LN_Curve& curve1, LN_Surface& surface)
{
	int degree0 = curve0.Degree;
	std::vector<double> knotVector0 = curve0.KnotVector;
	std::vector<XYZW> controlPoints0 = curve0.ControlPoints;

	VALIDATE_ARGUMENT(degree0 > 0, "degree0", "Degree must greater than zero.");
	VALIDATE_ARGUMENT(knotVector0.size() > 0, "knotVector0", "KnotVector size must greater than zero.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidKnotVector(knotVector0), "knotVector0", "KnotVector must be a nondecreasing sequence of real numbers.");
	VALIDATE_ARGUMENT(controlPoints0.size() > 0, "controlPoints0", "ControlPoints must contains one point at least.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidNurbs(degree0, knotVector0.size(), controlPoints0.size()), "controlPoints0", "Arguments must fit: m = n + p + 1");

	int degree1 = curve1.Degree;
	std::vector<double> knotVector1 = curve1.KnotVector;
	std::vector<XYZW> controlPoints1 = curve1.ControlPoints;

	VALIDATE_ARGUMENT(degree1 > 0, "degree1", "Degree must greater than zero.");
	VALIDATE_ARGUMENT(knotVector1.size() > 0, "knotVector1", "KnotVector size must greater than zero.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidKnotVector(knotVector1), "knotVector1", "KnotVector must be a nondecreasing sequence of real numbers.");
	VALIDATE_ARGUMENT(controlPoints1.size() > 0, "controlPoints1", "ControlPoints must contains one point at least.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidNurbs(degree1, knotVector1.size(), controlPoints1.size()), "controlPoints1", "Arguments must fit: m = n + p + 1");

	int k0Size = knotVector0.size();
	int k1Size = knotVector1.size();
	bool knotVectorCheck = MathUtils::IsAlmostEqualTo(knotVector0[0], knotVector1[0]) && MathUtils::IsAlmostEqualTo(knotVector0[k0Size - 1], knotVector1[k1Size - 1]);
	VALIDATE_ARGUMENT(knotVectorCheck, "knotVector0 & knotVector1", "Ensure that the two curves are defined on the same parameter range.");

	surface.DegreeU = std::max(degree0, degree1);
	surface.DegreeV = 1;

	std::vector<double> kv0 = knotVector0;
	std::vector<XYZW> cp0 = controlPoints0;
	if (degree0 < surface.DegreeU)
	{
		int times = surface.DegreeU - degree0;
		LN_Curve tc;
		NurbsCurve::ElevateDegree(curve0, times, tc);
		kv0 = tc.KnotVector;
		cp0 = tc.ControlPoints;
	}

	std::vector<double> kv1 = knotVector1;
	std::vector<XYZW> cp1 = controlPoints1;
	if (degree1 < surface.DegreeU)
	{
		int times = surface.DegreeU - degree1;
		LN_Curve tc;
		NurbsCurve::ElevateDegree(curve1, times, tc);
		kv1 = tc.KnotVector;
		cp1 = tc.ControlPoints;
	}
	
	surface.KnotVectorU = kv0;
	if (kv0 != kv1)
	{
		std::vector<double> insertedKnotElement0;
		std::vector<double> insertedKnotElement1;
		KnotVectorUtils::GetInsertedKnotElement(kv0, kv1, insertedKnotElement0, insertedKnotElement1);

		if (insertedKnotElement0.size() > 0)
		{
			LN_Curve tc;
			tc.Degree = surface.DegreeU;
			tc.KnotVector = kv0;
			tc.ControlPoints = cp0;

			LN_Curve newtc;
			NurbsCurve::RefineKnotVector(tc, insertedKnotElement0, newtc);
			surface.KnotVectorU = newtc.KnotVector;
			cp0 = newtc.ControlPoints;
		}
		if (insertedKnotElement1.size() > 0)
		{
			LN_Curve tc;
			tc.Degree = surface.DegreeU;
			tc.KnotVector = kv1;
			tc.ControlPoints = cp1;

			LN_Curve newtc;
			NurbsCurve::RefineKnotVector(tc, insertedKnotElement1, newtc);
			cp1 = newtc.ControlPoints;
		}
	}

	surface.KnotVectorV = { 0,0,1,1 };
	int size = cp0.size();
	std::vector<std::vector<XYZW>> controlPoints(size, std::vector<XYZW>(2));
	for (int i = 0; i < size; i++)
	{
		controlPoints[i][0] = cp0[i];
		controlPoints[i][1] = cp1[i];
	}
	surface.ControlPoints = controlPoints;
}

bool LNLib::NurbsSurface::CreateRevolvedSurface(const XYZ& origin, const XYZ& axis, double rad, const LN_Curve& profile, LN_Surface& surface)
{
	VALIDATE_ARGUMENT(!axis.IsZero(), "axis", "Axis must not be zero vector.");
	std::vector<double> knotVectorU;

	int narcs = 0;
	if (MathUtils::IsLessThanOrEqual(rad, Constants::Pi / 2.0))
	{
		narcs = 1;
		knotVectorU.resize(2 * narcs + 3 + 1);
	}
	else
	{
		if (MathUtils::IsLessThanOrEqual(rad, Constants::Pi))
		{
			narcs = 2;
			knotVectorU.resize(2 * narcs + 3 + 1);
			knotVectorU[3] = knotVectorU[4] = 0.5;
		}
		else if (MathUtils::IsLessThanOrEqual(rad, 3 * Constants::Pi / 2))
		{
			narcs = 3;
			knotVectorU.resize(2 * narcs + 3 + 1);
			knotVectorU[3] = knotVectorU[4] = 1.0 / 3.0;
			knotVectorU[5] = knotVectorU[6] = 2.0 / 3.0;
		}
		else
		{
			narcs = 4;
			knotVectorU.resize(2 * narcs + 3 + 1);
			knotVectorU[3] = knotVectorU[4] = 0.25;
			knotVectorU[5] = knotVectorU[6] = 0.5;
			knotVectorU[7] = knotVectorU[8] = 0.75;
		}
	}

	double dtheta = rad / narcs;
	int j = 3 + 2 * (narcs - 1);
	for (int i = 0; i < 3; i++)
	{
		knotVectorU[i] = 0.0;
		knotVectorU[j + i] = 1.0;
	}

	int n = 2 * narcs;
	double wm = cos(dtheta / 2.0);
	double angle = 0.0;
	std::vector<double> cosines(narcs + 1,0.0);
	std::vector<double> sines(narcs + 1,0.0);

	for (int i = 1; i <= narcs; i++)
	{
		angle += dtheta;
		cosines[i] = cos(angle);
		sines[i] = sin(angle);
	}

	int m = profile.ControlPoints.size() - 1;
	XYZ X, Y, O, P0, P2, T0, T2;
	double r = 0.0;
	int index = 0;

	int degreeU = 2;
	std::vector<std::vector<XYZW>> controlPoints(n + 1, std::vector<XYZW>(m + 1));

	for (int j = 0; j <= m; j++)
	{
		XYZW gp = profile.ControlPoints[j];
		XYZ p = gp.ToXYZ(true);

		XYZ O = Projection::PointToRay(origin, axis, p);
		X = p - O;

		r = X.Normalize().Length();
		Y = axis.CrossProduct(X);

		if (MathUtils::IsGreaterThan(r, 0.0))
		{
			X = X / r;
			Y = Y / r;
		}

		P0 = p;
		controlPoints[0][j] = XYZW(P0, gp[3]);

		T0 = Y;
		index = 0;
		angle = 0.0;

		for (int i = 1; i <= narcs; i++)
		{
			P2 = MathUtils::IsAlmostEqualTo(r,0.0)? O : O + r * cosines[i] * X + r * sines[i] * Y;
			controlPoints[index + 2][j] = XYZW(P2, gp[3]);
			T2 = -sines[i] * X + cosines[i] * Y;

			if (MathUtils::IsAlmostEqualTo(r, 0.0))
			{
				controlPoints[index + 1][j] = XYZW(O, wm * gp[3]);
			}
			else
			{
				XYZ intersectPoint;
				double param0;
				double param1;
				CurveCurveIntersectionType type = Intersection::ComputeRays(P0, T0, P2, T2, param0, param1, intersectPoint);
				if (type != CurveCurveIntersectionType::Intersecting)
				{
					return false;
				}
				controlPoints[index + 1][j] = XYZW(intersectPoint, wm * gp[3]);
			}
			
			index += 2;
			if (i < narcs)
			{
				P0 = P2;
				T0 = T2;
			}
		}
	}
	surface.DegreeU = degreeU;
	surface.DegreeV = profile.Degree;
	surface.KnotVectorU = knotVectorU;
	surface.KnotVectorV = profile.KnotVector;
	surface.ControlPoints = controlPoints;
	return true;
}

std::vector<std::vector<LNLib::XYZW>> LNLib::NurbsSurface::NonuniformScaling(const std::vector<std::vector<XYZW>>& controlPoints, double xFactor, double yFactor, double zFactor, const XYZ& referencePoint)
{
	std::vector<std::vector<XYZW>> result;
	int row = controlPoints.size();
	int col = controlPoints[0].size();

	if (referencePoint.IsAlmostEqualTo(XYZ(0, 0, 0)))
	{
		for (int i = 0; i < row; i++)
		{
			std::vector<XYZW> row;
			for (int j = 0; j < col; j++)
			{
				XYZW current = controlPoints[i][j];
				double weight = current.GetW();
				double x = current.GetWX() / weight;
				double y = current.GetWY() / weight;
				double z = current.GetWZ() / weight;

				row.emplace_back(XYZW(x * xFactor, y * yFactor, z * zFactor, weight));
			}
			result.emplace_back(row);
		}
	}
	else
	{
		for (int i = 0; i < row; i++)
		{
			std::vector<XYZW> row;
			for (int j = 0; j < col; j++)
			{
				XYZW current = controlPoints[i][j];
				double weight = current.GetW();
				double x = current.GetWX() / weight;
				double y = current.GetWY() / weight;
				double z = current.GetWZ() / weight;

				double newX = x * xFactor + (1 - xFactor) * referencePoint[0];
				double newY = y * yFactor + (1 - yFactor) * referencePoint[1];
				double newZ = z * zFactor + (1 - zFactor) * referencePoint[2];

				row.emplace_back(XYZW(newX, newY, newZ, weight));
			}
			result.emplace_back(row);
		}
	}
	return result;
}

void LNLib::NurbsSurface::GlobalInterpolation(const std::vector<std::vector<XYZ>>& throughPoints, int degreeU, int degreeV, LN_Surface& surface)
{
	VALIDATE_ARGUMENT(throughPoints.size() > 0, "throughPoints", "ThroughPoints row size must greater than zero.");
	VALIDATE_ARGUMENT(throughPoints[0].size() > 0, "throughPoints", "ThroughPoints column size must greater than zero.");
	VALIDATE_ARGUMENT(degreeU > 0, "degreeU", "DegreeU must greater than zero.");
	VALIDATE_ARGUMENT(degreeV > 0, "degreeV", "DegreeV must greater than zero.");

	std::vector<double> uk;
	std::vector<double> vl;

	Interpolation::GetSurfaceMeshParameterization(throughPoints, uk, vl);

	int rows = throughPoints.size();
	int cols = throughPoints[0].size();

	std::vector<std::vector<XYZW>> controlPoints(rows, std::vector<XYZW>(cols));
	std::vector<double> knotVectorU = Interpolation::AverageKnotVector(degreeU, uk);
	std::vector<double> knotVectorV = Interpolation::AverageKnotVector(degreeV, vl);

	for (int j = 0; j < cols; j++)
	{
		std::vector<XYZ> temp(rows);
		for (int i = 0; i < rows; i++)
		{
			temp[i] = throughPoints[i][j];
		}
		LN_Curve tc;
		NurbsCurve::GlobalInterpolation(degreeU, temp, tc);
		for (int i = 0; i < rows; i++)
		{
			controlPoints[i][j] = tc.ControlPoints[i];
		}
	}

	for (int i = 0; i < rows; i++)
	{
		std::vector<XYZ> temp(cols);
		for (int j = 0; j < cols; j++)
		{
			temp[j] = throughPoints[i][j];
		}
		LN_Curve tc;
		NurbsCurve::GlobalInterpolation(degreeV, temp, tc);
		for (int j = 0; j < cols; j++)
		{
			controlPoints[i][j] = tc.ControlPoints[j];
		}
	}
	surface.DegreeU = degreeU;
	surface.DegreeV = degreeV;
	surface.KnotVectorU = knotVectorU;
	surface.KnotVectorV = knotVectorV;
	surface.ControlPoints = controlPoints;
}

bool LNLib::NurbsSurface::BicubicLocalInterpolation(const std::vector<std::vector<XYZ>>& throughPoints, LN_Surface& surface)
{
	VALIDATE_ARGUMENT(throughPoints.size() > 0, "throughPoints", "ThroughPoints row size must greater than zero.");
	VALIDATE_ARGUMENT(throughPoints[0].size() > 0, "throughPoints", "ThroughPoints column size must greater than zero.");

	int degreeU = 3;
	int degreeV = 3;

	std::vector<double> knotVectorU;
	std::vector<double> knotVectorV;
	std::vector<std::vector<XYZW>> controlPoints;

	int row = throughPoints.size();
	int n = row - 1;
	int column = throughPoints[0].size();
	int m = column - 1;

	std::vector<std::vector<std::vector<XYZ>>> td(n + 1, std::vector<std::vector<XYZ>>(m+1, std::vector<XYZ>(3,XYZ())));

	std::vector<double> ub(n + 1, 0.0);
	std::vector<double> vb(m + 1, 0.0);

	std::vector<double> r(m + 1);
	std::vector<double> s(n + 1);

	double total = 0.0;
	for (int l = 0; l <= m; l++)
	{
		std::vector<XYZ> columnData = MathUtils::GetColumn(throughPoints, l);
		std::vector<XYZ> tvkl = Interpolation::ComputeTangent(columnData);

		r[l] = 0.0;
		for (int k = 0; k <= n; k++)
		{
			td[k][l][1] = tvkl[k];

			if (k > 0)
			{
				double d = throughPoints[k][l].Distance(throughPoints[k - 1][l]);
				ub[k] += d;
				r[l] += d;
			}
		}
		total += r[l];
	}
	for (int k = 1; k < n; k++)
	{
		ub[k] = ub[k - 1] + ub[k] / total;
	}
	ub[n] = 1.0;
	total = 0.0;

	for (int k = 0; k <= n; k++)
	{
		std::vector<XYZ> tukl = Interpolation::ComputeTangent(throughPoints[k]);
		s[k] = 0.0;
		for (int l = 0; l <= m; l++)
		{
			td[k][l][0] = tukl[l];

			if (l > 0)
			{
				double d = throughPoints[k][l].Distance(throughPoints[k][l - 1]);
				vb[l] += d;
				s[k] += d;
			}
		}
		total += s[k];
	}
	for (int l = 1; l < m; l++)
	{
		vb[l] = vb[l - 1] + vb[l] / total;
	}
	vb[m] = 1.0;
	total = 0.0;

	int kuSize = 2 * ub.size() + 2 + 2;
	knotVectorU.resize(kuSize);
	for (int i= 0; i < 4; i++)
	{
		knotVectorU[i] = knotVectorU[i + 1] = knotVectorU[i + 2] = knotVectorU[i + 3] = 0.0;
		knotVectorU[kuSize -1] = knotVectorU[kuSize - 2] = knotVectorU[kuSize - 3] = knotVectorU[kuSize - 4] = 1.0;
	}
	int ii = 4;
	for (int i = 1; i < ub.size() - 1; i++)
	{
		knotVectorU[ii] = ub[i];
		knotVectorU[ii + 1] = ub[i];
		ii += 2;
	}

	int kvSize = 2 * vb.size() + 2 + 2;
	knotVectorV.resize(kvSize);
	for (int i = 0; i < 4; i++)
	{
		knotVectorV[i] = knotVectorV[i + 1] = knotVectorV[i + 2] = knotVectorV[i + 3] = 0.0;
		knotVectorV[kvSize - 1] = knotVectorV[kvSize - 2] = knotVectorV[kvSize - 3] = knotVectorV[kvSize - 4] = 1.0;
	}
	
	ii = 4;
	for (int i = 1; i < vb.size() - 1; i++)
	{
		knotVectorV[ii] = vb[i];
		knotVectorV[ii + 1] = vb[i];
		ii += 2;
	}
	
	std::vector<std::vector<XYZ>> bezierControlPoints(3* row - 2, std::vector<XYZ>(3 * column - 2,XYZ(0,0,0)));
	for (int i = 0; i < row; i++)
	{
		int n = throughPoints[i].size();
		std::vector<XYZ> T;
		for (int c = 0; c < td[0].size(); c++)
		{
			T.emplace_back(td[i][c][1]);
		}

		std::vector<XYZ> temp(3 * n - 2, XYZ(0,0,0));
		for (int j = 0; j < n; j++)
		{
			temp[3 * i] = throughPoints[i][j];
		}
		for (int j = 0; j < n - 1; j++)
		{
			double a =  (ub[j + 1] - ub[j]) * Interpolation::GetTotalChordLength(throughPoints[i]);
			temp[3 * j + 1] = throughPoints[i][j] + a / 3.0 * T[j];
			temp[3 * j + 2] = throughPoints[i][j+1] - a / 3.0 * T[j];
		}
		bezierControlPoints[3 * i] = temp;
	}
	for (int i = 0; i < column; i++)
	{
		auto columnData = MathUtils::GetColumn(throughPoints, i);
		int n = columnData.size();
		std::vector<XYZ> T;
		for (int r = 0; r < td.size(); r++)
		{
			T.emplace_back(td[r][i][0]);
		}

		std::vector<XYZ> temp(3 * n - 2, XYZ(0,0,0));
		for (int j = 0; j < n; j++)
		{
			temp[3 * i] = columnData[j];
		}
		for (int j = 0; j < n - 1; j++)
		{
			double a = (vb[j + 1] - vb[j]) * Interpolation::GetTotalChordLength(columnData);
			temp[3 * j + 1] = columnData[j] + a / 3.0 * T[j];
			temp[3 * j + 2] = columnData[j + 1] - a / 3.0 * T[j];
		}
		for (int r = 0; r < bezierControlPoints.size(); r++)
		{
			bezierControlPoints[r][3 * i] = temp[r];
		}
	}


	for (int k = 1; k < n; k++)
	{
		double ak = (ub[k] - ub[k - 1]) / ((ub[k] - ub[k - 1]) + (ub[k + 1] - ub[k]));
		for (int l = 1; l < m; l++)
		{
			double bl = bl = (vb[l] - vb[l - 1]) / ((vb[l] - vb[l - 1]) + (vb[l + 1] - vb[l]));
			
			XYZ dvukl = (1 - ak) * (td[k][l][1] - td[k - 1][l][1]) / (ub[k] - ub[k - 1]) + ak * (td[k + 1][l][1] - td[k][l][1]) / (ub[k + 1] - ub[k]);
			XYZ duvkl = (1 - bl) * (td[k][l][0] - td[k][l - 1][0]) / (vb[l] - vb[l - 1]) + bl * (td[k][l + 1][0] - td[k][l][0]) / (vb[l + 1] - vb[l]);

			td[k][l][2] = (ak * duvkl + bl * dvukl) / (ak + bl);
		}
	}

	for (int k = 0; k < n; k++)
	{
		for (int l = 0; l < m; l++)
		{
			double gamma = (ub[k + 1] - ub[k]) * (vb[l + 1] - vb[l]) / 9.0;
			int ii = 3 * k;
			int jj = 3 * l;
			bezierControlPoints[ii + 1][jj + 1] =  gamma * td[k][l][2]         + bezierControlPoints[ii][jj + 1]     + bezierControlPoints[ii + 1][jj]     - bezierControlPoints[ii][jj];
			bezierControlPoints[ii + 2][jj + 1] = -gamma * td[k + 1][l][2]     + bezierControlPoints[ii + 3][jj + 1] - bezierControlPoints[ii + 3][jj]     + bezierControlPoints[ii + 2][jj];
			bezierControlPoints[ii + 1][jj + 2] = -gamma * td[k][l + 1][2]     + bezierControlPoints[ii + 1][jj + 3] - bezierControlPoints[ii][jj + 3]     + bezierControlPoints[ii][jj + 2];
			bezierControlPoints[ii + 2][jj + 2] =  gamma * td[k + 1][l + 1][2] + bezierControlPoints[ii + 2][jj + 3] + bezierControlPoints[ii + 3][jj + 2] - bezierControlPoints[ii + 3][jj + 3];
		} 
	}

	std::vector<std::vector<XYZ>> columnFilter;
	auto ind = GetIndex(column);
	for (int c = 0; c < ind.size(); c++)
	{
		auto columnData = MathUtils::GetColumn(bezierControlPoints, ind[c]);
		columnFilter.emplace_back(columnData);
	}
	std::vector<std::vector<XYZ>> Tcf;
	MathUtils::Transpose(columnFilter, Tcf);
	ind = GetIndex(row);
	std::vector<std::vector<XYZ>> rowFilter;
	for (int r = 0; r < ind.size(); r++)
	{
		rowFilter.emplace_back(Tcf[ind[r]]);
	}
	controlPoints = ControlPointsUtils::ToXYZW(rowFilter);

	surface.DegreeU = degreeU;
	surface.DegreeV = degreeV;
	surface.KnotVectorU = knotVectorU;
	surface.KnotVectorV = knotVectorV;
	surface.ControlPoints = controlPoints;
	return true;
}

bool LNLib::NurbsSurface::GlobalApproximation(const std::vector<std::vector<XYZ>>& throughPoints, int degreeU, int degreeV, int controlPointsRows, int controlPointsColumns, LN_Surface& surface)
{
	VALIDATE_ARGUMENT(throughPoints.size() > 0, "throughPoints", "ThroughPoints row size must greater than zero.");
	VALIDATE_ARGUMENT(throughPoints[0].size() > 0, "throughPoints", "ThroughPoints column size must greater than zero.");
	VALIDATE_ARGUMENT(degreeU > 0, "degreeU", "DegreeU must greater than zero.");
	VALIDATE_ARGUMENT(degreeV > 0, "degreeV", "DegreeV must greater than zero.");

	int rows = controlPointsRows;
	int n = rows - 1;
	int columns = controlPointsColumns;
	int m = columns - 1;

	std::vector<double> knotVectorU;
	std::vector<double> knotVectorV;
	std::vector<std::vector<XYZW>> controlPoints;

	std::vector<std::vector<XYZ>> tempControlPoints;
	for (int i = 0; i < rows; i++)
	{
		LN_Curve tc;
		bool result = NurbsCurve::LeastSquaresApproximation(degreeU, throughPoints[i], rows, tc);
		if (!result) return false;
		std::vector<XYZ> points = ControlPointsUtils::ToXYZ(tc.ControlPoints);
		tempControlPoints.emplace_back(points);
		knotVectorU = tc.KnotVector;
	}

	std::vector<std::vector<XYZ>> preControlPoints;
	std::vector<std::vector<XYZ>> tPoints;
 	for (int i = 0; i < columns; i++)
	{
		std::vector<XYZ> c = MathUtils::GetColumn(tempControlPoints, i);
		LN_Curve tc;
		bool result = NurbsCurve::LeastSquaresApproximation(degreeV, c, columns, tc);
		if (!result) return false;
		std::vector<XYZ> points = ControlPointsUtils::ToXYZ(tc.ControlPoints);
		tPoints.emplace_back(points);
		knotVectorV = tc.KnotVector;
	}
	MathUtils::Transpose(tPoints, preControlPoints);
	controlPoints = ControlPointsUtils::ToXYZW(preControlPoints);

	surface.DegreeU = degreeU;
	surface.DegreeV = degreeV;
	surface.KnotVectorU = knotVectorU;
	surface.KnotVectorV = knotVectorV;
	surface.ControlPoints = controlPoints;
	return true;
}

bool LNLib::NurbsSurface::CreateSwungSurface(const LN_Curve& profile, const LN_Curve& trajectory, double scale, LN_Surface& surface)
{
	int pDegree = profile.Degree;
	std::vector<double> pKnotVector = profile.KnotVector;
	std::vector<XYZW> pControlPoints = profile.ControlPoints;
	VALIDATE_ARGUMENT(pDegree > 0, "degree", "Degree must greater than zero.");
	VALIDATE_ARGUMENT(pKnotVector.size() > 0, "knotVector", "KnotVector size must greater than zero.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidKnotVector(pKnotVector), "knotVector", "KnotVector must be a nondecreasing sequence of real numbers.");
	VALIDATE_ARGUMENT(pControlPoints.size() > 0, "controlPoints", "ControlPoints must contains one point at least.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidNurbs(pDegree, pKnotVector.size(), pControlPoints.size()), "controlPoints", "Arguments must fit: m = n + p + 1");
	
	int tDegree = trajectory.Degree;
	std::vector<double> tKnotVector = trajectory.KnotVector;
	std::vector<XYZW> tControlPoints = trajectory.ControlPoints;
	VALIDATE_ARGUMENT(tDegree > 0, "degree", "Degree must greater than zero.");
	VALIDATE_ARGUMENT(tKnotVector.size() > 0, "knotVector", "KnotVector size must greater than zero.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidKnotVector(tKnotVector), "knotVector", "KnotVector must be a nondecreasing sequence of real numbers.");
	VALIDATE_ARGUMENT(tControlPoints.size() > 0, "controlPoints", "ControlPoints must contains one point at least.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidNurbs(tDegree, tKnotVector.size(), tControlPoints.size()), "controlPoints", "Arguments must fit: m = n + p + 1");
	
	int size = tControlPoints.size();
	int degreeU = pDegree;
	int degreeV = tDegree;
	std::vector<double> knotVectorU = pKnotVector;
	std::vector<double> knotVectorV = tKnotVector;
	std::vector<std::vector<XYZW>> controlPoints(pControlPoints.size(), std::vector<XYZW>(size));

	std::vector<double> vl(size);
	vl[0] = tKnotVector[tDegree];
	vl[size - 1] = tKnotVector[tKnotVector.size() - tDegree - 1];
	for (int k = 1; k < size - 1; k++)
	{
		vl[k] = 0.0;
		for (int i = k + 1; i < k + tDegree + 1; i++)
		{
			vl[k] += knotVectorV[i];
		}
		vl[k] /= tDegree;
	}

	std::vector<XYZ> tempB(vl.size());
	std::vector<XYZ> td = NurbsCurve::ComputeRationalCurveDerivatives(trajectory, 1, vl[0]);
	
	tempB[0] = td[1];
	if (MathUtils::IsAlmostEqualTo(td[1].GetY(),0.0))
	{
		tempB[0] = XYZ(0, 1, 0).CrossProduct(tempB[0]);
	}
	else
	{
		tempB[0] = XYZ(1, 0, 0).CrossProduct(tempB[0]);
	}
	tempB[0] = tempB[0].Normalize();

	for (int i = 1; i < vl.size(); i++)
	{
		td = NurbsCurve::ComputeRationalCurveDerivatives(trajectory, 1, vl[i]);
		XYZ ti = td[1].Normalize();
		XYZ bi = tempB[i - 1] - (tempB[i - 1] * ti) * ti;
		tempB[i] = bi.Normalize();
	}

	LN_Curve Bv;
	NurbsCurve::GlobalInterpolation(std::min(3,(int)(tempB.size()-1)), tempB, Bv, vl);

	for (int k = 0; k < size; k++) 
	{
		std::vector<XYZW> tempCp(pControlPoints.size());
		for (int i = 0; i < tempCp.size(); i++)
		{
			tempCp[i] = scale * pControlPoints[i];
		}		
		td = NurbsCurve::ComputeRationalCurveDerivatives(trajectory, 1, vl[k]);
		XYZ x = td[1].Normalize();
		XYZ z = NurbsCurve::GetPointOnCurve(Bv, vl[k]).Normalize();
		XYZ y = z.CrossProduct(x);

		Matrix4d R(y,x,z,XYZ(0,0,1));
		XYZ o = NurbsCurve::GetPointOnCurve(trajectory, vl[k]);
		Matrix4d tx = Matrix4d::CreateTranslation(o);
		Matrix4d A = tx.Multiply(R);
		for (int i = 0; i < tempCp.size(); i++)
		{
			controlPoints[i][k] = A.OfWeightedPoint(tempCp[i]);
		}
	}

	for (int i = 0; i < controlPoints.size(); i++) 
	{
		std::vector<XYZ> throughPoints(controlPoints[0].size());
		for (int k = 0; k < controlPoints[0].size(); ++k)
		{
			throughPoints[k] = controlPoints[i][k].ToXYZ(true);
		}
		LN_Curve R;
		NurbsCurve::GlobalInterpolation(degreeV, throughPoints, R, vl);
		for (int k = 0; k < controlPoints[0].size(); k++)
		{
			controlPoints[i][k] = R.ControlPoints[k];
		}	
	}

	surface.DegreeU = degreeU;
	surface.DegreeV = degreeV;
	surface.KnotVectorU = knotVectorU;
	surface.KnotVectorV = knotVectorV;
	surface.ControlPoints = controlPoints;
	return true;
}

void LNLib::NurbsSurface::CreateLoftSurface(const std::vector<LN_Curve>& sections, LN_Surface& surface)
{
	int degree_max = 0;
	for (int k = 0; k < sections.size(); k++)
	{
		LN_Curve current = sections[k];
		VALIDATE_ARGUMENT(current.Degree > 0, "degree", "Degree must greater than zero.");
		VALIDATE_ARGUMENT(current.KnotVector.size() > 0, "knotVector", "KnotVector size must greater than zero.");
		VALIDATE_ARGUMENT(ValidationUtils::IsValidKnotVector(current.KnotVector), "knotVector", "KnotVector must be a nondecreasing sequence of real numbers.");
		VALIDATE_ARGUMENT(current.ControlPoints.size() > 0, "controlPoints", "ControlPoints must contains one point at least.");
		VALIDATE_ARGUMENT(ValidationUtils::IsValidNurbs(current.Degree, current.KnotVector.size(), current.ControlPoints.size()), "controlPoints", "Arguments must fit: m = n + p + 1");

		if (degree_max < current.Degree)
		{
			degree_max = current.Degree;
		}
	}

	int size = sections.size();
	std::vector<LN_Curve> internals(size);
	std::vector<std::vector<XYZW>> curvesControlPoints(size);
	for (int k = 0; k < size; k++)
	{
		LN_Curve current = sections[k];
		if (degree_max > current.Degree)
		{
			int times = degree_max - current.Degree;
			LN_Curve tc;
			NurbsCurve::ElevateDegree(current, times, tc);
			current.Degree = degree_max;
			current.KnotVector = tc.KnotVector;
			current.ControlPoints = tc.ControlPoints;
		}
		curvesControlPoints.emplace_back(current.ControlPoints);
		internals.emplace_back(current);
	}

	int degreeU = degree_max;
	int degreeV = degree_max;
	
	std::vector<double> vl(size);
	vl[0] = 0;
	vl[size - 1] = 1;
	for (int k = 1; k <= size - 2; k++)
	{
		std::vector<XYZW> current = curvesControlPoints[k];

		int tsize = current.size();
		std::vector<XYZ> cps(tsize);
		for (int i = 0; i <= tsize; i++)
		{
			cps[i] = current[i].ToXYZ(true);
		}

		double average = 0.0;
		double length = Interpolation::GetTotalChordLength(cps);
		for (int i = 0; i <= tsize; i++)
		{
			double distance = curvesControlPoints[k][i].Distance(curvesControlPoints[k - 1][i]);
			average += distance / length;
		}
		vl[k] = vl[k - 1] + 1.0/(tsize+1) * average;
	}
	std::vector<double> knotVectorV = Interpolation::AverageKnotVector(degreeV, vl);

	std::vector<std::vector<XYZW>> controlPoints;
	int column = curvesControlPoints[0].size();
	for (int c = 0; c < column; c++)
	{
		std::vector<XYZ> temp(size);
		for (int k = 0; k < size; k++)
		{
			temp[k] = curvesControlPoints[k][c].ToXYZ(true);
		}
		LN_Curve tc;
		NurbsCurve::GlobalInterpolation(degreeV, temp, tc, vl);
		controlPoints.emplace_back(tc.ControlPoints);
	}

	surface.DegreeU = degreeU;
	surface.DegreeV = degreeV;
	surface.KnotVectorU = internals[0].KnotVector;
	surface.KnotVectorV = knotVectorV;
	surface.ControlPoints = controlPoints;
}

void LNLib::NurbsSurface::CreateSweepSurface(const LN_Curve& path, const std::vector<LN_Curve>& profiles, LN_Surface& surface)
{
	VALIDATE_ARGUMENT(path.Degree > 0, "degree", "Degree must greater than zero.");
	VALIDATE_ARGUMENT(path.KnotVector.size() > 0, "knotVector", "KnotVector size must greater than zero.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidKnotVector(path.KnotVector), "knotVector", "KnotVector must be a nondecreasing sequence of real numbers.");
	VALIDATE_ARGUMENT(path.ControlPoints.size() > 0, "controlPoints", "ControlPoints must contains one point at least.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidNurbs(path.Degree, path.KnotVector.size(), path.ControlPoints.size()), "controlPoints", "Arguments must fit: m = n + p + 1");

	for (int i = 0; i < profiles.size(); i++)
	{
		LN_Curve current = profiles[i];
		VALIDATE_ARGUMENT(current.Degree > 0, "degree", "Degree must greater than zero.");
		VALIDATE_ARGUMENT(current.KnotVector.size() > 0, "knotVector", "KnotVector size must greater than zero.");
		VALIDATE_ARGUMENT(ValidationUtils::IsValidKnotVector(current.KnotVector), "knotVector", "KnotVector must be a nondecreasing sequence of real numbers.");
		VALIDATE_ARGUMENT(current.ControlPoints.size() > 0, "controlPoints", "ControlPoints must contains one point at least.");
		VALIDATE_ARGUMENT(ValidationUtils::IsValidNurbs(current.Degree, current.KnotVector.size(), current.ControlPoints.size()), "controlPoints", "Arguments must fit: m = n + p + 1");
	}

	int profilesSize = profiles.size();
	double path_min = path.KnotVector[0];
	double path_max = path.KnotVector[path.KnotVector.size() - 1];
	double delta = (path_max - path_min) / profilesSize;

	std::vector<LN_Curve> sections(profilesSize);
	for (int i = 0; i < profilesSize; i++)
	{
		double param =  path_min + i * delta;
		XYZ point = NurbsCurve::GetPointOnCurve(path, param);
		std::vector<XYZ> ders = NurbsCurve::ComputeRationalCurveDerivatives(profiles[i], 1, param);
		XYZ tangent = ders[1];
		Matrix4d transform = Matrix4d::CreateTranslation(point);
		if (!tangent.IsZero())
		{
			XYZ binormal = NurbsCurve::Normal(profiles[i], CurveNormal::Binormal, param);
			double rad = binormal.AngleTo(tangent);
			Matrix4d rotation = Matrix4d::CreateRotationAtPoint(point, binormal, rad);
			transform = transform.Multiply(rotation);
		}
		LN_Curve newProfile;
		NurbsCurve::CreateTransformed(profiles[i], transform, newProfile);
		sections.emplace_back(newProfile);
	}

	CreateLoftSurface(sections, surface);
}

void LNLib::NurbsSurface::CreateGordonSurface(const std::vector<LN_Curve>& uCurves, const std::vector<LN_Curve>& vCurves, const std::vector<std::vector<XYZ>>& intersectionPoints, LN_Surface& surface)
{
	int degree_u_max = 0;
	for (int i = 0; i < uCurves.size(); i++)
	{
		LN_Curve current = uCurves[i];
		VALIDATE_ARGUMENT(current.Degree > 0, "degree", "Degree must greater than zero.");
		VALIDATE_ARGUMENT(current.KnotVector.size() > 0, "knotVector", "KnotVector size must greater than zero.");
		VALIDATE_ARGUMENT(ValidationUtils::IsValidKnotVector(current.KnotVector), "knotVector", "KnotVector must be a nondecreasing sequence of real numbers.");
		VALIDATE_ARGUMENT(current.ControlPoints.size() > 0, "controlPoints", "ControlPoints must contains one point at least.");
		VALIDATE_ARGUMENT(ValidationUtils::IsValidNurbs(current.Degree, current.KnotVector.size(), current.ControlPoints.size()), "controlPoints", "Arguments must fit: m = n + p + 1");

		if (current.Degree > degree_u_max)
		{
			degree_u_max = current.Degree;
		}
	}

	std::vector<LN_Curve> uInternals;
	for (int i = 0; i < uCurves.size(); i++)
	{
		LN_Curve current;
		NurbsCurve::Reparametrize(uCurves[i], 0, 1, current);

		if (degree_u_max > current.Degree)
		{
			int times = degree_u_max - current.Degree;
			LN_Curve tc;
			NurbsCurve::ElevateDegree(current, times, tc);
			uInternals.emplace_back(tc);
		}
		else
		{
			uInternals.emplace_back(current);
		}
	}

	std::vector<std::vector<double>> knotVectorsU;
	for (int i = 0; i < uInternals.size(); i++)
	{
		knotVectorsU.emplace_back(uInternals[i].KnotVector);
	}
	auto insertElements = KnotVectorUtils::GetInsertedKnotElements(knotVectorsU);
	for (int i = 0; i < uInternals.size(); i++)
	{
		auto insertElement = insertElements[i];
		if (insertElement.size() > 0)
		{
			LN_Curve tc;
			NurbsCurve::RefineKnotVector(uInternals[i], insertElement, tc);
			uInternals[i] = tc;
		}
	}

	int degree_v_max = 0;
	for (int i = 0; i < vCurves.size(); i++)
	{
		LN_Curve current = vCurves[i];
		VALIDATE_ARGUMENT(current.Degree > 0, "degree", "Degree must greater than zero.");
		VALIDATE_ARGUMENT(current.KnotVector.size() > 0, "knotVector", "KnotVector size must greater than zero.");
		VALIDATE_ARGUMENT(ValidationUtils::IsValidKnotVector(current.KnotVector), "knotVector", "KnotVector must be a nondecreasing sequence of real numbers.");
		VALIDATE_ARGUMENT(current.ControlPoints.size() > 0, "controlPoints", "ControlPoints must contains one point at least.");
		VALIDATE_ARGUMENT(ValidationUtils::IsValidNurbs(current.Degree, current.KnotVector.size(), current.ControlPoints.size()), "controlPoints", "Arguments must fit: m = n + p + 1");
	
		if (current.Degree > degree_v_max)
		{
			degree_v_max = current.Degree;
		}
	}

	std::vector<LN_Curve> vInternals;
	for (int i = 0; i < vCurves.size(); i++)
	{
		LN_Curve current;
		NurbsCurve::Reparametrize(vCurves[i], 0, 1, current);
		if (degree_v_max > current.Degree)
		{
			int times = degree_v_max - current.Degree;
			LN_Curve tc;
			NurbsCurve::ElevateDegree(current, times, tc);
			vInternals.emplace_back(tc);
		}
		else
		{
			vInternals.emplace_back(current);
		}
	}

	std::vector<std::vector<double>> knotVectorsV;
	for (int i = 0; i < vInternals.size(); i++)
	{
		knotVectorsV.emplace_back(vInternals[i].KnotVector);
	}
	insertElements = KnotVectorUtils::GetInsertedKnotElements(knotVectorsV);
	for (int i = 0; i < vInternals.size(); i++)
	{
		auto insertElement = insertElements[i];
		if (insertElement.size() > 0)
		{
			LN_Curve tc;
			NurbsCurve::RefineKnotVector(vInternals[i], insertElement, tc);
			vInternals[i] = tc;
		}
	}

	int rows = intersectionPoints.size();
	int columns = intersectionPoints[0].size();

	LN_Surface loftSurfaceV;
	CreateLoftSurface(uInternals, loftSurfaceV);
	LN_Surface ts;
	CreateLoftSurface(vInternals, ts);
	std::vector<std::vector<XYZW>> transposedControlPoints;
	MathUtils::Transpose(ts.ControlPoints, transposedControlPoints);
	LN_Surface loftSurfaceU;
	loftSurfaceU.DegreeU = ts.DegreeV;
	loftSurfaceU.DegreeV = ts.DegreeU;
	loftSurfaceU.KnotVectorU = ts.KnotVectorV;
	loftSurfaceU.KnotVectorV = ts.KnotVectorU;
	loftSurfaceU.ControlPoints = transposedControlPoints;

	int degreeU = std::min(columns - 1, degree_u_max);
	int degreeV = std::min(rows - 1, degree_v_max);
	GlobalInterpolation(intersectionPoints, degreeU, degreeV, ts);
	MathUtils::Transpose(ts.ControlPoints, transposedControlPoints);
	LN_Surface interpolatedSurface;
	Swap(ts, interpolatedSurface);

	{
		int lsu_degreeU = loftSurfaceU.DegreeU;
		int lsv_degreeU = loftSurfaceV.DegreeU;
		int inp_degreeU = interpolatedSurface.DegreeU;

		int degreeU = std::max(lsu_degreeU, std::max(lsv_degreeU, inp_degreeU));
		if (degreeU > lsu_degreeU)
		{
			int times = degreeU - lsu_degreeU;
			LN_Surface temp;
			ElevateDegree(loftSurfaceU, times, true, temp);
			loftSurfaceU = temp;
		}

		if (degreeU > lsv_degreeU)
		{
			int times = degreeU - lsv_degreeU;
			LN_Surface temp;
			ElevateDegree(loftSurfaceV, times, true, temp);
			loftSurfaceV = temp;
		}

		if (degreeU > inp_degreeU)
		{
			int times = degreeU - inp_degreeU;
			LN_Surface temp;
			ElevateDegree(interpolatedSurface, times, true, temp);
			interpolatedSurface = temp;
		}

		int lsu_degreeV = loftSurfaceU.DegreeV;
		int lsv_degreeV = loftSurfaceV.DegreeV;
		int inp_degreeV = interpolatedSurface.DegreeV;

		int degreeV = std::max(lsu_degreeV, std::max(lsv_degreeV, inp_degreeV));
		if (degreeV > lsu_degreeV)
		{
			int times = degreeV - lsu_degreeV;
			LN_Surface temp;
			ElevateDegree(loftSurfaceU, times, false, temp);
			loftSurfaceU = temp;
		}

		if (degreeV > lsv_degreeV)
		{
			int times = degreeV - lsv_degreeV;
			LN_Surface temp;
			ElevateDegree(loftSurfaceV, times, false, temp);
			loftSurfaceV = temp;
		}

		if (degreeV > inp_degreeV)
		{
			int times = degreeV - inp_degreeV;
			LN_Surface temp;
			ElevateDegree(interpolatedSurface, times, false, temp);
			interpolatedSurface = temp;
		}
	}

	{
		std::vector<std::vector<double>> knotVectorsU;
		knotVectorsU.emplace_back(loftSurfaceU.KnotVectorU);
		knotVectorsU.emplace_back(loftSurfaceV.KnotVectorU);
		knotVectorsU.emplace_back(interpolatedSurface.KnotVectorU);

		auto insertElements = KnotVectorUtils::GetInsertedKnotElements(knotVectorsU);
		if (insertElements[0].size() > 0)
		{
			LN_Surface temp;
			RefineKnotVector(loftSurfaceU, insertElements[0], true, temp);
			loftSurfaceU = temp;
		}
		if (insertElements[1].size() > 0)
		{
			LN_Surface temp;
			RefineKnotVector(loftSurfaceV, insertElements[1], true, temp);
			loftSurfaceV = temp;
		}
		if (insertElements[2].size() > 0)
		{
			LN_Surface temp;
			RefineKnotVector(interpolatedSurface, insertElements[2], true, temp);
			interpolatedSurface = temp;
		}
	}

	{
		std::vector<std::vector<double>> knotVectorsV;
		knotVectorsV.emplace_back(loftSurfaceU.KnotVectorV);
		knotVectorsV.emplace_back(loftSurfaceV.KnotVectorV);
		knotVectorsV.emplace_back(interpolatedSurface.KnotVectorV);

		auto insertElements = KnotVectorUtils::GetInsertedKnotElements(knotVectorsV);
		if (insertElements[0].size() > 0)
		{
			LN_Surface temp;
			RefineKnotVector(loftSurfaceU, insertElements[0], false, temp);
			loftSurfaceU = temp;
		}
		if (insertElements[1].size() > 0)
		{
			LN_Surface temp;
			RefineKnotVector(loftSurfaceV, insertElements[1], false, temp);
			loftSurfaceV = temp;
		}
		if (insertElements[2].size() > 0)
		{
			LN_Surface temp;
			RefineKnotVector(interpolatedSurface, insertElements[2], false, temp);
			interpolatedSurface = temp;
		}
	}

	surface.DegreeU = interpolatedSurface.DegreeU;
	surface.DegreeV = interpolatedSurface.DegreeV;
	surface.KnotVectorU = interpolatedSurface.KnotVectorU;
	surface.KnotVectorV = interpolatedSurface.KnotVectorV;

	int row = interpolatedSurface.ControlPoints.size();
	int column = interpolatedSurface.ControlPoints[0].size();

	std::vector<std::vector<XYZW>> controlPoints(row, std::vector<XYZW>(column));
	for (int i = 0; i < row; i++)
	{
		for (int j = 0; j < column; j++)
		{
			controlPoints[i][j] = loftSurfaceU.ControlPoints[i][j] + loftSurfaceV.ControlPoints[i][j] - interpolatedSurface.ControlPoints[i][j];
		}
	}
	surface.ControlPoints = controlPoints;
}

void LNLib::NurbsSurface::CreateCoonsSurface(const LN_Curve& curve0, const LN_Curve& curve1, const LN_Curve& curve2, const LN_Curve& curve3, LN_Surface& surface)
{
	std::vector<LN_Curve> nurbs(4);
	nurbs.emplace_back(curve0);
	nurbs.emplace_back(curve1);
	nurbs.emplace_back(curve2);
	nurbs.emplace_back(curve3);

	for (int i = 0; i < nurbs.size(); i++)
	{
		LN_Curve current = nurbs[i];
		VALIDATE_ARGUMENT(current.Degree > 0, "degree", "Degree must greater than zero.");
		VALIDATE_ARGUMENT(current.KnotVector.size() > 0, "knotVector", "KnotVector size must greater than zero.");
		VALIDATE_ARGUMENT(ValidationUtils::IsValidKnotVector(current.KnotVector), "knotVector", "KnotVector must be a nondecreasing sequence of real numbers.");
		VALIDATE_ARGUMENT(current.ControlPoints.size() > 0, "controlPoints", "ControlPoints must contains one point at least.");
		VALIDATE_ARGUMENT(ValidationUtils::IsValidNurbs(current.Degree, current.KnotVector.size(), current.ControlPoints.size()), "controlPoints", "Arguments must fit: m = n + p + 1");

		LN_Curve tc;
		NurbsCurve::Reparametrize(current, 0, 1, tc);
		nurbs[i] = tc;
	}

	auto n0 = nurbs[0]; auto n2 = nurbs[2];
	int degree_n0 = n0.Degree;
	int degree_n2 = n2.Degree;
	if (degree_n0 > degree_n2)
	{
		int times = degree_n0 - degree_n2;
		LN_Curve tc;
		NurbsCurve::ElevateDegree(n2, times, tc);
		n2 = tc;
	}
	if (degree_n2 > degree_n0)
	{
		int times = degree_n2 - degree_n0;
		LN_Curve tc;
		NurbsCurve::ElevateDegree(n0, times, tc);
		n0 = tc;
	}
	if (n0.KnotVector != n2.KnotVector)
	{
		std::vector<double> insertedKnotElement0;
		std::vector<double> insertedKnotElement2;
		KnotVectorUtils::GetInsertedKnotElement(n0.KnotVector, n2.KnotVector, insertedKnotElement0, insertedKnotElement2);

		if (insertedKnotElement0.size() > 0)
		{
			LN_Curve tc;
			NurbsCurve::RefineKnotVector(n0, insertedKnotElement0, tc);
			n0 = tc;
		}
		if (insertedKnotElement2.size() > 0)
		{
			LN_Curve tc;
			NurbsCurve::RefineKnotVector(n2, insertedKnotElement2, tc);
			n2 = tc;
		}
	}

	auto n1 = nurbs[1]; auto n3 = nurbs[3];
	int degree_n1 = n1.Degree;
	int degree_n3 = n3.Degree;
	if (degree_n1 > degree_n3)
	{
		int times = degree_n1 - degree_n3;
		LN_Curve tc;
		NurbsCurve::ElevateDegree(n3, times, tc);
		n3 = tc;
	}
	if (degree_n3 > degree_n1)
	{
		int times = degree_n3 - degree_n1;
		LN_Curve tc;
		NurbsCurve::ElevateDegree(n1, times, tc);
		n1 = tc;
	}
	if (n1.KnotVector != n3.KnotVector)
	{
		std::vector<double> insertedKnotElement1;
		std::vector<double> insertedKnotElement3;
		KnotVectorUtils::GetInsertedKnotElement(n1.KnotVector, n3.KnotVector, insertedKnotElement1, insertedKnotElement3);

		if (insertedKnotElement1.size() > 0)
		{
			LN_Curve tc;
			NurbsCurve::RefineKnotVector(n1, insertedKnotElement1, tc);
			n1 = tc;
		}
		if (insertedKnotElement3.size() > 0)
		{
			LN_Curve tc;
			NurbsCurve::RefineKnotVector(n3, insertedKnotElement3, tc);
			n3 = tc;
		}
	}

	{
		LN_Curve tc;
		NurbsCurve::Reverse(n0, tc);
		n0 = tc;
	}

	{
		LN_Curve tc;
		NurbsCurve::Reverse(n3, tc);
		n3 = tc;
	}

	LN_Surface ruledSurface0;
	{
		std::vector<std::vector<XYZW>> cps;
		CreateRuledSurface(n0, n2, ruledSurface0);
	}
	
	LN_Surface ruledSurface1;
	{
		LN_Surface ts;
		CreateRuledSurface(n1, n3, ts);
		Swap(ts, ruledSurface1);
	}

	LN_Surface bilinearSurface;
	{
		XYZ point1 = NurbsCurve::GetPointOnCurve(n0, n0.KnotVector[0]);
		XYZ point2 = NurbsCurve::GetPointOnCurve(n2, n2.KnotVector[0]);
		XYZ point3 = NurbsCurve::GetPointOnCurve(n0, n0.KnotVector[n0.KnotVector.size() - 1]);
		XYZ point4 = NurbsCurve::GetPointOnCurve(n2, n2.KnotVector[n2.KnotVector.size() - 1]);
		int degree_u;
		int degree_v;
		std::vector<double> kv_u;
		std::vector<double> kv_v;
		std::vector<std::vector<XYZW>> cps;
		CreateBilinearSurface(point1, point2, point3, point4, bilinearSurface);
	}

	{
		int rs0_degreeU = ruledSurface0.DegreeU;
		int rs1_degreeU = ruledSurface1.DegreeU;
		int bs_degreeU = bilinearSurface.DegreeU;

		int degreeU = std::max(rs0_degreeU, std::max(rs1_degreeU, bs_degreeU));
		if (degreeU > rs0_degreeU)
		{
			int times = degreeU - rs0_degreeU;
			LN_Surface temp;
			ElevateDegree(ruledSurface0, times, true, temp);
			ruledSurface0 = temp;
		}

		if (degreeU > rs1_degreeU)
		{
			int times = degreeU - rs1_degreeU;
			LN_Surface temp;
			ElevateDegree(ruledSurface1, times, true, temp);
			ruledSurface1 = temp;
		}

		if (degreeU > bs_degreeU)
		{
			int times = degreeU - bs_degreeU;
			LN_Surface temp;
			ElevateDegree(bilinearSurface, times, true, temp);
			bilinearSurface = temp;
		}

		int rs0_degreeV = ruledSurface0.DegreeV;
		int rs1_degreeV = ruledSurface1.DegreeV;
		int bs_degreeV = bilinearSurface.DegreeV;

		int degreeV = std::max(rs0_degreeV, std::max(rs1_degreeV, bs_degreeV));
		if (degreeV > rs0_degreeV)
		{
			int times = degreeV - rs0_degreeV;
			LN_Surface temp;
			ElevateDegree(ruledSurface0, times, false, temp);
			ruledSurface0 = temp;
		}

		if (degreeV > rs1_degreeV)
		{
			int times = degreeV - rs1_degreeV;
			LN_Surface temp;
			ElevateDegree(ruledSurface1, times, false, temp);
			ruledSurface1 = temp;
		}

		if (degreeV > bs_degreeV)
		{
			int times = degreeV - bs_degreeV;
			LN_Surface temp;
			ElevateDegree(bilinearSurface, times, false, temp);
			bilinearSurface = temp;
		}
	}
	
	{
		std::vector<std::vector<double>> knotVectorsU;
		knotVectorsU.emplace_back(ruledSurface0.KnotVectorU);
		knotVectorsU.emplace_back(ruledSurface1.KnotVectorU);
		knotVectorsU.emplace_back(bilinearSurface.KnotVectorU);

		auto insertElements = KnotVectorUtils::GetInsertedKnotElements(knotVectorsU);
		if (insertElements[0].size() > 0)
		{
			LN_Surface temp;
			RefineKnotVector(ruledSurface0, insertElements[0], true, temp);
			ruledSurface0 = temp;
		}
		if (insertElements[1].size() > 0)
		{
			LN_Surface temp;
			RefineKnotVector(ruledSurface1, insertElements[1], true, temp);
			ruledSurface1 = temp;
		}
		if (insertElements[2].size() > 0)
		{
			LN_Surface temp;
			RefineKnotVector(bilinearSurface, insertElements[2], true, temp);
			bilinearSurface = temp;
		}
	}

	{
		std::vector<std::vector<double>> knotVectorsV;
		knotVectorsV.emplace_back(ruledSurface0.KnotVectorV);
		knotVectorsV.emplace_back(ruledSurface1.KnotVectorV);
		knotVectorsV.emplace_back(bilinearSurface.KnotVectorV);

		auto insertElements = KnotVectorUtils::GetInsertedKnotElements(knotVectorsV);
		if (insertElements[0].size() > 0)
		{
			LN_Surface temp;
			RefineKnotVector(ruledSurface0, insertElements[0], false, temp);
			ruledSurface0 = temp;
		}
		if (insertElements[1].size() > 0)
		{
			LN_Surface temp;
			RefineKnotVector(ruledSurface1, insertElements[1], false, temp);
			ruledSurface1 = temp;
		}
		if (insertElements[2].size() > 0)
		{
			LN_Surface temp;
			RefineKnotVector(bilinearSurface, insertElements[2], false, temp);
			bilinearSurface = temp;
		}
	}

	surface.DegreeU = ruledSurface0.DegreeU;
	surface.DegreeV = ruledSurface0.DegreeV;
	surface.KnotVectorU = ruledSurface0.KnotVectorU;
	surface.KnotVectorV = ruledSurface0.KnotVectorV;

	int row = ruledSurface0.ControlPoints.size();
	int column = ruledSurface0.ControlPoints[0].size();

	std::vector<std::vector<XYZW>> controlPoints(row, std::vector<XYZW>(column));
	for (int i = 0; i < row; i++)
	{
		for (int j = 0; j < column; j++)
		{
			controlPoints[i][j] = ruledSurface0.ControlPoints[i][j] + ruledSurface1.ControlPoints[i][j] - bilinearSurface.ControlPoints[i][j];
		}
	}
	surface.ControlPoints = controlPoints;
}


