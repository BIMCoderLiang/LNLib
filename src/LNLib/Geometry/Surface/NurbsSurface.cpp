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
#include "MathUtils.h"
#include "NurbsCurve.h"
#include "BsplineSurface.h"
#include "Projection.h"
#include "Intersection.h"
#include "Interpolation.h"
#include "ValidationUtils.h"
#include "LNLibExceptions.h"
#include <algorithm>

namespace LNLib
{
	std::vector<std::vector<XYZ>> ToXYZ(const std::vector<std::vector<XYZW>>& surfacePoints)
	{
		int row = static_cast<int>(surfacePoints.size());
		int column = static_cast<int>(surfacePoints[0].size());

		std::vector<std::vector<XYZ>> result;
		result.resize(row);
		for (int i = 0; i < row; i++)
		{
			result[i].resize(column);
			for (int j = 0; j < column; j++)
			{
				result[i][j] = const_cast<XYZW&>(surfacePoints[i][j]).ToXYZ(true);
			}
		}
		return result;
	}

	std::vector<std::vector<XYZW>> ToXYZW(const std::vector<std::vector<XYZ>>& surfacePoints)
	{
		int row = static_cast<int>(surfacePoints.size());
		int column = static_cast<int>(surfacePoints[0].size());

		std::vector<std::vector<XYZW>> result;
		result.resize(row);
		for (int i = 0; i < row; i++)
		{
			result[i].resize(column);
			for (int j = 0; j < column; j++)
			{
				result[i][j] = XYZW(const_cast<XYZ&>(surfacePoints[i][j]),1);
			}
		}
		return result;
	}
}

LNLib::XYZ LNLib::NurbsSurface::GetPointOnSurface(int degreeU, int degreeV, const std::vector<double>& knotVectorU, const std::vector<double>& knotVectorV, UV uv, const std::vector<std::vector<XYZW>>& controlPoints)
{
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


std::vector<std::vector<LNLib::XYZ>> LNLib::NurbsSurface::ComputeRationalSurfaceDerivatives(int degreeU, int degreeV, int derivative, const std::vector<double>& knotVectorU, const std::vector<double>& knotVectorV, UV uv, const std::vector<std::vector<XYZW>>& controlPoints)
{
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

void LNLib::NurbsSurface::InsertKnot(int degree, const std::vector<double>& knotVector, const std::vector<std::vector<XYZW>>& controlPoints, double insertKnot, int times, bool isUDirection, std::vector<double>& insertedKnotVector, std::vector<std::vector<XYZW>>& updatedControlPoints)
{
	VALIDATE_ARGUMENT(degree > 0, "degree", "Degree must greater than zero.");
	VALIDATE_ARGUMENT(knotVector.size() > 0, "knotVector", "KnotVector size must greater than zero.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidKnotVector(knotVector), "knotVector", "KnotVector must be a nondecreasing sequence of real numbers.");
	VALIDATE_ARGUMENT(controlPoints.size() > 0, "controlPoints", "ControlPoints must contains one point at least.");
	VALIDATE_ARGUMENT(times > 0, "times", "Times must greater than zero.");
	if (isUDirection)
	{
		VALIDATE_ARGUMENT(ValidationUtils::IsValidNurbs(degree, knotVector.size(), controlPoints.size()), "controlPoints", "Arguments must fit: m = n + p + 1");
	}
	else
	{
		VALIDATE_ARGUMENT(ValidationUtils::IsValidNurbs(degree, knotVector.size(), controlPoints[0].size()), "controlPoints", "Arguments must fit: m = n + p + 1");
	}
	
	int knotSpanIndex = Polynomials::GetKnotSpanIndex(degree, knotVector, insertKnot);
	int multiplicity = Polynomials::GetKnotMultiplicity(knotVector, insertKnot);

	if (multiplicity == degree)
	{
		insertedKnotVector = knotVector;
		updatedControlPoints = controlPoints;
		return;
	}

	if ((times + multiplicity) > degree)
	{
		times = degree - multiplicity;
	}

	insertedKnotVector.resize(knotVector.size() + times);
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

	int rows = static_cast<int>(controlPoints.size());
	int columns = static_cast<int>(controlPoints[0].size());

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
	}
}

void LNLib::NurbsSurface::RefineKnotVector(int degreeU, int degreeV, const std::vector<double>& knotVectorU, const std::vector<double>& knotVectorV, const std::vector<std::vector<XYZW>>& controlPoints, std::vector<double>& insertKnotElements, bool isUDirection, std::vector<double>& insertedKnotVectorU, std::vector<double>& insertedKnotVectorV, std::vector<std::vector<XYZW>>& updatedControlPoints)
{
	VALIDATE_ARGUMENT(degreeU > 0, "degreeU", "Degree must greater than zero.");
	VALIDATE_ARGUMENT(degreeV > 0, "degreeV", "Degree must greater than zero.");
	VALIDATE_ARGUMENT(knotVectorU.size() > 0, "knotVectorU", "KnotVector size must greater than zero.");
	VALIDATE_ARGUMENT(knotVectorV.size() > 0, "knotVectorV", "KnotVector size must greater than zero.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidKnotVector(knotVectorU), "knotVectorU", "KnotVector must be a nondecreasing sequence of real numbers.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidKnotVector(knotVectorV), "knotVectorV", "KnotVector must be a nondecreasing sequence of real numbers.");
	VALIDATE_ARGUMENT(controlPoints.size() > 0, "controlPoints", "ControlPoints must contains one point at least.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidNurbs(degreeU, knotVectorU.size(), controlPoints.size()), "controlPoints", "Arguments must fit: m = n + p + 1");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidNurbs(degreeV, knotVectorV.size(), controlPoints[0].size()), "controlPoints", "Arguments must fit: m = n + p + 1");
	VALIDATE_ARGUMENT(insertKnotElements.size() > 0, "insertKnotElements", "insertKnotElements size must greater than zero.");

	std::vector<double> tempKnotVector;
	std::vector<std::vector<XYZW>> tempControlPoints;

	if (isUDirection)
	{
		std::vector<std::vector<XYZW>> transposedControlPoints;
		MathUtils::Transpose(controlPoints, transposedControlPoints);
		std::vector<XYZW> temp;
		for (int i = 0; i < transposedControlPoints.size(); i++)
		{
			NurbsCurve::RefineKnotVector(degreeU, knotVectorU, transposedControlPoints[i], insertKnotElements, tempKnotVector, temp);
			tempControlPoints.emplace_back(temp);
		}
		insertedKnotVectorU = tempKnotVector;
		insertedKnotVectorV = knotVectorV;
		MathUtils::Transpose(tempControlPoints, updatedControlPoints);
	}
	else
	{
		std::vector<XYZW> temp;
		for (int i = 0; i < controlPoints.size(); i++)
		{
			NurbsCurve::RefineKnotVector(degreeU, knotVectorU, controlPoints[i], insertKnotElements, tempKnotVector, temp);
			tempControlPoints.emplace_back(temp);
		}
		insertedKnotVectorU = knotVectorU;
		insertedKnotVectorV = tempKnotVector;
		updatedControlPoints = tempControlPoints;
	}
}

std::vector<std::vector<std::vector<LNLib::XYZW>>> LNLib::NurbsSurface::DecomposeToBeziers(int degreeU, int degreeV, const std::vector<double>& knotVectorU, const std::vector<double>& knotVectorV, const std::vector<std::vector<XYZW>>& controlPoints)
{
	VALIDATE_ARGUMENT(degreeU > 0, "degreeU", "Degree must greater than zero.");
	VALIDATE_ARGUMENT(degreeV > 0, "degreeV", "Degree must greater than zero.");
	VALIDATE_ARGUMENT(knotVectorU.size() > 0, "knotVectorU", "KnotVector size must greater than zero.");
	VALIDATE_ARGUMENT(knotVectorV.size() > 0, "knotVectorV", "KnotVector size must greater than zero.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidKnotVector(knotVectorU), "knotVectorU", "KnotVector must be a nondecreasing sequence of real numbers.");
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

void LNLib::NurbsSurface::RemoveKnot(int degreeU, int degreeV, const std::vector<double>& knotVectorU, const std::vector<double>& knotVectorV, const std::vector<std::vector<XYZW>>& controlPoints, double removeKnot, int times, bool isUDirection, std::vector<double>& restKnotVectorU, std::vector<double>& restKnotVectorV, std::vector<std::vector<XYZW>>& updatedControlPoints)
{
	VALIDATE_ARGUMENT(degreeU > 0, "degreeU", "Degree must greater than zero.");
	VALIDATE_ARGUMENT(degreeV > 0, "degreeV", "Degree must greater than zero.");
	VALIDATE_ARGUMENT(knotVectorU.size() > 0, "knotVectorU", "KnotVector size must greater than zero.");
	VALIDATE_ARGUMENT(knotVectorV.size() > 0, "knotVectorV", "KnotVector size must greater than zero.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidKnotVector(knotVectorU), "knotVectorU", "KnotVector must be a nondecreasing sequence of real numbers.");
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
	if (isUDirection)
	{
		std::vector<std::vector<XYZW>> transposedControlPoints;
		MathUtils::Transpose(controlPoints, transposedControlPoints);
		std::vector<XYZW> temp;
		for (int i = 0; i < transposedControlPoints.size(); i++)
		{
			NurbsCurve::RemoveKnot(degreeU, knotVectorU, transposedControlPoints[i], removeKnot, times, tempKnotVector, temp);
			tempControlPoints.emplace_back(temp);
		}
		restKnotVectorU = tempKnotVector;
		restKnotVectorV = knotVectorV;
		MathUtils::Transpose(tempControlPoints, updatedControlPoints);
	}
	else
	{
		std::vector<XYZW> temp;
		for (int i = 0; i < controlPoints.size(); i++)
		{
			NurbsCurve::RemoveKnot(degreeU, knotVectorU, controlPoints[i], removeKnot, times, tempKnotVector, temp);
			tempControlPoints.emplace_back(temp);
		}
		restKnotVectorU = knotVectorU;
		restKnotVectorV = tempKnotVector;
		updatedControlPoints = tempControlPoints;
	}
}

void LNLib::NurbsSurface::ElevateDegree(int degreeU, int degreeV, const std::vector<double>& knotVectorU, const std::vector<double>& knotVectorV, const std::vector<std::vector<XYZW>>& controlPoints, int times, bool isUDirection, std::vector<double>& updatedKnotVectorU, std::vector<double>& updatedKnotVectorV, std::vector<std::vector<XYZW>>& updatedControlPoints)
{
	VALIDATE_ARGUMENT(degreeU > 0, "degreeU", "Degree must greater than zero.");
	VALIDATE_ARGUMENT(degreeV > 0, "degreeV", "Degree must greater than zero.");
	VALIDATE_ARGUMENT(knotVectorU.size() > 0, "knotVectorU", "KnotVector size must greater than zero.");
	VALIDATE_ARGUMENT(knotVectorV.size() > 0, "knotVectorV", "KnotVector size must greater than zero.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidKnotVector(knotVectorU), "knotVectorU", "KnotVector must be a nondecreasing sequence of real numbers.");
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
		std::vector<XYZW> temp;
		for (int i = 0; i < transposedControlPoints.size(); i++)
		{
			NurbsCurve::ElevateDegree(degreeU, knotVectorU, transposedControlPoints[i], times, tempKnotVector, temp);
			tempControlPoints.emplace_back(temp);
		}
		updatedKnotVectorU = tempKnotVector;
		updatedKnotVectorV = knotVectorV;
		MathUtils::Transpose(tempControlPoints, updatedControlPoints);
	}
	else
	{
		std::vector<XYZW> temp;
		for (int i = 0; i < controlPoints.size(); i++)
		{
			NurbsCurve::ElevateDegree(degreeU, knotVectorU, controlPoints[i], times, tempKnotVector, temp);
			tempControlPoints.emplace_back(temp);
		}
		updatedKnotVectorU = knotVectorU;
		updatedKnotVectorV = tempKnotVector;
		updatedControlPoints = tempControlPoints;
	}
}

bool LNLib::NurbsSurface::ReduceDegree(int degreeU, int degreeV, const std::vector<double>& knotVectorU, const std::vector<double>& knotVectorV, const std::vector<std::vector<XYZW>>& controlPoints, bool isUDirection, std::vector<double>& updatedKnotVectorU, std::vector<double>& updatedKnotVectorV, std::vector<std::vector<XYZW>>& updatedControlPoints)
{
	VALIDATE_ARGUMENT(ValidationUtils::IsValidDegreeReduction(degreeU), "degreeU", "Degree must greater than one.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidDegreeReduction(degreeV), "degreeV", "Degree must greater than one.");
	VALIDATE_ARGUMENT(knotVectorU.size() > 0, "knotVectorU", "KnotVector size must greater than zero.");
	VALIDATE_ARGUMENT(knotVectorV.size() > 0, "knotVectorV", "KnotVector size must greater than zero.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidKnotVector(knotVectorU), "knotVectorU", "KnotVector must be a nondecreasing sequence of real numbers.");
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
		std::vector<XYZW> temp;
		for (int i = 0; i < transposedControlPoints.size(); i++)
		{
			bool result = NurbsCurve::ReduceDegree(degreeU, knotVectorU, transposedControlPoints[i], tempKnotVector, temp);
			if (!result)
				return false;
			tempControlPoints.emplace_back(temp);
		}
		updatedKnotVectorU = tempKnotVector;
		updatedKnotVectorV = knotVectorV;
		MathUtils::Transpose(tempControlPoints, updatedControlPoints);
	}
	else
	{
		std::vector<XYZW> temp;
		for (int i = 0; i < controlPoints.size(); i++)
		{
			bool result = NurbsCurve::ReduceDegree(degreeU, knotVectorU, controlPoints[i],  tempKnotVector, temp);
			if (!result)
				return false;
			tempControlPoints.emplace_back(temp);
		}
		updatedKnotVectorU = knotVectorU;
		updatedKnotVectorV = tempKnotVector;
		updatedControlPoints = tempControlPoints;
	}
	return true;
}

void LNLib::NurbsSurface::EquallyTessellate(int degreeU, int degreeV, const std::vector<double>& knotVectorU, const std::vector<double>& knotVectorV, const std::vector<std::vector<XYZW>>& controlPoints, std::vector<XYZ>& tessellatedPoints, std::vector<UV>& correspondingKnots)
{
	VALIDATE_ARGUMENT(degreeU > 0, "degreeU", "Degree must greater than zero.");
	VALIDATE_ARGUMENT(degreeV > 0, "degreeV", "Degree must greater than zero.");
	VALIDATE_ARGUMENT(knotVectorU.size() > 0, "knotVectorU", "KnotVector size must greater than zero.");
	VALIDATE_ARGUMENT(knotVectorV.size() > 0, "knotVectorV", "KnotVector size must greater than zero.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidKnotVector(knotVectorU), "knotVectorU", "KnotVector must be a nondecreasing sequence of real numbers.");
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
			tessellatedPoints.emplace_back(GetPointOnSurface(degreeU, degreeV, knotVectorU, knotVectorV, uv, controlPoints));
		}
	}

	correspondingKnots.emplace_back(UV(knotVectorU[knotVectorU.size() - 1], knotVectorV[knotVectorV.size() - 1]));
	tessellatedPoints.emplace_back(const_cast<XYZW&>(controlPoints[controlPoints.size() - 1][controlPoints[0].size()-1]).ToXYZ(true));
}

LNLib::UV LNLib::NurbsSurface::GetParamOnSurface(int degreeU, int degreeV, const std::vector<double>& knotVectorU, const std::vector<double>& knotVectorV, const std::vector<std::vector<XYZW>>& controlPoints, const XYZ& givenPoint)
{
	VALIDATE_ARGUMENT(degreeU > 0, "degreeU", "Degree must greater than zero.");
	VALIDATE_ARGUMENT(degreeV > 0, "degreeV", "Degree must greater than zero.");
	VALIDATE_ARGUMENT(knotVectorU.size() > 0, "knotVectorU", "KnotVector size must greater than zero.");
	VALIDATE_ARGUMENT(knotVectorV.size() > 0, "knotVectorV", "KnotVector size must greater than zero.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidKnotVector(knotVectorU), "knotVectorU", "KnotVector must be a nondecreasing sequence of real numbers.");
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
	EquallyTessellate(degreeU, degreeV, knotVectorU, knotVectorV, controlPoints, tessellatedPoints, correspondingKnots);
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
		std::vector<std::vector<XYZ>> derivatives = ComputeRationalSurfaceDerivatives(degreeU, degreeV, 2, knotVectorU, knotVectorV, param, controlPoints);
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

bool LNLib::NurbsSurface::GetUVTangent(int degreeU, int degreeV, const std::vector<double>& knotVectorU, const std::vector<double>& knotVectorV, const std::vector<std::vector<XYZW>>& controlPoints, const UV param, const XYZ& tangent, UV& uvTangent)
{
	VALIDATE_ARGUMENT(degreeU > 0, "degreeU", "Degree must greater than zero.");
	VALIDATE_ARGUMENT(degreeV > 0, "degreeV", "Degree must greater than zero.");
	VALIDATE_ARGUMENT(knotVectorU.size() > 0, "knotVectorU", "KnotVector size must greater than zero.");
	VALIDATE_ARGUMENT(knotVectorV.size() > 0, "knotVectorV", "KnotVector size must greater than zero.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidKnotVector(knotVectorU), "knotVectorU", "KnotVector must be a nondecreasing sequence of real numbers.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidKnotVector(knotVectorV), "knotVectorV", "KnotVector must be a nondecreasing sequence of real numbers.");
	VALIDATE_ARGUMENT(controlPoints.size() > 0, "controlPoints", "ControlPoints must contains one point at least.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidNurbs(degreeU, knotVectorU.size(), controlPoints.size()), "controlPoints", "Arguments must fit: m = n + p + 1");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidNurbs(degreeV, knotVectorV.size(), controlPoints[0].size()), "controlPoints", "Arguments must fit: m = n + p + 1");

	std::vector<std::vector<XYZ>> derivatives = ComputeRationalSurfaceDerivatives(degreeU, degreeV, 1, knotVectorU, knotVectorV, param, controlPoints);
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

void LNLib::NurbsSurface::CreateBilinearSurface(const XYZ& point1, const XYZ& point2, const XYZ& point3, const XYZ& point4, int& degreeU, int& degreeV, std::vector<double>& knotVectorU, std::vector<double>& knotVectorV, std::vector<std::vector<XYZW>>& controlPoints)
{
	int degree = 3;
	degreeU = degreeV = degree;

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
}

bool LNLib::NurbsSurface::CreateCylindricalSurface(const XYZ& origin, const XYZ& xAxis, const XYZ& yAxis, double startRad, double endRad, double radius, double height, int& degreeU, int& degreeV, std::vector<double>& knotVectorU, std::vector<double>& knotVectorV, std::vector<std::vector<XYZW>>& controlPoints)
{
	VALIDATE_ARGUMENT(!xAxis.IsZero(), "xAxis", "xAxis must not be zero vector.");
	VALIDATE_ARGUMENT(!yAxis.IsZero(), "yAxis", "yAxis must not be zero vector.")
	VALIDATE_ARGUMENT(MathUtils::IsGreaterThan(endRad, startRad), "endRad", "endRad must greater than startRad.");
	VALIDATE_ARGUMENT(MathUtils::IsGreaterThan(radius, 0.0), "radius", "Radius must greater than zero.");
	VALIDATE_ARGUMENT(MathUtils::IsGreaterThan(height, 0.0), "height", "Height must greater than zero.");

	XYZ nX = const_cast<XYZ&>(xAxis).Normalize();
	XYZ nY = const_cast<XYZ&>(yAxis).Normalize();

	int arcDegree;
	std::vector<double> arcKnotVector;
	std::vector<XYZW> arcControlPoints;
	bool isCreated = NurbsCurve::CreateArc(origin, nX, nY, startRad, endRad, radius, radius, arcDegree, arcKnotVector, arcControlPoints);
	if (!isCreated) return false;

	XYZ axis = nX.CrossProduct(nY);
	XYZ translation = height * axis;
	XYZ halfTranslation = 0.5 * height * axis;

	int size = arcControlPoints.size();
	controlPoints.resize(3,std::vector<XYZW>(size));

	for (int i = 0; i < size; i++)
	{
		double w = arcControlPoints[i].GetW();

		controlPoints[2][i] = XYZW(arcControlPoints[i].ToXYZ(true), w);
		controlPoints[1][i] = XYZW(halfTranslation + arcControlPoints[i].ToXYZ(true), w);
		controlPoints[0][i] = XYZW(translation + arcControlPoints[i].ToXYZ(true), w);
	}

	degreeU = 2;
	degreeV = arcDegree;
	knotVectorU = { 0,0,0,1,1,1 };
	knotVectorV = arcKnotVector;

	return true;
}

void LNLib::NurbsSurface::CreateRuledSurface(int degree0, const std::vector<double>& knotVector0, const std::vector<XYZW> controlPoints0, int degree1, const std::vector<double>& knotVector1, const std::vector<XYZW>& controlPoints1, int& degreeU, int& degreeV, std::vector<double>& knotVectorU, std::vector<double>& knotVectorV, std::vector<std::vector<XYZW>>& controlPoints)
{
	VALIDATE_ARGUMENT(degree0 > 0, "degree0", "Degree must greater than zero.");
	VALIDATE_ARGUMENT(knotVector0.size() > 0, "knotVector0", "KnotVector size must greater than zero.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidKnotVector(knotVector0), "knotVector0", "KnotVector must be a nondecreasing sequence of real numbers.");
	VALIDATE_ARGUMENT(controlPoints0.size() > 0, "controlPoints0", "ControlPoints must contains one point at least.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidNurbs(degree0, knotVector0.size(), controlPoints0.size()), "controlPoints0", "Arguments must fit: m = n + p + 1");

	VALIDATE_ARGUMENT(degree1 > 0, "degree1", "Degree must greater than zero.");
	VALIDATE_ARGUMENT(knotVector1.size() > 0, "knotVector1", "KnotVector size must greater than zero.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidKnotVector(knotVector1), "knotVector1", "KnotVector must be a nondecreasing sequence of real numbers.");
	VALIDATE_ARGUMENT(controlPoints1.size() > 0, "controlPoints1", "ControlPoints must contains one point at least.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidNurbs(degree1, knotVector1.size(), controlPoints1.size()), "controlPoints1", "Arguments must fit: m = n + p + 1");

	int k0Size = knotVector0.size();
	int k1Size = knotVector1.size();
	bool knotVectorCheck = MathUtils::IsAlmostEqualTo(knotVector0[0], knotVector1[0]) && MathUtils::IsAlmostEqualTo(knotVector0[k0Size - 1], knotVector1[k1Size - 1]);
	VALIDATE_ARGUMENT(knotVectorCheck, "knotVector0 & knotVector1", "Ensure that the two curves are defined on the same parameter range.");

	degreeU = std::max(degree0, degree1);
	degreeV = 1;

	std::vector<double> kv0 = knotVector0;
	std::vector<XYZW> cp0 = controlPoints0;
	if (degree0 < degreeU)
	{
		std::vector<double> updatedKnotVector0;
		std::vector<XYZW> updatedControlPoints0;
		int times = degreeU - degree0;
		NurbsCurve::ElevateDegree(degree0, kv0, cp0, times, updatedKnotVector0, updatedControlPoints0);
		kv0 = updatedKnotVector0;
		cp0 = updatedControlPoints0;
	}

	std::vector<double> kv1 = knotVector1;
	std::vector<XYZW> cp1 = controlPoints1;
	if (degree1 < degreeU)
	{
		std::vector<double> updatedKnotVector1;
		std::vector<XYZW> updatedControlPoints1;
		int times = degreeU - degree1;
		NurbsCurve::ElevateDegree(degree1, kv1, cp1, times, updatedKnotVector1, updatedControlPoints1);
		kv1 = updatedKnotVector1;
		cp1 = updatedControlPoints1;
	}
	
	knotVectorU = kv0;
	if (kv0 != kv1)
	{
		std::vector<double> insertedKnotElement0;
		std::vector<double> insertedKnotElement1;
		Polynomials::GetInsertedKnotElement(kv0, kv1, insertedKnotElement0, insertedKnotElement1);

		if (insertedKnotElement0.size() > 0)
		{
			std::vector<double> updatedKnotVector0;
			std::vector<XYZW> updatedControlPoints0;
			NurbsCurve::RefineKnotVector(degreeU, kv0, cp0, insertedKnotElement0, updatedKnotVector0, updatedControlPoints0);
			knotVectorU = updatedKnotVector0;
			cp0 = updatedControlPoints0;
		}
		if (insertedKnotElement1.size() > 0)
		{
			std::vector<double> updatedKnotVector1;
			std::vector<XYZW> updatedControlPoints1;
			NurbsCurve::RefineKnotVector(degreeU, kv1, cp1, insertedKnotElement1, updatedKnotVector1, updatedControlPoints1);
			cp1 = updatedControlPoints1;
		}
	}

	knotVectorV = { 0,0,1,1 };
	int size = cp0.size();
	controlPoints.resize(size, std::vector<XYZW>(2));
	for (int i = 0; i < size; i++)
	{
		controlPoints[i][0] = cp0[i];
		controlPoints[i][1] = cp1[i];
	}
}

bool LNLib::NurbsSurface::CreateRevolvedSurface(const XYZ& origin, const XYZ& axis, double rad, const std::vector<XYZW>& generatrixControlPoints, int& degreeU, std::vector<double>& knotVectorU, std::vector<std::vector<XYZW>>& controlPoints)
{
	VALIDATE_ARGUMENT(!axis.IsZero(), "axis", "Axis must not be zero vector.");

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

	int m = generatrixControlPoints.size() - 1;
	XYZ X, Y, O, P0, P2, T0, T2;
	double r = 0.0;
	int index = 0;

	degreeU = 2;
	controlPoints.resize(n + 1, std::vector<XYZW>(m + 1));

	for (int j = 0; j <= m; j++)
	{
		XYZW gp = generatrixControlPoints[j];
		XYZ p = gp.ToXYZ(true);

		XYZ O = Projection::PointToLine(origin, axis, p);
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

void LNLib::NurbsSurface::GlobalInterpolation(const std::vector<std::vector<XYZ>>& throughPoints, int degreeU, int degreeV, std::vector<double>& knotVectorU, std::vector<double>& knotVectorV, std::vector<std::vector<XYZW>>& controlPoints)
{
	std::vector<double> uk;
	std::vector<double> vl;

	Interpolation::GetSurfaceMeshParameterization(throughPoints, uk, vl);

	int rows = throughPoints.size();
	int cols = throughPoints[0].size();

	controlPoints.resize(rows, std::vector<XYZW>(cols));
	knotVectorU = Interpolation::ComputeKnotVector(degreeU, uk);
	knotVectorV = Interpolation::ComputeKnotVector(degreeV, vl);

	for (int j = 0; j < cols; j++)
	{
		std::vector<XYZ> temp(rows);
		for (int i = 0; i < rows; i++)
		{
			temp[i] = throughPoints[i][j];
		}
		std::vector<XYZW> cps;
		NurbsCurve::GlobalInterpolation(degreeU, temp, knotVectorU, cps);
		for (int i = 0; i < rows; i++)
		{
			controlPoints[i][j] = cps[i];
		}
	}

	for (int i = 0; i < rows; i++)
	{
		std::vector<XYZ> temp(cols);
		for (int j = 0; j < cols; j++)
		{
			temp[j] = throughPoints[i][j];
		}
		std::vector<XYZW> cps;
		NurbsCurve::GlobalInterpolation(degreeV, temp, knotVectorV, cps);
		for (int j = 0; j < cols; j++)
		{
			controlPoints[i][j] = cps[j];
		}
	}
}

bool LNLib::NurbsSurface::BicubicLocalInterpolation(const std::vector<std::vector<XYZ>>& throughPoints, std::vector<double>& knotVectorU, std::vector<double>& knotVectorV, std::vector<std::vector<XYZW>>& controlPoints)
{
	int degreeU = 3;
	int degreeV = 3;

	int row = throughPoints.size();
	int n = row - 1;
	int column = throughPoints[0].size();
	int m = column - 1;

	std::vector<std::vector<std::vector<XYZ>>> td(n + 1, std::vector<std::vector<XYZ>>(m+1, std::vector<XYZ>(3)));

	std::vector<double> ub(n + 1, 0.0);
	std::vector<double> vb(m + 1, 0.0);

	std::vector<double> r(m + 1);
	std::vector<double> s(n + 1);

	double total = 0.0;
	for (int l = 0; l <= m; l++)
	{
		std::vector<XYZ> columnData;
		MathUtils::GetColumn(throughPoints, l, columnData);

		std::vector<XYZ> tvkl;
		bool hasTangents = Interpolation::ComputeTangent(columnData, tvkl);
		if (!hasTangents) return false;

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
		std::vector<XYZ> tukl;
		bool hasTangents = Interpolation::ComputeTangent(throughPoints[k], tukl);
		if (!hasTangents) return false;

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

	knotVectorU.resize(2 * (degreeU + 1) + 2 * (n - 1));
	for (int i = 0; i <= degreeU; i++)
	{
		knotVectorU[i] = 0;
		knotVectorU[ub.size() - 1 - i] = 1;
	}
	for (int i = 1; i < n; i = i + 2)
	{
		knotVectorU[degreeU + i] = knotVectorU[degreeU + (i + 1)] = ub[i] ;
	}

	knotVectorV.resize(2 * (degreeV + 1) + 2 * (m - 1));
	for (int i = 0; i <= degreeV; i++)
	{
		knotVectorV[i] = 0;
		knotVectorV[vb.size() - 1 - i] = 1;
	}
	for (int i = 1; i < m; i = i + 2)
	{
		knotVectorV[degreeV + i] = knotVectorV[degreeV + (i + 1)] = vb[i];
	}

	std::vector<std::vector<XYZW>> bcp;
	std::vector<std::vector<XYZW>> tcp;
	for (int i = 0; i < row; i++)
	{
		std::vector<double> temp;
		NurbsCurve::CubicLocalInterpolation(throughPoints[i], temp, tcp[i]);
	}	

	std::vector<std::vector<XYZW>> transpose;
	for (int j = 0; j < tcp[0].size(); j++)
	{
		std::vector<XYZ> columnData;
		MathUtils::GetColumn(ToXYZ(tcp), j, columnData);
		std::vector<double> temp;
		NurbsCurve::CubicLocalInterpolation(columnData, temp, transpose[j]);
	}

	MathUtils::Transpose(transpose, bcp);
	std::vector<std::vector<XYZ>> bezierControlPoints = ToXYZ(bcp);

	for (int k = 0; k <= n; k++)
	{
		double ak = 0.0;
		if (k > 0)
		{
			ak = (ub[k] - ub[k - 1]) / ((ub[k] - ub[k - 1]) + (ub[k + 1] - ub[k]));
		}
		for (int l = 0; l <= m; l++)
		{
			double bl = 0.0;
			if (k == l == 0)
			{
				td[0][0][2] = XYZ(0,0,0);
				continue;
			}
			if (l > 0)
			{
				bl = (vb[l] - vb[l - 1]) / ((vb[l] - vb[l - 1]) + (vb[l + 1] - vb[l]));
			}
			
			XYZ dvukl = (1 - ak) * (td[k][l][1] - td[k - 1][l][1]) / (ub[k] - ub[k - 1]) + ak * (td[k + 1][l][1] - td[k][l][1]) / (ub[k + 1] - ub[k]);
			XYZ duvkl = (1 - bl) * (td[k][l][0] - td[k][l - 1][0]) / (vb[l] - vb[l - 1]) + bl * (td[k][l + 1][0] - td[k][l][0]) / (vb[l + 1] - vb[l]);

			td[k][l][2] = (ak * duvkl + bl * dvukl) / (ak + bl);
		}
	}

	for (int k = 0; k < n; k++)
	{
		for (int l = 0; l < m; l++)
		{
			double gamma = (ub[k + 1] - ub[k]) * (vb[l + 1] - vb[l]) / 9;
			bezierControlPoints[3 * k + 1][3 * l + 1] =  gamma * td[k][l][2]         + bezierControlPoints[3 * k][3 * l + 1]     + bezierControlPoints[3 * k + 1][3 * l]     - bezierControlPoints[3 * k][3 * l];
			bezierControlPoints[3 * k + 2][3 * l + 1] = -gamma * td[k + 1][l][2]     + bezierControlPoints[3 * k + 3][3 * l + 1] + bezierControlPoints[3 * k + 3][3 * l]     - bezierControlPoints[3 * k + 2][3 * l];
			bezierControlPoints[3 * k + 1][3 * l + 2] = -gamma * td[k][l + 1][2]     + bezierControlPoints[3 * k + 1][3 * l + 3] + bezierControlPoints[3 * k][3 * l + 3]     - bezierControlPoints[3 * k][3 * l + 2];
			bezierControlPoints[3 * k + 2][3 * l + 2] =  gamma * td[k + 1][l + 1][2] + bezierControlPoints[3 * k + 2][3 * l + 3] + bezierControlPoints[3 * k + 3][3 * l + 2] - bezierControlPoints[3 * k + 3][3 * l + 3];
		} 
	}

	controlPoints = ToXYZW(bezierControlPoints);
	return true;
}

void LNLib::NurbsSurface::GlobalSurfaceApproximation(const std::vector<std::vector<XYZ>>& throughPoints, unsigned int degreeU, unsigned int degreeV, int controlPointsRows, int controlPointsColumns, std::vector<double>& knotVectorU, std::vector<double>& knotVectorV, std::vector<std::vector<XYZW>>& controlPoints)
{
	int rows = controlPointsRows;
	int n = rows - 1;
	int columns = controlPointsColumns;
	int m = columns - 1;

	std::vector<double> uk;
	std::vector<double> vl;

	Interpolation::GetSurfaceMeshParameterization(throughPoints, uk, vl);

	int sizeU = static_cast<int>(throughPoints.size());
	int r = sizeU - 1;
	int sizeV = static_cast<int>(throughPoints[0].size());
	int s = sizeV - 1;

	knotVectorU = Interpolation::ComputeKnotVector(degreeU, sizeU, rows, uk);
	knotVectorV = Interpolation::ComputeKnotVector(degreeV, sizeV, columns, vl);

	std::vector<std::vector<double>> Nu;
	for (int i = 1; i < r; i++)
	{
		std::vector<double> temp;
		for (int j = 1; j < n; j++)
		{
			temp.emplace_back(Polynomials::OneBasisFunction(j, degreeU, knotVectorU, uk[i]));
		}
		Nu.emplace_back(temp);
	}
	std::vector<std::vector<double>> NTu;
	MathUtils::Transpose(Nu, NTu);
	std::vector<std::vector<double>> NTNu = MathUtils::MatrixMultiply(NTu, Nu);
	std::vector<std::vector<double>> NTNul;
	std::vector<std::vector<double>> NTNuu;
	MathUtils::LUDecomposition(NTNu, NTNul, NTNuu);

	/*std::vector<std::vector<XYZ>> tempControlPoints;
	tempControlPoints.resize(rows);
	for (int j = 0; j < sizeV; j++)
	{
		tempControlPoints[0][j] = throughPoints[0][j];
		tempControlPoints[n][j] = throughPoints[r][j];

		XYZ Q0 = tempControlPoints[0][j];
		XYZ Qm = tempControlPoints[n][j];

		std::vector<XYZ> Rku;
		for (int i = 1; i < r; i++)
		{
			double N0p = Polynomials::OneBasisFunction(0, degreeU, knotVectorU, uk[i]);
			double Nnp = Polynomials::OneBasisFunction(n, degreeU, knotVectorU, uk[i]);

			Rku[i] = throughPoints[i][j] - N0p * Q0 - Nnp * Qm;
		}

		std::vector<XYZ> R;
		for (int i = 1; i < n; i++)
		{
			XYZ Rk_temp;
			for (int k = 1; k < r; k++)
			{
				double Np = Polynomials::OneBasisFunction(i, degreeU, knotVectorU, uk[k]);
				Rk_temp += Np* Rku[k];
			}
			R[i] = Rk_temp;
		}

		for (int i = 0; i < 3; i++)
		{
			std::vector<double> rhs;
			for (int k = 0; k < n; k++)
			{
				rhs[k] = R[k][i];
			}
			std::vector<double> column = MathUtils::ForwardSubstitution(NTNul, rhs);
			std::vector<double> sol = MathUtils::BackwardSubstitution(NTNuu, column);

			for (int k = 1; k < n; k++)
			{
				tempControlPoints[k][j][i] = sol[k - 1];
			}
		}
	}
	
	std::vector<std::vector<double>> Nv;
	for (int i = 1; i < s; i++)
	{
		std::vector<double> temp;
		for (int j = 1; j < m; j++)
		{
			temp.emplace_back(Polynomials::OneBasisFunction(j, degreeV, knotVectorV, vl[i]));
		}
		Nv.emplace_back(temp);
	}
	std::vector<std::vector<double>> NTv;
	MathUtils::Transpose(Nv, NTv);
	std::vector<std::vector<double>> NTNv = MathUtils::MatrixMultiply(NTv, Nv);
	std::vector<std::vector<double>> NTNvl;
	std::vector<std::vector<double>> NTNvu;
	MathUtils::LUDecomposition(NTNv, NTNvl, NTNvu);

	std::vector<std::vector<XYZ>> P;
	P.resize(rows);
	for (int i = 0; i < rows; i++)
	{
		P[i].resize(columns);
	}

	for (int i = 0; i < rows; i++)
	{
		P[i][0] = tempControlPoints[i][0];
		P[i][m] = tempControlPoints[i][s];

		XYZ Q0 = P[i][0];
		XYZ Qm = P[i][m];

		std::vector<XYZ> Rkv;
		for (int j = 1; j < s; j++)
		{
			double N0p = Polynomials::OneBasisFunction(0, degreeV, knotVectorV, vl[j]);
			double Nnp = Polynomials::OneBasisFunction(m, degreeV, knotVectorV, vl[j]);

			Rkv[i] = tempControlPoints[i][j] - N0p * Q0 - Nnp * Qm;
		}

		std::vector<XYZ> R;
		for (int j = 1; j < m; j++)
		{
			XYZ Rk_temp;
			for (int k = 1; k < r; k++)
			{
				double Np = Polynomials::OneBasisFunction(i, degreeV, knotVectorV, vl[k]);
				Rk_temp += Np * Rkv[k];
			}
			R[j] = Rk_temp;
		}

		for (int j = 0; j < 3; j++)
		{
			std::vector<double> rhs;
			for (int k = 0; k < m; k++)
			{
				rhs[k] = R[i][k];
			}
			std::vector<double> column = MathUtils::ForwardSubstitution(NTNul, rhs);
			std::vector<double> sol = MathUtils::BackwardSubstitution(NTNuu, column);

			for (int k = 1; k < m; k++)
			{
				P[i][k][j] = sol[k - 1];
			}
		}
	}

	controlPoints = ToXYZW(P);*/
}


