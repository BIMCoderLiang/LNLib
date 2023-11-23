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

	std::vector<std::vector<XYZ>> ToXYZ(const std::vector<std::vector<XYZW>>& surfacePoints)
	{
		int row = surfacePoints.size();
		int column = surfacePoints[0].size();

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
		int row = surfacePoints.size();
		int column = surfacePoints[0].size();

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
		KnotVectorUtils::GetInsertedKnotElement(kv0, kv1, insertedKnotElement0, insertedKnotElement1);

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
	VALIDATE_ARGUMENT(throughPoints.size() > 0, "throughPoints", "ThroughPoints row size must greater than zero.");
	VALIDATE_ARGUMENT(throughPoints[0].size() > 0, "throughPoints", "ThroughPoints column size must greater than zero.");
	VALIDATE_ARGUMENT(degreeU > 0, "degreeU", "DegreeU must greater than zero.");
	VALIDATE_ARGUMENT(degreeV > 0, "degreeV", "DegreeV must greater than zero.");

	std::vector<double> uk;
	std::vector<double> vl;

	Interpolation::GetSurfaceMeshParameterization(throughPoints, uk, vl);

	int rows = throughPoints.size();
	int cols = throughPoints[0].size();

	controlPoints.resize(rows, std::vector<XYZW>(cols));
	knotVectorU = Interpolation::AverageKnotVector(degreeU, uk);
	knotVectorV = Interpolation::AverageKnotVector(degreeV, vl);

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
	VALIDATE_ARGUMENT(throughPoints.size() > 0, "throughPoints", "ThroughPoints row size must greater than zero.");
	VALIDATE_ARGUMENT(throughPoints[0].size() > 0, "throughPoints", "ThroughPoints column size must greater than zero.");

	int degreeU = 3;
	int degreeV = 3;

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
	controlPoints = ToXYZW(rowFilter);
	return true;
}

bool LNLib::NurbsSurface::GlobalApproximation(const std::vector<std::vector<XYZ>>& throughPoints, int degreeU, int degreeV, int controlPointsRows, int controlPointsColumns, std::vector<double>& knotVectorU, std::vector<double>& knotVectorV, std::vector<std::vector<XYZW>>& controlPoints)
{
	VALIDATE_ARGUMENT(throughPoints.size() > 0, "throughPoints", "ThroughPoints row size must greater than zero.");
	VALIDATE_ARGUMENT(throughPoints[0].size() > 0, "throughPoints", "ThroughPoints column size must greater than zero.");
	VALIDATE_ARGUMENT(degreeU > 0, "degreeU", "DegreeU must greater than zero.");
	VALIDATE_ARGUMENT(degreeV > 0, "degreeV", "DegreeV must greater than zero.");

	int rows = controlPointsRows;
	int n = rows - 1;
	int columns = controlPointsColumns;
	int m = columns - 1;

	std::vector<std::vector<XYZ>> tempControlPoints;
	for (int i = 0; i < rows; i++)
	{
		std::vector<double> tempKv;
		std::vector<XYZW> tempCps;
		bool result = NurbsCurve::LeastSquaresApproximation(degreeU, throughPoints[i], rows, knotVectorU, tempCps);
		if (!result) return false;
		std::vector<XYZ> points = ControlPointsUtils::ToXYZ(tempCps);
		tempControlPoints.emplace_back(points);
	}

	std::vector<std::vector<XYZ>> preControlPoints;
	std::vector<std::vector<XYZ>> tPoints;
 	for (int i = 0; i < columns; i++)
	{
		std::vector<double> tempKv;
		std::vector<XYZW> tempCps;
		std::vector<XYZ> c = MathUtils::GetColumn(tempControlPoints, i);
		bool result = NurbsCurve::LeastSquaresApproximation(degreeV, c, columns, knotVectorV, tempCps);
		if (!result) return false;
		std::vector<XYZ> points = ControlPointsUtils::ToXYZ(tempCps);
		tPoints.emplace_back(points);
	}
	MathUtils::Transpose(tPoints, preControlPoints);
	controlPoints = ToXYZW(preControlPoints);
	return true;
}

bool LNLib::NurbsSurface::CreateLoftSurface(const std::vector<LN_Curve>& sections, LN_Surface& surface)
{
	int degree_max = 0;
	for (int i = 0; i < sections.size(); i++)
	{
		LN_Curve current = sections[i];
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

	for (int i = 0; i < sections.size(); i++)
	{
		LN_Curve current = sections[i];
		if (degree_max > current.Degree)
		{
			int times = degree_max - current.Degree;
			std::vector<double> newKv;
			std::vector<XYZW> newCps;
			NurbsCurve::ElevateDegree(current.Degree, current.KnotVector, current.ControlPoints, times, newKv, newCps);
			current.Degree = degree_max;
			current.KnotVector = newKv; 
			current.ControlPoints = newCps;
		}
	}

	// to be continued....
	return true;
}

bool LNLib::NurbsSurface::CreateSweepSurface(const LN_Curve& path, const std::vector<LN_Curve>& profiles, LN_Surface& surface)
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
	std::vector<double> parametersAlongPath(profilesSize);
	double delta = (path_max - path_min) / profilesSize;
	for (int i = 0; i < profilesSize; i++)
	{
		parametersAlongPath[i] = path_min + i * delta;
	}

	// to be continued....
	return true;
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
	}

	auto n0 = nurbs[0]; auto n2 = nurbs[2];
	int degree_n0 = n0.Degree;
	int degree_n2 = n2.Degree;
	if (degree_n0 > degree_n2)
	{
		int times = degree_n0 - degree_n2;
		std::vector<double> newKv;
		std::vector<XYZW> newCps;
		NurbsCurve::ElevateDegree(n2.Degree, n2.KnotVector, n2.ControlPoints, times, newKv, newCps);
		n2.Degree = degree_n0;
		n2.KnotVector = newKv; 
		n2.ControlPoints = newCps;
	}
	if (degree_n2 > degree_n0)
	{
		int times = degree_n2 - degree_n0;
		std::vector<double> newKv;
		std::vector<XYZW> newCps;
		NurbsCurve::ElevateDegree(n0.Degree, n0.KnotVector, n0.ControlPoints, times, newKv, newCps);
		n0.Degree = degree_n2;
		n0.KnotVector = newKv; 
		n0.ControlPoints = newCps;
	}
	if (n0.KnotVector != n2.KnotVector)
	{
		std::vector<double> insertedKnotElement0;
		std::vector<double> insertedKnotElement2;
		KnotVectorUtils::GetInsertedKnotElement(n0.KnotVector, n2.KnotVector, insertedKnotElement0, insertedKnotElement2);

		if (insertedKnotElement0.size() > 0)
		{
			std::vector<double> updatedKnotVector0;
			std::vector<XYZW> updatedControlPoints0;
			NurbsCurve::RefineKnotVector(n0.Degree, n0.KnotVector, n0.ControlPoints, insertedKnotElement0, updatedKnotVector0, updatedControlPoints0);
			n0.KnotVector = updatedKnotVector0;
			n0.ControlPoints = updatedControlPoints0;
		}
		if (insertedKnotElement2.size() > 0)
		{
			std::vector<double> updatedKnotVector2;
			std::vector<XYZW> updatedControlPoints2;
			NurbsCurve::RefineKnotVector(n2.Degree, n2.KnotVector, n2.ControlPoints, insertedKnotElement2, updatedKnotVector2, updatedControlPoints2);
			n2.KnotVector = updatedKnotVector2;
			n2.ControlPoints = updatedControlPoints2;
		}
	}

	auto n1 = nurbs[1]; auto n3 = nurbs[3];
	int degree_n1 = n1.Degree;
	int degree_n3 = n3.Degree;
	if (degree_n1 > degree_n3)
	{
		int times = degree_n1 - degree_n3;
		std::vector<double> newKv;
		std::vector<XYZW> newCps;
		NurbsCurve::ElevateDegree(n3.Degree, n3.KnotVector, n3.ControlPoints, times, newKv, newCps);
		n3.Degree = degree_n1;
		n3.KnotVector = newKv; 
		n3.ControlPoints = newCps;
	}
	if (degree_n3 > degree_n1)
	{
		int times = degree_n3 - degree_n1;
		std::vector<double> newKv;
		std::vector<XYZW> newCps;
		NurbsCurve::ElevateDegree(n1.Degree, n1.KnotVector, n1.ControlPoints, times, newKv, newCps);
		n1.Degree = degree_n3;
		n1.KnotVector = newKv; 
		n1.ControlPoints = newCps;
	}
	if (n1.KnotVector != n3.KnotVector)
	{
		std::vector<double> insertedKnotElement1;
		std::vector<double> insertedKnotElement3;
		KnotVectorUtils::GetInsertedKnotElement(n1.KnotVector, n3.KnotVector, insertedKnotElement1, insertedKnotElement3);

		if (insertedKnotElement1.size() > 0)
		{
			std::vector<double> updatedKnotVector1;
			std::vector<XYZW> updatedControlPoints1;
			NurbsCurve::RefineKnotVector(n1.Degree, n1.KnotVector, n1.ControlPoints, insertedKnotElement1, updatedKnotVector1, updatedControlPoints1);
			n1.KnotVector = updatedKnotVector1; n1.ControlPoints = updatedControlPoints1;
		}
		if (insertedKnotElement3.size() > 0)
		{
			std::vector<double> updatedKnotVector3;
			std::vector<XYZW> updatedControlPoints3;
			NurbsCurve::RefineKnotVector(n3.Degree, n3.KnotVector, n3.ControlPoints, insertedKnotElement3, updatedKnotVector3, updatedControlPoints3);
			n3.KnotVector = updatedKnotVector3; n3.ControlPoints = updatedControlPoints3;
		}
	}

	{
		std::vector<double> updatedKnotVector0;
		std::vector<XYZW> updatedControlPoints0;
		NurbsCurve::Reverse(n0.KnotVector, n0.ControlPoints, updatedKnotVector0, updatedControlPoints0);
		n0.KnotVector = updatedKnotVector0; n0.ControlPoints = updatedControlPoints0;
	}

	{
		std::vector<double> updatedKnotVector3;
		std::vector<XYZW> updatedControlPoints3;
		NurbsCurve::Reverse(n3.KnotVector, n3.ControlPoints, updatedKnotVector3, updatedControlPoints3);
		n3.KnotVector = updatedKnotVector3; n3.ControlPoints = updatedControlPoints3;
	}

	LN_Surface ruledSurface0;
	{
		int degree_u;
		int degree_v;
		std::vector<double> kv_u;
		std::vector<double> kv_v;
		std::vector<std::vector<XYZW>> cps;
		CreateRuledSurface(n0.Degree, n0.KnotVector, n0.ControlPoints, n2.Degree, n2.KnotVector, n2.ControlPoints, degree_u, degree_v, kv_u, kv_v, cps);
		ruledSurface0.DegreeU = degree_u; ruledSurface0.DegreeV = degree_v;
		ruledSurface0.KnotVectorU = kv_u; ruledSurface0.KnotVectorV = kv_v;
		ruledSurface0.ControlPoints = cps;
	}
	
	LN_Surface ruledSurface1;
	{
		int degree_u;
		int degree_v;
		std::vector<double> kv_u;
		std::vector<double> kv_v;
		std::vector<std::vector<XYZW>> cps;
		CreateRuledSurface(n1.Degree, n1.KnotVector, n1.ControlPoints, n3.Degree, n3.KnotVector, n3.ControlPoints, degree_u, degree_v, kv_u, kv_v, cps);
		ruledSurface1.DegreeU = degree_v; ruledSurface1.DegreeV = degree_u;
		ruledSurface1.KnotVectorU = kv_v; ruledSurface1.KnotVectorV = kv_u;
		std::vector<std::vector<XYZW>> transposedControlPoints;
		MathUtils::Transpose(cps, transposedControlPoints);
		ruledSurface1.ControlPoints = transposedControlPoints;
	}

	LN_Surface bilinearSurface;
	{
		XYZ point1 = NurbsCurve::GetPointOnCurve(n0.Degree, n0.KnotVector, n0.KnotVector[0], n0.ControlPoints);
		XYZ point2 = NurbsCurve::GetPointOnCurve(n2.Degree, n2.KnotVector, n2.KnotVector[0], n2.ControlPoints);
		XYZ point3 = NurbsCurve::GetPointOnCurve(n0.Degree, n0.KnotVector, n0.KnotVector[n0.KnotVector.size() - 1], n0.ControlPoints);
		XYZ point4 = NurbsCurve::GetPointOnCurve(n2.Degree, n2.KnotVector, n2.KnotVector[n2.KnotVector.size() - 1], n2.ControlPoints);
		int degree_u;
		int degree_v;
		std::vector<double> kv_u;
		std::vector<double> kv_v;
		std::vector<std::vector<XYZW>> cps;
		CreateBilinearSurface(point1, point2, point3, point4, degree_u, degree_v, kv_u, kv_v, cps);
		bilinearSurface.DegreeU = degree_u; bilinearSurface.DegreeV = degree_v;
		bilinearSurface.KnotVectorU = kv_u; bilinearSurface.KnotVectorV = kv_v;
		bilinearSurface.ControlPoints = cps;
	}

	{
		int rs0_degreeU = ruledSurface0.DegreeU;
		int rs1_degreeU = ruledSurface1.DegreeU;
		int bs_degreeU = bilinearSurface.DegreeU;

		int degreeU = std::max(rs0_degreeU, std::max(rs1_degreeU, bs_degreeU));
		if (degreeU > rs0_degreeU)
		{
			int times = degreeU - rs0_degreeU;
			std::vector<double> kv_u;
			std::vector<double> kv_v;
			std::vector<std::vector<XYZW>> cps;
			ElevateDegree(ruledSurface0.DegreeU, ruledSurface0.DegreeV, ruledSurface0.KnotVectorU, ruledSurface0.KnotVectorV, ruledSurface0.ControlPoints, times, true, kv_u, kv_v, cps);
			ruledSurface0.DegreeU = degreeU;
			ruledSurface0.KnotVectorU = kv_u;
			ruledSurface0.ControlPoints = cps;
		}

		if (degreeU > rs1_degreeU)
		{
			int times = degreeU - rs1_degreeU;
			std::vector<double> kv_u;
			std::vector<double> kv_v;
			std::vector<std::vector<XYZW>> cps;
			ElevateDegree(ruledSurface1.DegreeU, ruledSurface1.DegreeV, ruledSurface1.KnotVectorU, ruledSurface1.KnotVectorV, ruledSurface1.ControlPoints, times, true, kv_u, kv_v, cps);
			ruledSurface1.DegreeU = degreeU;
			ruledSurface1.KnotVectorU = kv_u;
			ruledSurface1.ControlPoints = cps;
		}

		if (degreeU > bs_degreeU)
		{
			int times = degreeU - bs_degreeU;
			std::vector<double> kv_u;
			std::vector<double> kv_v;
			std::vector<std::vector<XYZW>> cps;
			ElevateDegree(bilinearSurface.DegreeU, bilinearSurface.DegreeV, bilinearSurface.KnotVectorU, bilinearSurface.KnotVectorV, bilinearSurface.ControlPoints, times, true, kv_u, kv_v, cps);
			bilinearSurface.DegreeU = degreeU;
			bilinearSurface.KnotVectorU = kv_u;
			bilinearSurface.ControlPoints = cps;
		}

		int rs0_degreeV = ruledSurface0.DegreeV;
		int rs1_degreeV = ruledSurface1.DegreeV;
		int bs_degreeV = bilinearSurface.DegreeV;

		int degreeV = std::max(rs0_degreeV, std::max(rs1_degreeV, bs_degreeV));
		if (degreeV > rs0_degreeV)
		{
			int times = degreeV - rs0_degreeV;
			std::vector<double> kv_u;
			std::vector<double> kv_v;
			std::vector<std::vector<XYZW>> cps;
			ElevateDegree(ruledSurface0.DegreeU, ruledSurface0.DegreeV, ruledSurface0.KnotVectorU, ruledSurface0.KnotVectorV, ruledSurface0.ControlPoints, times, false, kv_u, kv_v, cps);
			ruledSurface0.DegreeV = degreeV;
			ruledSurface0.KnotVectorV = kv_v;
			ruledSurface0.ControlPoints = cps;
		}

		if (degreeV > rs1_degreeV)
		{
			int times = degreeV - rs1_degreeV;
			std::vector<double> kv_u;
			std::vector<double> kv_v;
			std::vector<std::vector<XYZW>> cps;
			ElevateDegree(ruledSurface1.DegreeU, ruledSurface1.DegreeV, ruledSurface1.KnotVectorU, ruledSurface1.KnotVectorV, ruledSurface1.ControlPoints, times, false, kv_u, kv_v, cps);
			ruledSurface1.DegreeV = degreeV;
			ruledSurface1.KnotVectorV = kv_v;
			ruledSurface1.ControlPoints = cps;
		}

		if (degreeV > bs_degreeV)
		{
			int times = degreeV - bs_degreeV;
			std::vector<double> kv_u;
			std::vector<double> kv_v;
			std::vector<std::vector<XYZW>> cps;
			ElevateDegree(bilinearSurface.DegreeU, bilinearSurface.DegreeV, bilinearSurface.KnotVectorU, bilinearSurface.KnotVectorV, bilinearSurface.ControlPoints, times, false, kv_u, kv_v, cps);
			bilinearSurface.DegreeV = degreeV;
			bilinearSurface.KnotVectorV = kv_v;
			bilinearSurface.ControlPoints = cps;
		}
	}
	
	{
		if (ruledSurface0.KnotVectorU != ruledSurface1.KnotVectorU)
		{
			std::vector<double> insert1;
			std::vector<double> insert2;
			KnotVectorUtils::GetInsertedKnotElement(ruledSurface0.KnotVectorU, ruledSurface1.KnotVectorU, insert1, insert2);

			if (insert1.size() > 0)
			{
				std::vector<double> kv_u;
				std::vector<double> kv_v;
				std::vector<std::vector<XYZW>> cps;
				RefineKnotVector(ruledSurface0.DegreeU, ruledSurface0.DegreeV, ruledSurface0.KnotVectorU, ruledSurface0.KnotVectorV, ruledSurface0.ControlPoints, insert1, true, kv_u, kv_v, cps);
				ruledSurface0.KnotVectorU = kv_u;
				ruledSurface0.ControlPoints = cps;
			}
			if (insert2.size() > 0)
			{
				std::vector<double> kv_u;
				std::vector<double> kv_v;
				std::vector<std::vector<XYZW>> cps;
				RefineKnotVector(ruledSurface1.DegreeU, ruledSurface1.DegreeV, ruledSurface1.KnotVectorU, ruledSurface1.KnotVectorV, ruledSurface1.ControlPoints, insert2, true, kv_u, kv_v, cps);
				ruledSurface1.KnotVectorU = kv_u;
				ruledSurface1.ControlPoints = cps;
			}
		}

		if (ruledSurface1.KnotVectorU != bilinearSurface.KnotVectorU)
		{
			std::vector<double> insert1;
			std::vector<double> insert2;
			KnotVectorUtils::GetInsertedKnotElement(ruledSurface1.KnotVectorU, bilinearSurface.KnotVectorU, insert1, insert2);

			if (insert1.size() > 0)
			{
				std::vector<double> kv_u1;
				std::vector<double> kv_v1;
				std::vector<std::vector<XYZW>> cps1;
				RefineKnotVector(ruledSurface1.DegreeU, ruledSurface1.DegreeV, ruledSurface1.KnotVectorU, ruledSurface1.KnotVectorV, ruledSurface1.ControlPoints, insert1, true, kv_u1, kv_v1, cps1);
				ruledSurface1.KnotVectorU = kv_u1;
				ruledSurface1.ControlPoints = cps1;

				std::vector<double> kv_u0;
				std::vector<double> kv_v0;
				std::vector<std::vector<XYZW>> cps0;
				RefineKnotVector(ruledSurface0.DegreeU, ruledSurface0.DegreeV, ruledSurface0.KnotVectorU, ruledSurface0.KnotVectorV, ruledSurface0.ControlPoints, insert1, true, kv_u0, kv_v0, cps0);
				ruledSurface0.KnotVectorU = kv_u0;
				ruledSurface0.ControlPoints = cps0;
			}
			if (insert2.size() > 0)
			{
				std::vector<double> kv_u;
				std::vector<double> kv_v;
				std::vector<std::vector<XYZW>> cps;
				RefineKnotVector(bilinearSurface.DegreeU, bilinearSurface.DegreeV, bilinearSurface.KnotVectorU, bilinearSurface.KnotVectorV, bilinearSurface.ControlPoints, insert2, true, kv_u, kv_v, cps);
				bilinearSurface.KnotVectorU = kv_u;
				bilinearSurface.ControlPoints = cps;
			}
		}

		if (ruledSurface0.KnotVectorV != ruledSurface1.KnotVectorV)
		{
			std::vector<double> insert1;
			std::vector<double> insert2;
			KnotVectorUtils::GetInsertedKnotElement(ruledSurface0.KnotVectorV, ruledSurface1.KnotVectorV, insert1, insert2);

			if (insert1.size() > 0)
			{
				std::vector<double> kv_u;
				std::vector<double> kv_v;
				std::vector<std::vector<XYZW>> cps;
				RefineKnotVector(ruledSurface0.DegreeU, ruledSurface0.DegreeV, ruledSurface0.KnotVectorU, ruledSurface0.KnotVectorV, ruledSurface0.ControlPoints, insert1, false, kv_u, kv_v, cps);
				ruledSurface0.KnotVectorV = kv_v;
				ruledSurface0.ControlPoints = cps;
			}
			if (insert2.size() > 0)
			{
				std::vector<double> kv_u;
				std::vector<double> kv_v;
				std::vector<std::vector<XYZW>> cps;
				RefineKnotVector(ruledSurface1.DegreeU, ruledSurface1.DegreeV, ruledSurface1.KnotVectorU, ruledSurface1.KnotVectorV, ruledSurface1.ControlPoints, insert2, false, kv_u, kv_v, cps);
				ruledSurface1.KnotVectorV = kv_v;
				ruledSurface1.ControlPoints = cps;
			}
		}

		if (ruledSurface1.KnotVectorV != bilinearSurface.KnotVectorV)
		{
			std::vector<double> insert1;
			std::vector<double> insert2;
			KnotVectorUtils::GetInsertedKnotElement(ruledSurface1.KnotVectorV, bilinearSurface.KnotVectorV, insert1, insert2);

			if (insert1.size() > 0)
			{
				std::vector<double> kv_u1;
				std::vector<double> kv_v1;
				std::vector<std::vector<XYZW>> cps1;
				RefineKnotVector(ruledSurface1.DegreeU, ruledSurface1.DegreeV, ruledSurface1.KnotVectorU, ruledSurface1.KnotVectorV, ruledSurface1.ControlPoints, insert1, false, kv_u1, kv_v1, cps1);
				ruledSurface1.KnotVectorV = kv_v1;
				ruledSurface1.ControlPoints = cps1;

				std::vector<double> kv_u0;
				std::vector<double> kv_v0;
				std::vector<std::vector<XYZW>> cps0;
				RefineKnotVector(ruledSurface0.DegreeU, ruledSurface0.DegreeV, ruledSurface0.KnotVectorU, ruledSurface0.KnotVectorV, ruledSurface0.ControlPoints, insert1, false, kv_u0, kv_v0, cps0);
				ruledSurface0.KnotVectorV = kv_v0;
				ruledSurface0.ControlPoints = cps0;
			}
			if (insert2.size() > 0)
			{
				std::vector<double> kv_u;
				std::vector<double> kv_v;
				std::vector<std::vector<XYZW>> cps;
				RefineKnotVector(bilinearSurface.DegreeU, bilinearSurface.DegreeV, bilinearSurface.KnotVectorU, bilinearSurface.KnotVectorV, bilinearSurface.ControlPoints, insert2, false, kv_u, kv_v, cps);
				bilinearSurface.KnotVectorV = kv_v;
				bilinearSurface.ControlPoints = cps;
			}
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


