/*
 * Author:
 * 2023/06/08 - Yuqing Liang (BIMCoder Liang)
 * bim.frankliang@foxmail.com
 * 微信公众号：BIMCoder梁老师
 *
 * Use of this source code is governed by a GPL-3.0 license license that can be found in
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
#include "ValidationUtils.h"

void LNLib::NurbsSurface::GetPointOnSurface(const std::vector<std::vector<XYZW>>& controlPoints, const std::vector<double>& knotVectorU, const std::vector<double>& knotVectorV, unsigned int degreeU, unsigned int degreeV, UV uv, XYZ& point)
{
	XYZW sw = XYZW();
	BsplineSurface::GetPointOnSurface(controlPoints, knotVectorU, knotVectorV, degreeU, degreeV, uv, sw);
	point = sw.ToXYZ(true);
}


void LNLib::NurbsSurface::ComputeRationalSurfaceDerivatives(const std::vector<std::vector<XYZW>>& controlPoints, const std::vector<double>& knotVectorU, const std::vector<double>& knotVectorV, unsigned int degreeU, unsigned int degreeV, UV uv, unsigned int derivative, std::vector<std::vector<XYZ>>& derivatives)
{
	
	std::vector<std::vector<XYZW>> ders;
	BsplineSurface::ComputeDerivatives(controlPoints, knotVectorU, knotVectorV, degreeU, degreeV, uv, derivative, ders);

	std::vector<std::vector<XYZ>> Aders;
	for (int i = 0; i <= static_cast<int>(degreeU); i++)
	{
		for (int j = 0; j <= static_cast<int>(degreeV); j++)
		{
			Aders[i][j] = ders[i][j].ToXYZ(false);
		}
	}

	std::vector<std::vector<double>> wders;
	for (int i = 0; i <= static_cast<int>(degreeU); i++)
	{
		for (int j = 0; j <= static_cast<int>(degreeV); j++)
		{
			wders[i][j] = ders[i][j].GetW();
		}
	}

	for (int k = 0; k <= static_cast<int>(derivative); k++)
	{
		for (int l = 0; l <= static_cast<int>(derivative - k); l++)
		{
			XYZ v = Aders[k][l];
			for ( int j = 1; j <= l; j++)
			{
				v = v - MathUtils::Binomial(l, j) * wders[0][j] * derivatives[k][l - j];
			}

			for (int i = 1; i <= k; i++)
			{
				v = v - MathUtils::Binomial(k, i) * wders[i][0] * derivatives[k - i][l];

				XYZ v2 = XYZ();
				for (int j = 1; j <= l; j++)
				{
					v2 = v2 + MathUtils::Binomial(l, j) * wders[i][j] * derivatives[k - i][l - j];
				}

				v = v - MathUtils::Binomial(k, i) * v2;
			}

			derivatives[k][l] = v / wders[0][0];
		}
	}
}

void LNLib::NurbsSurface::InsertKnot(const std::vector<std::vector<XYZW>>& controlPoints, const std::vector<double>& knotVector, unsigned int degree, double insertKnot, unsigned int times, bool isUDirection, std::vector<double>& insertedKnotVector, std::vector<std::vector<XYZW>>& updatedControlPoints)
{
	int n = static_cast<int>(knotVector.size()) - degree - 2;
	int knotSpanIndex = Polynomials::GetKnotSpanIndex(n, degree, insertKnot, knotVector);
	int multiplicity = Polynomials::GetKnotMultiplicity(insertKnot, knotVector);

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

	for (int i = 1; i <= static_cast<int>(times); i++)
	{
		insertedKnotVector[knotSpanIndex + i] = insertKnot;
	}

	for (int i = knotSpanIndex + 1; i < static_cast<int>(knotVector.size()); i++)
	{
		insertedKnotVector[i + times] = knotVector[i];
	}

	std::vector<std::vector<double>> alpha;
	alpha.resize(degree - multiplicity);
	for (int i = 0; i < static_cast<int>(degree) - multiplicity; i++)
	{
		alpha[i].resize(times + 1);
	}

	for (int j = 1; j <= static_cast<int>(times); j++)
	{
		int L = knotSpanIndex - degree + j;
		for (int i = 0; i <= static_cast<int>(degree) - j - multiplicity; i++)
		{
			alpha[i][j] = (insertKnot - knotVector[L + i]) / (knotVector[i + knotSpanIndex + 1] - knotVector[L + i]);
		}
	}

	std::vector<XYZW> temp;
	temp.resize(degree + 1);

	int controlPointsRows = static_cast<int>(controlPoints.size());
	int controlPointsColumns = static_cast<int>(controlPoints[0].size());

	if (isUDirection)
	{
		updatedControlPoints.resize(controlPointsRows + times);
		for (int i = 0; i < controlPointsRows + static_cast<int>(times); i++)
		{
			updatedControlPoints[i].resize(controlPointsColumns);
		}

		for (int col = 0; col < controlPointsColumns; col++)
		{
			for (int i = 0; i <= knotSpanIndex - static_cast<int>(degree); i++)
			{
				updatedControlPoints[i][col] = controlPoints[i][col];
			}

			for (int i = knotSpanIndex - multiplicity; i < controlPointsRows; i++)
			{
				updatedControlPoints[i + times][col] = controlPoints[i][col];
			}

			for (int i = 0; i < static_cast<int>(degree) - multiplicity + 1; i++)
			{
				temp[i] = controlPoints[knotSpanIndex - degree + i][col];
			}


			int L = 0;
			for (int j = 1; j <= static_cast<int>(times); j++)
			{
				L = knotSpanIndex - degree + j;
				for (int i = 0; i <= static_cast<int>(degree) - j - multiplicity; i++)
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
		updatedControlPoints.resize(controlPointsRows);
		for (int i = 0; i < controlPointsRows; i++)
		{
			updatedControlPoints[i].resize(controlPointsColumns + times);
		}

		for (int row = 0; row < controlPointsRows; row++)
		{
			for (int i = 0; i <= knotSpanIndex - static_cast<int>(degree); i++)
			{
				updatedControlPoints[row][i] = controlPoints[row][i];
			}

			for (int i = knotSpanIndex - multiplicity; i < controlPointsColumns; i++)
			{
				updatedControlPoints[row][i+times] = controlPoints[row][i];
			}

			for (int i = 0; i < static_cast<int>(degree) - multiplicity + 1; i++)
			{
				temp[i] = controlPoints[row][knotSpanIndex - degree + i];
			}


			int L = 0;
			for (int j = 1; j <= static_cast<int>(times); j++)
			{
				L = knotSpanIndex - degree + j;
				for (int i = 0; i <= static_cast<int>(degree) - j - multiplicity; i++)
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

void LNLib::NurbsSurface::RefineKnotVector(const std::vector<std::vector<XYZW>>& controlPoints, const std::vector<double>& knotVectorU, const std::vector<double>& knotVectorV, unsigned int degreeU, unsigned int degreeV, std::vector<double>& insertKnotElements, bool isUDirection, std::vector<double>& insertedKnotVectorU, std::vector<double>& insertedKnotVectorV, std::vector<std::vector<XYZW>>& updatedControlPoints)
{
	std::vector<double> knots;
	unsigned int degree;
	std::vector<std::vector<XYZW>> tempControlPoints;
	std::vector<std::vector<XYZW>> newControlPoints;

	if (isUDirection)
	{
		MathUtils::Transpose(controlPoints, tempControlPoints);
		int size = static_cast<int>(knotVectorU.size());
		knots.resize(size);
		for (int i = 0; i < size; i++)
		{
			knots[i] = knotVectorU[i];
		}
		degree = degreeU;
	}
	else
	{
		tempControlPoints = controlPoints;
		int size = static_cast<int>(knotVectorV.size());
		knots.resize(size);
		for (int i = 0; i < size; i++)
		{
			knots[i] = knotVectorV[i];
		}
		degree = degreeV;
	}

	std::vector<double> insertedKnotVector;
	
	for (int i = 0; i < static_cast<int>(tempControlPoints.size()); i++)
	{
		insertedKnotVector.clear();
		std::vector<XYZW> updatedCurveControlPoints;
		NurbsCurve::RefineKnotVector(degree, knots, tempControlPoints[i], insertKnotElements, insertedKnotVector, updatedCurveControlPoints);
		newControlPoints[i] = updatedCurveControlPoints;
	}

	if (isUDirection)
	{
		insertedKnotVectorU = insertedKnotVector;
		insertedKnotVectorV = knotVectorV;
		MathUtils::Transpose(newControlPoints, updatedControlPoints);
	}
	else
	{
		insertedKnotVectorU = knotVectorU;
		insertedKnotVectorV = insertedKnotVector;
		updatedControlPoints = newControlPoints;
	}
}

void LNLib::NurbsSurface::ToBezierPatches(const std::vector<std::vector<XYZW>>& controlPoints, const std::vector<double>& knotVectorU, const std::vector<double>& knotVectorV, unsigned int degreeU, unsigned int degreeV, int& bezierPatchesCount, std::vector<std::vector<std::vector<XYZW>>>& decomposedControlPoints)
{
	int controlPointsRow = static_cast<int>(controlPoints.size());
	int conrolPointsColumn = static_cast<int>(controlPoints[0].size());

	std::vector<std::vector<std::vector<XYZW>>> temp;
	temp.resize(conrolPointsColumn);

	int ubezierCurvesCount = 0;
	for (int col = 0; col < conrolPointsColumn; col++)
	{
		std::vector<XYZW> uDirectionPoints;
		uDirectionPoints.resize(controlPointsRow);
		MathUtils::GetColumn(controlPoints, col, uDirectionPoints);

		ubezierCurvesCount = 0;
		std::vector<std::vector<XYZW>> decomposedUPoints;
		NurbsCurve::ToBezierCurves(degreeU, knotVectorU, uDirectionPoints, ubezierCurvesCount, decomposedUPoints);

		temp.emplace_back(decomposedUPoints);
	}

	int vbezierCurvesCount = 0;
	for (int i = 0; i < ubezierCurvesCount; i++)
	{
		int row = static_cast<int>(temp[0][i].size());
		for (int r = 0; r < row; r++)
		{
			std::vector<XYZW> vDirectionPoints;
			for (int j = 0; j < conrolPointsColumn; j++)
			{
				std::vector<XYZW> segement = temp[j][i];
				vDirectionPoints.emplace_back(segement[r]);
			}

			vbezierCurvesCount = 0;
			std::vector<std::vector<XYZW>> decomposedVPoints;
			NurbsCurve::ToBezierCurves(degreeV, knotVectorV, vDirectionPoints, vbezierCurvesCount, decomposedVPoints);

			for (int v = 0; v < vbezierCurvesCount; v++)
			{
				decomposedControlPoints[i * vbezierCurvesCount + v][r] = decomposedVPoints[v];
			}			
		}
	}

	bezierPatchesCount = ubezierCurvesCount * vbezierCurvesCount;
}

void LNLib::NurbsSurface::RemoveKnot(const std::vector<std::vector<XYZW>>& controlPoints, const std::vector<double>& knotVectorU, const std::vector<double>& knotVectorV, unsigned int degreeU, unsigned int degreeV, double removeKnot, unsigned int times, bool isUDirection, std::vector<double>& restKnotVectorU, std::vector<double>& restKnotVectorV, std::vector<std::vector<XYZW>>& updatedControlPoints)
{
	int controlPointsRow = static_cast<int>(controlPoints.size());
	int conrolPointsColumn = static_cast<int>(controlPoints[0].size());
	if (isUDirection)
	{
		std::vector<std::vector<XYZW>> temp;
		std::vector<double> uKnotVector;
		for (int v = 0; v < conrolPointsColumn; v++)
		{
			std::vector<XYZW> uControlPoints;
			MathUtils::GetColumn(controlPoints, v, uControlPoints);
			uKnotVector.clear();
			std::vector<XYZW> uUpdatedControlPoints;
			NurbsCurve::RemoveKnot(degreeU, knotVectorU, uControlPoints, removeKnot, times, uKnotVector, uUpdatedControlPoints);
			temp.emplace_back(uUpdatedControlPoints);
		}

		MathUtils::Transpose(temp, updatedControlPoints);
		restKnotVectorU = uKnotVector;
		restKnotVectorV = knotVectorV;

	}
	else
	{
		std::vector<double> vKnotVector;
		for (int u = 0; u < controlPointsRow; u++)
		{
			vKnotVector.clear();
			std::vector<XYZW> vUpdatedControlPoints;
			NurbsCurve::RemoveKnot(degreeV, knotVectorV, controlPoints[u], removeKnot, times, vKnotVector, vUpdatedControlPoints);
			updatedControlPoints.emplace_back(vUpdatedControlPoints);
		}

		restKnotVectorU = knotVectorU;
		restKnotVectorV = vKnotVector;
	}
}



void LNLib::NurbsSurface::ElevateDegree(const std::vector<std::vector<XYZW>>& controlPoints, const std::vector<double>& knotVectorU, const std::vector<double>& knotVectorV, unsigned int degreeU, unsigned int degreeV, unsigned int times, bool isUDirection, std::vector<double>& updatedKnotVectorU, std::vector<double>& updatedKnotVectorV, std::vector<std::vector<XYZW>>& updatedControlPoints)
{
	std::vector<double> knots;
	unsigned int degree;
	std::vector<std::vector<XYZW>> tempControlPoints;
	std::vector<std::vector<XYZW>> newControlPoints;

	if (isUDirection)
	{
		MathUtils::Transpose(controlPoints, tempControlPoints);
		int size = static_cast<int>(knotVectorU.size());
		knots.resize(size);
		for (int i = 0; i < size; i++)
		{
			knots[i] = knotVectorU[i];
		}
		degree = degreeU;
	}
	else
	{
		tempControlPoints = controlPoints;
		int size = static_cast<int>(knotVectorV.size());
		knots.resize(size);
		for (int i = 0; i < size; i++)
		{
			knots[i] = knotVectorV[i];
		}
		degree = degreeV;
	}

	std::vector<double> tempKnotVector;

	for (int i = 0; i < static_cast<int>(tempControlPoints.size()); i++)
	{
		tempKnotVector.clear();
		std::vector<XYZW> updatedCurveControlPoints;
		NurbsCurve::ElevateDegree(degree, knots, tempControlPoints[i], times, tempKnotVector, updatedCurveControlPoints);
		newControlPoints[i] = updatedCurveControlPoints;
	}

	if (isUDirection)
	{
		updatedKnotVectorU = tempKnotVector;
		updatedKnotVectorV = knotVectorV;
		MathUtils::Transpose(newControlPoints, updatedControlPoints);
	}
	else
	{
		updatedKnotVectorU = knotVectorU;
		updatedKnotVectorV = tempKnotVector;
		updatedControlPoints = newControlPoints;
	}
}

bool LNLib::NurbsSurface::ReduceDegree(const std::vector<std::vector<XYZW>>& controlPoints, const std::vector<double>& knotVectorU, const std::vector<double>& knotVectorV, unsigned int degreeU, unsigned int degreeV, bool isUDirection, std::vector<double>& updatedKnotVectorU, std::vector<double>& updatedKnotVectorV, std::vector<std::vector<XYZW>>& updatedControlPoints)
{
	std::vector<double> knots;
	unsigned int degree;
	std::vector<std::vector<XYZW>> tempControlPoints;
	std::vector<std::vector<XYZW>> newControlPoints;

	if (isUDirection)
	{
		MathUtils::Transpose(controlPoints, tempControlPoints);
		int size = static_cast<int>(knotVectorU.size());
		knots.resize(size);
		for (int i = 0; i < size; i++)
		{
			knots[i] = knotVectorU[i];
		}
		degree = degreeU;
	}
	else
	{
		tempControlPoints = controlPoints;
		int size = static_cast<int>(knotVectorV.size());
		knots.resize(size);
		for (int i = 0; i < size; i++)
		{
			knots[i] = knotVectorV[i];
		}
		degree = degreeV;
	}

	std::vector<double> tempKnotVector;

	for (int i = 0; i < static_cast<int>(tempControlPoints.size()); i++)
	{
		tempKnotVector.clear();
		std::vector<XYZW> updatedCurveControlPoints;
		bool result  = NurbsCurve::ReduceDegree(degree, knots, tempControlPoints[i], tempKnotVector, updatedCurveControlPoints);
		if (!result) return false;
		newControlPoints[i] = updatedCurveControlPoints;
	}

	if (isUDirection)
	{
		updatedKnotVectorU = tempKnotVector;
		updatedKnotVectorV = knotVectorV;
		MathUtils::Transpose(newControlPoints, updatedControlPoints);
	}
	else
	{
		updatedKnotVectorU = knotVectorU;
		updatedKnotVectorV = tempKnotVector;
		updatedControlPoints = newControlPoints;
	}
	return true;
}

LNLib::UV LNLib::NurbsSurface::GetParamOnSurface(const std::vector<std::vector<XYZW>>& controlPoints, const std::vector<double>& knotVectorU, const std::vector<double>& knotVectorV, unsigned int degreeU, unsigned int degreeV, const XYZ& givenPoint)
{
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

	//to do... tessllate surface.

	int counters = 0;
	while (counters < maxIterations)
	{
		std::vector<std::vector<XYZ>> derivatives;
		ComputeRationalSurfaceDerivatives(controlPoints,knotVectorU,knotVectorV,degreeU,degreeV,param,2,derivatives);
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

		//to do... solve J*delta = k
		UV temp = UV(0, 0);

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


