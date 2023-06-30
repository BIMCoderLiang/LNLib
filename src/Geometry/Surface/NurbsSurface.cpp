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
#include "Projection.h"
#include "Intersection.h"
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

	int samplesU = static_cast<int>(controlPoints.size() * degreeU);
	int samplesV = static_cast<int>(controlPoints[0].size() * degreeV);
	double spanU = (maxUParam - minUParam) / (samplesU - 1);
	double spanV = (maxVParam - minVParam) / (samplesV - 1);
	for (int i = 0; i < samplesU - 1; i++)
	{
		double currentU = minUParam + spanU * i;
		double nextU = minUParam + spanU * (i + 1);
		for (int j = 0; j < samplesV; j++)
		{
			double v = minVParam + spanV * j;
			
			XYZ currentPoint;
			NurbsSurface::GetPointOnSurface(controlPoints, knotVectorU, knotVectorV, degreeU, degreeV, UV(currentU, v), currentPoint);
			
			XYZ nextPoint;
			NurbsSurface::GetPointOnSurface(controlPoints, knotVectorU, knotVectorV, degreeU, degreeV, UV(nextU, v), nextPoint);

			XYZ vector1 = currentPoint - givenPoint;
			XYZ vector2 = nextPoint - currentPoint;
			double dot = vector1.DotProduct(vector2);

			XYZ projectPoint;
			UV project = UV(Constants::DoubleEpsilon, Constants::DoubleEpsilon);

			if (dot < 0)
			{
				projectPoint = currentPoint;
				project = UV(currentU,v);
			}
			else if (dot > 1)
			{
				projectPoint = nextPoint;
				project = UV(nextU, v);
			}
			else
			{
				projectPoint = currentPoint + dot * vector1.Normalize();
				project = UV(currentU + (nextU - currentU) * dot, v);
			}

			double distance = (givenPoint - projectPoint).Length();
			if (distance < minValue)
			{
				minValue = distance;
				param = project;
			}
		}
	}

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

bool LNLib::NurbsSurface::GetUVTangent(const std::vector<std::vector<XYZW>>& controlPoints, const std::vector<double>& knotVectorU, const std::vector<double>& knotVectorV, unsigned int degreeU, unsigned int degreeV, const UV param, const XYZ& tangent, UV& uvTangent)
{
	std::vector<std::vector<XYZ>> derivatives;
	ComputeRationalSurfaceDerivatives(controlPoints, knotVectorU, knotVectorV, degreeU, degreeV, param, 1, derivatives);
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

void LNLib::NurbsSurface::ReverseU(const std::vector<std::vector<XYZW>>& controlPoints, const std::vector<double>& knotVectorU, std::vector<double>& reversedKnotVectorU, std::vector<std::vector<XYZW>>& reversedControlPoints)
{
	NurbsCurve::ReverseKnotVector(knotVectorU, reversedKnotVectorU);
	int row = static_cast<int>(controlPoints.size());
	int column = static_cast<int>(controlPoints[0].size());

	reversedControlPoints.resize(row);
	for (int i = 0; i < row; i++)
	{
		reversedControlPoints[i].resize(column);
		for (int j = 0; j < column; j++)
		{
			reversedControlPoints[i][j] = controlPoints[i][j];
		}
	}

	std::reverse(reversedControlPoints.begin(), reversedControlPoints.end());
}

void LNLib::NurbsSurface::ReverseV(const std::vector<std::vector<XYZW>>& controlPoints, const std::vector<double>& knotVectorV, std::vector<double>& reversedKnotVectorV, std::vector<std::vector<XYZW>>& reversedControlPoints)
{
	NurbsCurve::ReverseKnotVector(knotVectorV, reversedKnotVectorV);
	int row = static_cast<int>(controlPoints.size());

	reversedControlPoints.resize(row);
	for (int i = 0; i < row; i++)
	{
		std::vector<XYZW> rowReversed;
		NurbsCurve::ReverseControlPoints(controlPoints[i], rowReversed);

		reversedControlPoints.emplace_back(rowReversed);
	}
}


void LNLib::NurbsSurface::CreateBilinearSurface(const XYZ& point0, const XYZ& point1, const XYZ& point2, const XYZ& point3, std::vector<double>& knotVectorU, std::vector<double>& knotVectorV, std::vector<std::vector<XYZW>>& controlPoints)
{
	int degree = 3;

	for (int i = 0; i <= degree; i++) 
	{
		std::vector<XYZW> row;
		for (int j = 0; j <= degree; j++) 
		{
			double l = 1.0 - i / degree;
			XYZ inter12 = l * point0 + (1 - l) * point1;
			XYZ inter34 = l * point3 + (1 - l) * point2;

			XYZ res = inter12 * (1- j/degree) + (j/degree) * inter34;
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
	XYZ xTemp = xAxis;
	XYZ yTemp = yAxis;

	XYZ nX = xTemp.Normalize();
	XYZ nY = yTemp.Normalize();

	int arcDegree;
	std::vector<double> arcKnotVector;
	std::vector<XYZW> arcControlPoints;

	bool isCreated = NurbsCurve::CreateArc(origin, nX, nY, radius, radius, startRad, endRad, arcDegree, arcKnotVector, arcControlPoints);
	if (!isCreated) return false;
	XYZ axis = nX.CrossProduct(nY);
	XYZ translation = height * axis;
	XYZ halfTranslation = 0.5 * height * axis;

	int size = static_cast<int>(arcControlPoints.size());

	controlPoints.resize(3);
	for (int i = 0; i < 3; i++)
	{
		controlPoints[i].resize(size);
	}

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

bool LNLib::NurbsSurface::CreateRuledSurface(int degree0, const std::vector<double>& knotVector0, const std::vector<XYZW> controlPoints0, int degree1, const std::vector<double>& knotVector1, const std::vector<XYZW>& controlPoints1, int& degreeU, int& degreeV, std::vector<double>& knotVectorU, std::vector<double>& knotVectorV, std::vector<std::vector<XYZW>>& controlPoints)
{
	int k0Size = static_cast<int>(knotVector0.size());
	int k1Size = static_cast<int>(knotVector1.size());
	if (k0Size == k1Size &&
		!MathUtils::IsAlmostEqualTo(knotVector0[0], knotVector1[0]) ||
		!MathUtils::IsAlmostEqualTo(knotVector0[k0Size - 1], knotVector1[k1Size - 1]))
	{
		return false;
	}

	degreeU = std::max(degree0, degree1);
	degreeV = 1;

	std::vector<double> updatedKnotVector0;
	std::vector<XYZW> updatedControlPoints0;
	std::vector<double> updatedKnotVector1;
	std::vector<XYZW> updatedControlPoints1;
	if (degree0 < degreeU)
	{
		int times = degreeU - degree0;
		NurbsCurve::ElevateDegree(degree0, knotVector0, controlPoints0, times, updatedKnotVector0, updatedControlPoints0);
	}
	else
	{
		updatedKnotVector0 = knotVector0;
		updatedControlPoints0 = controlPoints0;
	}
	if (degree1 < degreeU)
	{
		int times = degreeU - degree1;
		NurbsCurve::ElevateDegree(degree1, knotVector1, controlPoints1, times, updatedKnotVector1, updatedControlPoints1);
	}
	else
	{
		updatedKnotVector1 = knotVector1;
		updatedControlPoints1 = controlPoints1;
	}

	bool isSameKnotVector = true;
	for (int i = 0; i < k0Size; i++)
	{
		if (!MathUtils::IsAlmostEqualTo(knotVector0[i], knotVector1[i]))
		{
			isSameKnotVector = false;
			break;
		}
	}
	
	std::vector<XYZW> finalControlPoints0;
	std::vector<XYZW> finalControlPoints1;
	if (!isSameKnotVector)
	{
		std::vector<double> insertedKnotElement0;
		std::vector<double> insertedKnotElement1;

		Polynomials::GetInsertedKnotElement(updatedKnotVector0, updatedKnotVector1, insertedKnotElement0, insertedKnotElement1);

		std::vector<double> tempKnotVector0;
		std::vector<double> tempKnotVector1;
		NurbsCurve::RefineKnotVector(degreeU, updatedKnotVector0, updatedControlPoints0, insertedKnotElement0, tempKnotVector0, finalControlPoints0);
		NurbsCurve::RefineKnotVector(degreeV, updatedKnotVector1, updatedControlPoints1, insertedKnotElement1, tempKnotVector1, finalControlPoints1);

		knotVectorU = tempKnotVector0;
	}
	else
	{
		knotVectorU = knotVector0;
		finalControlPoints0 = updatedControlPoints0;
		finalControlPoints1 = updatedControlPoints1;
	}

	knotVectorV = { 0,0,1,1 };
	int size = static_cast<int>(finalControlPoints0.size());
	controlPoints.resize(size);
	for (int i = 0; i < size; i++)
	{
		controlPoints[i].resize(2);
		controlPoints[i][0] = finalControlPoints0[i];
		controlPoints[i][1] = finalControlPoints1[i];
	}

	return true;
}

bool LNLib::NurbsSurface::CreateRevolvedSurface(const XYZ& origin, const XYZ& axis, double rad, const std::vector<XYZW>& generatrixCurve, int& degreeU, std::vector<double>& knotVectorU, std::vector<std::vector<XYZW>>& controlPoints)
{
	int narcs = 0;
	if (MathUtils::IsLessThanOrEqual(rad, Constants::Pi / 2))
	{
		narcs = 1;
	}
	else
	{
		if (MathUtils::IsLessThanOrEqual(rad, Constants::Pi))
		{
			narcs = 2;
			knotVectorU.resize(6);
			knotVectorU[3] = knotVectorU[4] = 0.5;
		}
		else if (MathUtils::IsLessThanOrEqual(rad, 3 * Constants::Pi / 2))
		{
			narcs = 3;
			knotVectorU.resize(8);
			knotVectorU[3] = knotVectorU[4] = 1.0 / 3.0;
			knotVectorU[5] = knotVectorU[6] = 2.0 / 3.0;
		}
		else
		{
			narcs = 4;
			knotVectorU.resize(10);
			knotVectorU[3] = knotVectorU[4] = 0.25;
			knotVectorU[5] = knotVectorU[6] = 0.5;
			knotVectorU[7] = knotVectorU[8] = 0.75;
		}
	}

	double dtheta = rad / narcs;
	int j = 3 + 2 * (narcs - 1);
	for (int i = 0; i < 3; j++,i++)
	{
		knotVectorU[i] = 0.0;
		knotVectorU[j] = 1.0;
	}

	int n = 2 * narcs;
	double wm = cos(dtheta / 2.0);
	double angle = 0.0;
	std::vector<double> cosines(narcs + 1), sines(narcs + 1);

	for (int i = 1; i <= narcs; i++)
	{
		angle = angle + dtheta;
		cosines[i] = cos(angle);
		sines[i] = sin(angle);
	}

	int m = static_cast<int>(generatrixCurve.size());
	XYZ X, Y, O, P0, P2, T0, T2;
	double r = 0.0;
	int index = 0;

	controlPoints.resize(n + 1);
	for (int i = 0; i <= n; i++)
	{
		controlPoints[i].resize(m + 1);
	}

	for (int j = 0; j <= m; j++)
	{
		XYZW gp = generatrixCurve[j];
		XYZ p = XYZ(gp[0], gp[1], gp[2]) / gp[3];

		Projection::PointToLine(origin, axis, p, O);
		X = p - O;

		r = X.Normalize().Length();
		Y = axis.CrossProduct(X);

		P0 = p;
		controlPoints[0][j] = XYZW(P0, gp[3]);

		T0 = Y;
		index = 0;
		angle = 0.0;

		for (int i = 1; i <= narcs; i++)
		{
			P2 = O + r * cosines[i] * X + r * sines[i] * Y;
			controlPoints[index + 2][j] = XYZW(P2, gp[3]);
			T2 = -sines[i] * X + cosines[i] * Y;

			XYZ intersectPoint;
			double param0;
			double param1;
			CurveCurveIntersectionType type = Intersection::ComputeRays(P0, T0, P2, T2, param0, param1, intersectPoint);
			if (type != CurveCurveIntersectionType::Intersecting)
			{
				return false;
			}
			controlPoints[index + 1][j] = XYZW(intersectPoint, gp[3]);

			index += 2;
			if (i < narcs)
			{
				P0 = P2;
				T0 = T2;
			}
		}
	}

	degreeU = 2;

	return true;
}


