/*
 * Author:
 * 2023/07/04 - Yuqing Liang (BIMCoder Liang)
 * bim.frankliang@foxmail.com
 * 微信公众号：BIMCoder梁老师
 *
 * Use of this source code is governed by a GPL-3.0 license that can be found in
 * the LICENSE file.
 */

#include "Interpolation.h"
#include "XYZ.h"
#include "MathUtils.h"
#include "Polynomials.h"
#include <algorithm>

namespace LNLib
{
	XYZ GetQk(const std::vector<XYZ>& throughPoints, unsigned int index)
	{
		return throughPoints[index] - throughPoints[index - 1];
	}

	double GetAk(const XYZ& qk_1, const XYZ& qk, const XYZ& qk1, const XYZ& qk2)
	{
		return (qk_1.CrossProduct(qk)).Length() / ((qk_1.CrossProduct(qk).Length()) + (qk1.CrossProduct(qk2)).Length());
	}

	XYZ GetVk(const XYZ& qk_1, const XYZ& qk, const XYZ& qk1, const XYZ& qk2)
	{
		double ak = GetAk(qk_1,qk,qk1,qk2);
		return ((1 - ak) * qk + ak * qk1).Normalize();
	}
}

double LNLib::Interpolation::GetTotalChordLength(const std::vector<XYZ>& throughPoints)
{
	int size = static_cast<int>(throughPoints.size());
	int n = size - 1;

	double result = 0.0;
	for (int i = 1; i <= n; i++)
	{
		result += throughPoints[i].Distance(throughPoints[i - 1]);
	}
	return result;
}

std::vector<double> LNLib::Interpolation::GetChordParameterization(const std::vector<XYZ>& throughPoints)
{
	int size = static_cast<int>(throughPoints.size());
	int n = size - 1;

	std::vector<double> uk;
	uk.resize(size,0.0);
	uk[n] = 1.0;

	double d = GetTotalChordLength(throughPoints);
	for (int i = 1; i <= n - 1; i++)
	{
		uk[i] = uk[i - 1] + (throughPoints[i].Distance(throughPoints[i - 1])) / d;
	}
	return uk;
}

void LNLib::Interpolation::ComputeKnotVector(unsigned int degree, int pointsCount, const std::vector<double> params, std::vector<double>& knotVector)
{
	std::vector<double> uk = params;

	int size = pointsCount;
	int n = size - 1;
	int m = n + degree + 1;

	knotVector.resize(m + 1, 0.0);
	for (int i = m - degree; i <= m; i++)
	{
		knotVector[i] = 1.0;
	}

	for (int j = 1; j <= static_cast<int>(n - degree + 1); j++)
	{
		double temp = 0.0;
		for (int i = j; i <= static_cast<int>(j + degree - 1); i++)
		{
			temp += uk[i];
		}
		knotVector[j + degree] = (1.0 / degree) * temp;
	}
}

void LNLib::Interpolation::ComputeKnotVector(unsigned int degree, int pointsCount, int controlPointsCount, const std::vector<double> params, std::vector<double>& knotVector)
{
	for (int i = 0; i <= degree; i++)
	{
		knotVector[i] = 0;
	}
	double d = pointsCount / (controlPointsCount - degree);
	for (int j = 1; j < controlPointsCount - degree; j++)
	{
		int i = static_cast<int>(j * d);
		double alpha = (j * d) - i;
		double temp = (1.0 - alpha) * params[i - 1] + (alpha * params[i]);
		knotVector.emplace_back(temp);
	}
	for (int i = 0; i <= degree; i++)
	{
		knotVector.emplace_back(1.0);
	}
}

std::vector<std::vector<double>> LNLib::Interpolation::MakeInterpolationMatrix(unsigned int degree, int dataCount, const std::vector<double>& params, const std::vector<double>& knotVector)
{
	int n = dataCount - 1;
	std::vector<std::vector<double>> A;
	A.resize(dataCount);
	for (int i = 0; i <= n; i++)
	{
		A[i].resize(dataCount);
	}

	for (int i = 0; i <= n; i++)
	{
		int spanIndex = Polynomials::GetKnotSpanIndex(n, degree, params[i], knotVector);
		std::vector<double> basis;
		Polynomials::BasisFunctions(spanIndex, degree, knotVector[i], knotVector, basis);

		for (int k = 0; k < static_cast<int>(degree); k++)
		{
			A[i][spanIndex - degree + k] = basis[k];
		}
	}
	return A;
}

std::vector<LNLib::XYZ> LNLib::Interpolation::ComputerMatrixMultiplyPoints(std::vector<std::vector<double>> matrix, std::vector<XYZ> points)
{
	int row = static_cast<int>(matrix.size());
	int column = static_cast<int>(matrix[0].size());
	int size = static_cast<int>(points.size());
	std::vector<XYZ> result(row);
	if (!(column == size))
	{
		return result;
	}
	for (int i = 0; i <= row; i++)
	{
		XYZ temp = XYZ(0, 0, 0);
		for (int j = 0; j <= column; j++)
		{
			temp += matrix[i][j] * points[j];
		}
		result.emplace_back(temp);
	}
	return result;
}

std::vector<LNLib::XYZ> LNLib::Interpolation::ComputerControlPointsByLUDecomposition(const std::vector<std::vector<double>>& matrix, const std::vector<XYZ>& data)
{
	int size = static_cast<int>(data.size());
	int n = size - 1;

	std::vector<XYZ> tempControlPoints;
	tempControlPoints.resize(size);

	std::vector<std::vector<double>> matrixL;
	std::vector<std::vector<double>> matrixU;
	MathUtils::LUDecomposition(matrix, matrixL, matrixU);

	for (int i = 0; i < 3; i++)
	{
		std::vector<double> rhs;
		for (int j = 0; j <= n; j++)
		{
			rhs[j] = data[j][i];
		}
		std::vector<double> column = MathUtils::ForwardSubstitution(matrixL, rhs);
		std::vector<double> sol = MathUtils::BackwardSubstitution(matrixU, column);

		for (int j = 0; j <= n; j++)
		{
			tempControlPoints[j][i] = sol[j];
		}
	}
	return tempControlPoints;
}

void LNLib::Interpolation::ComputerKnotVectorForTangents(unsigned int degree, const std::vector<double>& params, const std::vector<int>& derivativeIndices, std::vector<double>& knotVector)
{
	int paramsCount = static_cast<int>(params.size());
	int derCount = static_cast<int>(derivativeIndices.size());

	int startIndex;
	int endIndex;

	int internalDerCount = derCount;

	if (derivativeIndices[0] == 0)
	{
		startIndex = 0;
		--internalDerCount;
	}
	else
	{
		startIndex = 1;
	}

	if (derivativeIndices[derCount - 1] == paramsCount - 1)
	{
		endIndex = paramsCount - degree;
		--internalDerCount;
	}
	else
	{
		endIndex = paramsCount - degree - 1;
	}

	int averageKnotsCount = endIndex - startIndex + 1 + 2;
	std::vector<double> averageKnotVector(averageKnotsCount);

	averageKnotVector[0] = params[0];

	int newKnotIndex = 0;
	for (int i = startIndex; i <= endIndex; i++)
	{
		double sum = 0.0;
		for (int t = 0; t < static_cast<int>(degree); t++)
		{
			sum += params[i + t];
		}
		double mean = sum / degree;
		averageKnotVector[++newKnotIndex] = mean;
	}
	averageKnotVector[++newKnotIndex] = params[paramsCount - 1];

	newKnotIndex = -1;

	std::vector<double> tempParams(paramsCount + 1);
	std::vector<double> derKnotVector(internalDerCount);

	for (int i = 0; i < averageKnotsCount - 1; ++i)
	{
		tempParams[0] = averageKnotVector[i];
		std::vector<double> paramIndices(paramsCount);
		int tempParamIndex = 1;
		for (int j = 0; j < paramsCount; ++j)
		{
			double parameter = params[j];
			if (parameter >= averageKnotVector[i] && parameter < averageKnotVector[i + 1])
			{
				paramIndices[tempParamIndex] = j;
				tempParams[tempParamIndex++] = parameter;
			}
		}
		tempParams[tempParamIndex] = averageKnotVector[i + 1];

		for (int j = 0; j <= tempParamIndex - 2; ++j)
		{
			for (int k = 0; k < derivativeIndices.size(); ++k)
			{
				if (paramIndices[j + 1] == derivativeIndices[k] &&
					paramIndices[j + 1] != 0 &&
					paramIndices[j + 1] != paramsCount - 1)
				{
					double sum = 0.0;
					for (int t = 0; t < static_cast<int>(degree); t++)
					{
						sum += params[j + t];
					}
					double mean = sum / degree;
					derKnotVector[++newKnotIndex] = mean;
					break;
				}
			}
		}
	}

	std::vector<double> tempKnotVector(averageKnotVector.size() + derKnotVector.size());
	std::merge(averageKnotVector.data(), averageKnotVector.data() + averageKnotVector.size(),
		derKnotVector.data(), derKnotVector.data() + derKnotVector.size(),
		tempKnotVector.data());

	int numKnots = paramsCount + derCount + degree + 1;
	knotVector.resize(numKnots);

	for (int i = 0; i < static_cast<int>(degree); i++)
	{
		knotVector[i] = tempKnotVector[0];
		knotVector[numKnots - 1 - i] = tempKnotVector[tempKnotVector.size() - 1];
	}

	for (int i = 0; i < static_cast<int>(tempKnotVector.size()); i++)
	{
		knotVector[degree + i] = tempKnotVector[i];
	}
}



void LNLib::Interpolation::GetSurfaceMeshParameterization(const std::vector<std::vector<XYZ>>& throughPoints, std::vector<double>& paramVectorU, std::vector<double>& paramVectorV)
{
	int sizeU = static_cast<int>(throughPoints.size());
	int sizeV = static_cast<int>(throughPoints[0].size());

	int n = sizeU - 1;
	int m = sizeV - 1;

	int num = m + 1;
	paramVectorU.resize(sizeU, 0.0);
	paramVectorU[n] = 1.0;
	
	for (int l = 0; l <= m; l++)
	{
		double total = 0.0;
		std::vector<double> cds;
		cds.resize(n + 1, 0.0);

		for (int k = 1; k <= n; k++)
		{
			cds[k] = throughPoints[k][l].Distance(throughPoints[k - 1][l]);
			total += cds[k];
		}

		if (MathUtils::IsAlmostEqualTo(total, 0.0))
		{
			num = num - 1;
		}
		else
		{
			double d = 0.0;
			for (int k = 1; k < n; k++)
			{
				d += cds[k];
				paramVectorU[k] = paramVectorU[k] + d / total;
			}
		}
	}
	for (int k = 1; k < n; k++)
	{
		paramVectorU[k] = paramVectorU[k] / num;
	}

	num = n + 1;
	paramVectorV.resize(sizeV, 0.0);
	paramVectorV[m] = 1.0;

	for (int k = 0; k <= n; k++)
	{
		double total = 0.0;
		std::vector<double> cds;
		cds.resize(m + 1, 0.0);

		for (int l = 1; l <= m; l++)
		{
			cds[l] = throughPoints[k][l].Distance(throughPoints[k][l - 1]);
			total += cds[l];
		}

		if (MathUtils::IsAlmostEqualTo(total, 0.0))
		{
			num = num - 1;
		}
		else
		{
			double d = 0.0;
			for (int l = 1; l < m; l++)
			{
				d += cds[l];
				paramVectorV[l] = paramVectorU[l] + d / total;
			}
		}
	}
	for (int l = 1; l < m; l++)
	{
		paramVectorV[l] = paramVectorV[l] / num;
	}
}

bool LNLib::Interpolation::ComputerTangent(const std::vector<XYZ>& throughPoints, std::vector<XYZ>& tangents)
{
	int size = static_cast<int>(throughPoints.size());
	if (size < 5) return false;
	int n = size - 1;
	tangents.resize(size);
	for (int k = 2; k <= n - 2 ; k++)
	{
		LNLib::XYZ qk_1 = GetQk(throughPoints, k - 1);
		LNLib::XYZ qk = GetQk(throughPoints, k);
		LNLib::XYZ qk1 = GetQk(throughPoints, k+1);
		LNLib::XYZ qk2 = GetQk(throughPoints, k+2);

		tangents[k] = GetVk(qk_1, qk, qk1, qk2);
	}

	LNLib::XYZ q_1 = 2 * GetQk(throughPoints, 0) - GetQk(throughPoints, 1);
	LNLib::XYZ q0 = 2*GetQk(throughPoints, 1) - GetQk(throughPoints, 2);
	LNLib::XYZ qn1 = 2 * GetQk(throughPoints, n) - GetQk(throughPoints, n-1);
	LNLib::XYZ qn2 = 2 * GetQk(throughPoints, n+1) - GetQk(throughPoints, n);

	tangents[0] = GetVk(q_1, q0, GetQk(throughPoints, 1), GetQk(throughPoints, 2));
	tangents[1] = GetVk(q0, GetQk(throughPoints, 1), GetQk(throughPoints, 2), GetQk(throughPoints, 3));
	tangents[n-1] = GetVk(GetQk(throughPoints, n - 2), GetQk(throughPoints, n-1), GetQk(throughPoints, n),qn1);
	tangents[n] = GetVk(GetQk(throughPoints, n-1), GetQk(throughPoints, n),qn1,qn2);
	return true;
}
