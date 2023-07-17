/*
 * Author:
 * 2023/07/04 - Yuqing Liang (BIMCoder Liang)
 * bim.frankliang@foxmail.com
 * 微信公众号：BIMCoder梁老师
 *
 * Use of this source code is governed by a GPL-3.0 license license that can be found in
 * the LICENSE file.
 */

#include "Interpolation.h"
#include "XYZ.h"
#include "MathUtils.h"
#include "Polynomials.h"
#include <algorithm>

namespace LNLib
{
	XYZ& GetQk(const std::vector<XYZ>& throughPoints, unsigned int index)
	{
		return throughPoints[index] - throughPoints[index - 1];
	}

	double GetAk(const XYZ& qk_1, const XYZ& qk, const XYZ& qk1, const XYZ& qk2)
	{
		return (qk_1.CrossProduct(qk)).Length() / ((qk_1.CrossProduct(qk).Length()) + (qk1.CrossProduct(qk2)).Length());
	}

	XYZ& GetVk(const XYZ& qk_1, const XYZ& qk, const XYZ& qk1, const XYZ& qk2)
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

void LNLib::Interpolation::ComputeKnotVector(unsigned int degree, const int pointsCount, const std::vector<double> params, std::vector<double>& knotVector)
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

bool LNLib::Interpolation::LUDecomposition(const std::vector<std::vector<double>>& matrix, std::vector<std::vector<double>>& matrixL, std::vector<std::vector<double>>& matrixU)
{
	int row = static_cast<int>(matrix.size());
	if (row <= 0) return false;
	int column = static_cast<int>(matrix[0].size());
	if (row != column) return false;

	matrixL.resize(row);
	for (int i = 0; i < row; i++)
	{
		matrixL[i].resize(row, 0.0);
	}

	matrixU.resize(row);
	for (int i = 0; i < row; i++)
	{
		matrixU[i].resize(row, 0.0);
	}

	for (int i = 0; i < row; i++)
	{
		for (int k = 0; k < row; k++)
		{
			double temp = 0.0;
			for (int j = 0; j < i; j++)
			{
				temp += matrixL[i][j] * matrixU[j][k];
			}
			matrixU[i][k] = matrix[i][k] - temp;

			if (i == k)
			{
				matrixL[i][i] = 1.0;
			}
			else
			{
				temp = 0.0;
				for (int j = 0; j < i; j++)
				{
					temp += matrixL[k][j] * matrixU[j][i];
				}
				matrixL[k][i] = matrix[k][i] - temp;

				if (MathUtils::IsAlmostEqualTo(matrixU[i][i], 0.0))
				{
					matrixL[k][i] = 0.0;
				}
				else
				{
					matrixL[k][i] /= matrixU[i][i];
				}
			}
		}
	}

	return true;
}

std::vector<double> LNLib::Interpolation::ForwardSubstitution(const std::vector<std::vector<double>>& matrixL, const std::vector<double>& column)
{
	int size = static_cast<int>(column.size());

	std::vector<double> result;
	result.resize(size, 0.0);

	result[0] = column[0] / matrixL[0][0];
	for (int i = 1; i < size; i++)
	{
		double temp = 0.0;
		for (int j = 0; j < i; j++)
		{
			temp += matrixL[i][j] * result[j];
		}

		result[i] = (column[i] - temp)/ matrixL[i][i];
	}
	return result;
}

std::vector<double> LNLib::Interpolation::BackwardSubstitution(const std::vector<std::vector<double>>& matrixU, const std::vector<double>& column)
{
	int size = static_cast<int>(column.size());
	int n = size - 1;
	std::vector<double> result;
	result.resize(size, 0.0);

	result[n] = column[n] / matrixU[n][n];
	for (int i = size - 2; i >= 0; i--)
	{
		double temp = 0.0;
		for (int j = i; j < size; j++)
		{
			temp += matrixU[i][j] * result[j];
		}

		result[i] = (column[i] - temp) / matrixU[i][i];
	}
	return result;
}

std::vector<LNLib::XYZ> LNLib::Interpolation::GetSolvedMatrix(const std::vector<std::vector<double>>& matrix, const std::vector<XYZ>& data)
{
	int size = static_cast<int>(data.size());
	int n = size - 1;

	std::vector<XYZ> tempControlPoints;
	tempControlPoints.resize(size);

	std::vector<std::vector<double>> matrixL;
	std::vector<std::vector<double>> matrixU;
	Interpolation::LUDecomposition(matrix, matrixL, matrixU);

	for (int i = 0; i < 3; i++)
	{
		std::vector<double> rhs;
		for (int j = 0; j <= n; j++)
		{
			rhs[j] = data[j][i];
		}
		std::vector<double> column = Interpolation::ForwardSubstitution(matrixL, rhs);
		std::vector<double> sol = Interpolation::BackwardSubstitution(matrixU, column);

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

bool ComputerTangent(const std::vector<LNLib::XYZ>& throughPoints, std::vector<LNLib::XYZ>& tangents)
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
