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
		double ak = GetAk(qk_1, qk, qk1, qk2);
		return ((1 - ak) * qk + ak * qk1).Normalize();
	}
}

double LNLib::Interpolation::GetTotalChordLength(const std::vector<XYZ>& throughPoints)
{
	int n = throughPoints.size() - 1;
	double length = 0.0;
	for (int i = 1; i <= n; i++)
	{
		length += throughPoints[i].Distance(throughPoints[i - 1]);
	}
	return length;
}

std::vector<double> LNLib::Interpolation::GetChordParameterization(const std::vector<XYZ>& throughPoints)
{
	int size = throughPoints.size();
	int n = size - 1;

	std::vector<double> uk(size,0.0);
	uk[n] = 1.0;

	double d = GetTotalChordLength(throughPoints);
	for (int i = 1; i <= n - 1; i++)
	{
		uk[i] = uk[i - 1] + (throughPoints[i].Distance(throughPoints[i - 1])) / d;
	}
	return uk;
}

double LNLib::Interpolation::GetCentripetalLength(const std::vector<XYZ>& throughPoints)
{
	int size = throughPoints.size();
	int n = size - 1;

	double length = 0.0;
	for (int i = 1; i <= n; i++)
	{
		length += sqrt(throughPoints[i].Distance(throughPoints[i - 1]));
	}
	return length;
}

std::vector<double> LNLib::Interpolation::GetCentripetalParameterization(const std::vector<XYZ>& throughPoints)
{
	int size = throughPoints.size();
	int n = size - 1;

	std::vector<double> uk(size, 0.0);
	uk[n] = 1.0;

	double d = GetCentripetalLength(throughPoints);
	for (int i = 1; i <= n - 1; i++)
	{
		uk[i] = uk[i - 1] + sqrt(throughPoints[i].Distance(throughPoints[i - 1])) / d;
	}
	return uk;
}

std::vector<double> LNLib::Interpolation::GetChordParameterization(const std::vector<XYZ>& throughPoints, int startIndex, int endIndex)
{
	int size = endIndex - startIndex + 1;
	std::vector<double> uk;
	uk.resize(size, 0.0);

	double length = 0.0;
	for (int i = startIndex; i <= endIndex; i++)
	{
		length += throughPoints[i].Distance(throughPoints[i - 1]);
	}

	for (int i = startIndex; i <= endIndex; i++)
	{
		uk[i] = uk[i - 1] + (throughPoints[i].Distance(throughPoints[i - 1])) / length;
	}
	return uk;
}

std::vector<double> LNLib::Interpolation::ComputeKnotVector(int degree, const std::vector<double> params)
{
	std::vector<double> knotVector;
	std::vector<double> uk = params;

	int size = params.size();
	int n = size - 1;
	int m = n + degree + 1;

	knotVector.resize(m + 1, 0.0);
	for (int i = m - degree; i <= m; i++)
	{
		knotVector[i] = 1.0;
	}

	for (int j = 1; j <= n - degree; j++)
	{
		double sum = 0.0;
		for (int i = j; i <= j + degree - 1; i++)
		{
			sum += uk[i];
		}
		knotVector[j + degree] = (1.0 / degree) * sum;
	}
	return knotVector;
}

std::vector<double> LNLib::Interpolation::ComputeKnotVector(unsigned int degree, int pointsCount, int controlPointsCount, const std::vector<double> params)
{
	std::vector<double>  knotVector;
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
	return knotVector;
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

bool LNLib::Interpolation::GetSurfaceMeshParameterization(const std::vector<std::vector<XYZ>>& throughPoints, std::vector<double>& paramsU, std::vector<double>& paramsV)
{
	int n = throughPoints.size();
	int m = throughPoints[0].size();

	std::vector<double> cds(std::max(n, m),0.0);
	paramsU.resize(n,0.0);
	paramsV.resize(m,0.0);

	int num = m;
	
	for (int l = 0; l < m; l++)
	{
		double total = 0.0;
		for (int k = 1; k < n; k++)
		{
			cds[k] = throughPoints[k][l].Distance(throughPoints[k - 1][l]);
			total += cds[k];
		}

		if (MathUtils::IsAlmostEqualTo(total, 0.0))
		{
			num--;
		}
		else
		{
			double d = 0.0;
			for (int k = 1; k < n; k++)
			{
				d += cds[k];
				paramsU[k] = paramsU[k] + d / total;
			}
		}
	}
	if (num == 0)
	{
		return false;
	}

	for (int k = 1; k < n - 1; k++)
	{
		paramsU[k] = paramsU[k] / num;
	}
	paramsU[n - 1] = 1.0;

	num = n;

	for (int k = 0; k < n; k++) 
	{
		double total = 0.0;
		for (int l = 1; l < m; l++) 
		{
			cds[l] = throughPoints[k][l].Distance(throughPoints[k][l-1]);
			total += cds[l];
		}
		if (MathUtils::IsAlmostEqualTo(total, 0.0))
		{
			num--;
		}	
		else 
		{
			double d = 0.0;
			for (int l = 1; l < m; l++) 
			{
				d += cds[l];
				paramsV[l] += d / total;
			}
		}
	}

	if (num == 0)
	{
		return false;
	}
	for (int l = 1; l < m - 1; l++)
	{
		paramsV[l] = paramsV[l] / num;
	}
	paramsV[m - 1] = 1.0;

	return true;
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
