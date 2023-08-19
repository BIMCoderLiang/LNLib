/*
 * Author:
 * 2023/06/08 - Yuqing Liang (BIMCoder Liang)
 * bim.frankliang@foxmail.com
 * 微信公众号：BIMCoder梁老师
 *
 * Use of this source code is governed by a GPL-3.0 license that can be found in
 * the LICENSE file.
 */

#include "Polynomials.h"
#include "UV.h"
#include "ValidationUtils.h"
#include "MathUtils.h"
#include "LNLibExceptions.h"
#include <algorithm>

using namespace LNLib;

namespace LNLib
{
	struct CustomDoubleEqual
	{
		bool operator()(const double& d1, const double& d2) const
		{
			return MathUtils::IsAlmostEqualTo(d1, d2);
		}
	};
	std::unordered_map<double, int, std::hash<double>, CustomDoubleEqual> GetKnotMultiplicityMap(const std::vector<double>& knotVector)
	{
		std::unordered_map<double, int, std::hash<double>, CustomDoubleEqual> result;

		for (int i = 0; i < static_cast<int>(knotVector.size()); i++)
		{
			auto got = result.find(i);
			if (got == result.end())
			{
				result.insert({ i, LNLib::Polynomials::GetKnotMultiplicity(knotVector[i], knotVector) });
			}
		}
		return result;
	}
}



double Polynomials::Horner(const std::vector<double>& coefficients, unsigned int degree, double paramT)
{
	VALIDATE_ARGUMENT(degree > 0, "degree", "Degree must greater than zero.");
	VALIDATE_ARGUMENT(degree < static_cast<int>(coefficients.size()), "degree", "Degree must less than coefficients size.");
	double result = coefficients[degree];
	for (int i = degree - 1; i >= 0; i--)
	{
		result = result * paramT + coefficients[i];
	}
	return result;
}

double Polynomials::Bernstein(unsigned int index, unsigned int degree, double paramT)
{
	if (index < 0 || 
		index > degree)
	{
		return 0.0;
	}
	if (index == 0 ||
		index == degree)
	{
		return 1.0;
	}

	VALIDATE_ARGUMENT(degree > 0, "degree", "Degree must greater than zero.");
	VALIDATE_ARGUMENT_RANGE(paramT, 0.0, 1.0);
	
	std::vector<double> temp(degree + 1,0.0);
	temp[degree - index] = 1.0;
	double t1 = 1.0 - paramT;

	for (unsigned int k = index; k <= degree; k++)
	{
		for (unsigned int j = degree; j >= k; j--)
		{
			temp[j] = t1 * temp[j] + paramT * temp[j - 1];
		}
	}
	return temp[degree];
}

std::vector<double> Polynomials::AllBernstein(unsigned int degree, double paramT)
{
	VALIDATE_ARGUMENT(degree > 0, "degree", "Degree must greater than zero.");
	VALIDATE_ARGUMENT_RANGE(paramT, 0.0, 1.0);

	std::vector<double> bernsteinArray(degree + 1);
	bernsteinArray[0] = 1.0;

	double t1 = 1.0 - paramT;
	for (unsigned int j = 1; j <= degree; j++)
	{
		double saved = 0.0;
		for (unsigned int k = 0; k < j; k++)
		{
			double temp = bernsteinArray[k];
			bernsteinArray[k] = saved + t1 * temp;
			saved = paramT * temp;
		}
		bernsteinArray[j] = saved;
	}
	return bernsteinArray;
}

double Polynomials::Horner(const std::vector<std::vector<double>>& coefficients, unsigned int n, unsigned int m, UV& uv)
{
	std::vector<double> temp;
	temp.resize(n + 1);

	for (unsigned int i = 0; i <= n; i++)
	{
		temp[i] = Horner(coefficients[i], m, uv.GetV());
	}
	return Horner(temp, n, uv.GetU());
}

int LNLib::Polynomials::GetKnotSpanIndex(unsigned int n, unsigned int degree, double paramT, const std::vector<double>& knotVector)
{

	if (MathUtils::IsGreaterThan(paramT, knotVector[n + 1]))
	{
		return n;
	}

	if (MathUtils::IsLessThan(paramT, knotVector[degree]))
	{
		return degree;
	}

	int low = degree;
	int high = n + 1;
	int mid = static_cast<int>(floor((low + high) / 2.0));

	while (paramT < knotVector[mid] || paramT >= knotVector[mid + 1])
	{
		if (paramT < knotVector[mid])
			high = mid;
		else
			low = mid;
		mid = static_cast<int>(floor((low + high) / 2.0));
	}
	return mid;
}

int LNLib::Polynomials::GetKnotMultiplicity(double knot, const std::vector<double>& knotVector)
{
	int size = static_cast<int>(knotVector.size());
	int multi = 0;

	for (int index = 0; index < size; index++)
	{
		if (MathUtils::IsAlmostEqualTo(knot, knotVector[index]))
		{
			multi++;
		}	
	}

	return multi;
}

void LNLib::Polynomials::GetInsertedKnotElement(const std::vector<double> knotVector0, const std::vector<double> knotVector1, std::vector<double>& insertElements0, std::vector<double>& insertElements1)
{
	std::unordered_map<double, int, std::hash<double>, CustomDoubleEqual> map0 = GetKnotMultiplicityMap(knotVector0);
	std::unordered_map<double, int, std::hash<double>, CustomDoubleEqual> map1 = GetKnotMultiplicityMap(knotVector1);

	for (auto it = map0.begin(); it != map0.end(); ++it)
	{
		double key0 = it->first;
		int count0 = it->second;

		auto got = map1.find(key0);
		if (got == map1.end())
		{
			for (int i = 0; i < count0; i++)
			{
				insertElements1.emplace_back(key0);
			}
		}
		else
		{
			int count1 = got->second;
			if (count0 > count1)
			{
				int times = count0 - count1;
				for (int j = 0; j < times; j++)
				{
					insertElements1.emplace_back(key0);
				}
			}
			else
			{
				int times = count1 - count0;
				for (int j = 0; j < times; j++)
				{
					insertElements0.emplace_back(key0);
				}
			}
		}
	}

	for (auto it = map1.begin(); it != map1.end(); ++it)
	{
		double key1 = it->first;
		int count1 = it->second;

		auto got = map0.find(key1);
		if (got == map0.end())
		{
			for (int i = 0; i < count1; i++)
			{
				insertElements0.emplace_back(key1);
			}
		}
	}

	std::sort(insertElements0.begin(), insertElements0.end());
	std::sort(insertElements1.begin(), insertElements1.end());
}


void LNLib::Polynomials::BasisFunctions(unsigned int spanIndex, unsigned int degree, double paramT, const std::vector<double>& knotVector, std::vector<double>& basisFunctions)
{
	basisFunctions.resize(degree + 1);
	double saved = 0.0;
	double temp = 0.0;
	
	std::vector<double> left;
	left.resize(degree + 1);
	std::vector<double> right;
	right.resize(degree + 1);

	basisFunctions[0] = 1.0;

	for (unsigned int j = 1; j <= degree; j++)
	{		
		left[j] = paramT - knotVector[spanIndex + 1 - j];
		right[j] = knotVector[spanIndex + j] - paramT;

		saved = 0.0;

		for (unsigned int r = 0; r < j; r++)
		{
			temp = basisFunctions[r] / (right[r + 1] + left[j - r]);
			basisFunctions[r] = saved + right[r + 1] * temp;
			saved = left[j - r] * temp;
		}
		basisFunctions[j] = saved;
	}
}

void LNLib::Polynomials::BasisFunctionsDerivatives(unsigned int spanIndex, unsigned int degree, double paramT,  unsigned int derivative, const std::vector<double>& knotVector, std::vector<std::vector<double>>& derivatives)
{
	derivatives.resize(derivative+1);
	for (int i = 0; i <= static_cast<int>(derivative); i++)
	{
		derivatives.resize(degree + 1);
	}

	std::vector<std::vector<double>> ndu;
	ndu.resize(degree + 1);
	for (unsigned int i = 0; i <= degree; i++)
	{
		ndu[i].resize(degree + 1);
	}

	ndu[0][0] = 1.0;

	std::vector<double> left;
	left.resize(degree + 1);
	std::vector<double> right;
	right.resize(degree + 1);

	double saved = 0.0;
	double temp = 0.0;

	for (int j = 1; j <= static_cast<int>(degree); j++)
	{
		saved = 0.0;

		left[j] = paramT - knotVector[spanIndex + 1 - j];
		right[j] = knotVector[spanIndex + j] - paramT;

		for (int r = 0; r < j; r++)
		{
			ndu[j][r] = right[r + 1] + left[j - r];
			temp = ndu[r][j - 1] / ndu[j][r];

			ndu[r][j] = saved + right[r + 1] * temp;
			saved = left[j - r] * temp;
		}
		ndu[j][j] = saved;
	}

	for (int j = 0; j <= static_cast<int>(degree); j++)
	{
		derivatives[0][j] = ndu[j][degree];
	}


	std::vector<std::vector<double>> a;
	a.resize(2);
	for (int i = 0; i < 2; i++)
	{
		a[i].resize(degree + 1);
	}

	for (int r = 0; r <= static_cast<int>(degree); r++)
	{
		int s1 = 0; 
		int s2 = 1;
		a[0][0] = 1.0;

		for (int k = 1; k <= static_cast<int>(derivative); k++)
		{
			double d = 0.0;
			int rk = r - k;
			int pk = static_cast<int>(degree - k);

			if (r >= k)
			{
				a[s2][0] = a[s1][0] / ndu[pk + 1][rk];
				d = a[s2][0] * ndu[rk][pk];
			}

			int j1 = 0;
			int j2 = 0;

			if (rk >= -1)
				j1 = 1;
			else
				j1 = -rk;

			if (r - 1 <= pk)
				j2 = k - 1;
			else
				j2 = degree - r;

			for (int j = j1; j <= j2; j++)
			{
				a[s2][j] = (a[s1][j] - a[s1][j - 1]) / ndu[pk + 1][rk + j];
				d += a[s2][j] * ndu[rk + j][pk];
			}
			if (r <= pk)
			{
				a[s2][k] = -a[s1][k - 1] / ndu[pk + 1][r];
				d += a[s2][k] * ndu[r][pk];
			}
			derivatives[k][r] = d;

			int temp = s1; 
			s1 = s2; 
			s2 = temp;
		}
	}

	int r = static_cast<int>(degree);
	for (int k = 1; k <= static_cast<int>(derivative); k++)
	{
		for (int j = 0; j <= static_cast<int>(degree); j++)
		{
			derivatives[k][j] *= r;
		}
		r *= static_cast<int>(degree - k);
	}
}

double LNLib::Polynomials::OneBasisFunction(unsigned int index, unsigned int degree, const std::vector<double>& knotVector, double paramT)
{
	int m = static_cast<int>(knotVector.size()) - 1;
	if ((index == 0 &&
		MathUtils::IsAlmostEqualTo(paramT, knotVector[0])) ||
		(index == m - degree - 1 &&
		MathUtils::IsAlmostEqualTo(paramT, knotVector[m]))
		)
	{
		return 1.0;
	}
	if (MathUtils::IsLessThan(paramT, knotVector[index]) ||
		MathUtils::IsGreaterThanOrEqual(paramT, knotVector[index + degree + 1]))
	{
		return 0.0;
	}

	std::vector<double> temp;
	temp.resize(degree + 1);

	for (int j = 0; j <= static_cast<int>(degree); j++)
	{
		temp[j] = MathUtils::IsGreaterThanOrEqual(paramT, knotVector[index + j]) && 
					MathUtils::IsLessThan(paramT, knotVector[index + j + 1]) ? 
						1.0 : 0.0;
	}
	for (int k = 1; k <= static_cast<int>(degree); k++)
	{
		double saved;
		if (MathUtils::IsAlmostEqualTo(temp[0], 0.0))
		{
			saved = 0.0;
		}
		else
		{
			saved = ((paramT - knotVector[index] * temp[0])) / (knotVector[index + k] - knotVector[index]);
		}
		for (int j = 0; j < static_cast<int>(degree) - k + 1; j++)
		{
			double knotLeft = knotVector[index + j + 1];
			double knotRight = knotVector[index + j + k + 1];
			if (MathUtils::IsAlmostEqualTo(temp[j + 1], 0.0))
			{
				temp[j] = saved;
				saved = 0.0;
			}
			else
			{
				double tempValue = temp[j + 1] / (knotRight - knotLeft);
				temp[j] = saved + (knotRight - paramT) * tempValue;
				saved = (paramT - knotLeft) * tempValue;
			}
		}
	}
	return temp[0];
}

void LNLib::Polynomials::OneBasisFunctionDerivative(unsigned int index, unsigned int degree, const std::vector<double>& knotVector, double paramT, unsigned int derivative, std::vector<double>& derivatives)
{
	derivatives.resize(derivative + 1);

	if (MathUtils::IsLessThan(paramT, knotVector[index]) ||
		MathUtils::IsGreaterThanOrEqual(paramT, knotVector[index + degree + 1]))
	{
		for (unsigned int k = 0; k <= derivative; k++)
		{
			derivatives[k] = 0.0;
		}
		return;
	}

	std::vector<std::vector<double>> basisFunctions;
	basisFunctions.resize(degree + 1);

	for (int i = 0; i <= static_cast<int>(degree); i++)
	{
		basisFunctions[i].resize(degree + 1);
	}

	for (int j = 0; j <= static_cast<int>(degree); j++)
	{
		basisFunctions[j][0] = MathUtils::IsGreaterThanOrEqual(paramT, knotVector[index + j]) &&
								MathUtils::IsLessThan(paramT, knotVector[index + j + 1]) ?
								  1.0 : 0.0;
	}
	for (int k = 1; k <= static_cast<int>(degree); k++)
	{
		double saved;
		if (MathUtils::IsAlmostEqualTo(basisFunctions[0][k - 1],0.0))
		{
			saved = 0.0;
		}
		else
		{
			saved = ((paramT - knotVector[index] * basisFunctions[0][k - 1])) / (knotVector[index + k] - knotVector[index]);
		}
		for (int j = 0; j < static_cast<int>(degree) - k + 1; j++)
		{
			double knotLeft = knotVector[index + j + 1];
			double knotRight = knotVector[index + j + k + 1];

			if (MathUtils::IsAlmostEqualTo(basisFunctions[j + 1][k - 1], 0.0))
			{
				basisFunctions[j][k] = saved;
				saved = 0.0;
			}
			else
			{
				double temp = basisFunctions[j + 1][k - 1] / (knotRight - knotLeft);
				basisFunctions[j][k] = saved + (knotRight - paramT) * temp;
				saved = (paramT - knotLeft) * temp;
			}
		}
	}

	derivatives[0] = basisFunctions[0][degree];

	std::vector<double> basisFunctionsDerivative;
	basisFunctionsDerivative.resize(derivative);

	for (int k = 1; k <= static_cast<int>(derivative); k++)
	{
		for (int j = 0; j <= k; j++)
		{
			basisFunctionsDerivative[j] = basisFunctions[j][degree - k];
		}
		for (int jj = 1; jj <= k; jj++)
		{
			double saved;
			if (MathUtils::IsAlmostEqualTo(basisFunctionsDerivative[0], 0.0))
			{
				saved = 0.0;
			}
			else
			{
				saved = basisFunctionsDerivative[0] / (knotVector[index + degree - k + jj]) - knotVector[index];
			}

			for (int j = 0; j < k - jj + 1; j++)
			{
				double knotLeft = knotVector[index + j + 1];
				double knotRight = knotVector[index + j + degree + jj + 1];

				if (MathUtils::IsAlmostEqualTo(basisFunctionsDerivative[j + 1], 0.0))
				{
					basisFunctionsDerivative[j] = (degree - k + jj) * saved;
					saved = 0.0;
				}
				else
				{
					double temp = basisFunctionsDerivative[j + 1] / (knotRight - knotLeft);
					basisFunctionsDerivative[j] = (degree - k + jj) * (saved - temp);
					saved = temp;
				}
			}
		}

		derivatives[k] = basisFunctionsDerivative[0];
	}
}

void LNLib::Polynomials::BezierToPowerMatrix(unsigned int degree, std::vector<std::vector<double>>& matrix)
{
	matrix.resize(degree+1);
	for (int i = 0; i < static_cast<int>(degree); i++)
	{
		matrix[i].resize(degree+1);
		for (int j = i + 1; j <= static_cast<int>(degree); j++)
		{
			matrix[i][j] = 0.0;
		}
	}

	matrix[0][0] = matrix[degree][degree] = 1.0;
	matrix[degree][0] = degree % 2 == 0 ? -1.0 : 1.0;

	double sign = -1.0;
	for (int i = 1; i < static_cast<int>(degree); i++)
	{
		matrix[i][i] = MathUtils::Binomial(degree,i);
		matrix[i][0] = matrix[degree][degree - 1] = sign * matrix[i][i];
		sign = -sign;
	}

	int k1 = (degree + 1) / 2;
	int pk = degree - 1;
	for (int k = 1; k < k1; k++)
	{
		sign = -1.0;
		for (int j = k + 1; j <= pk; j++)
		{
			matrix[j][k] = matrix[pk][degree - j] = sign * MathUtils::Binomial(degree, k) * MathUtils::Binomial(degree - k, j - k);
			sign = -sign;
		}
		pk = pk - 1;
	}
}

void LNLib::Polynomials::PowerToBezierMatrix(unsigned int degree, const std::vector<std::vector<double>>& matrix, std::vector<std::vector<double>>& inverseMatrix)
{
	inverseMatrix.resize(degree + 1);
	for (int i = 0; i < static_cast<int>(degree); i++)
	{
		inverseMatrix[i].resize(degree + 1);
		for (int j = i + 1; j <= static_cast<int>(degree); j++)
		{
			inverseMatrix[i][j] = 0.0;
		}
	}

	for (int i = 0; i <= static_cast<int>(degree); i++)
	{
		inverseMatrix[i][0] = inverseMatrix[degree][i] = 1.0;
		inverseMatrix[i][i] = 1.0 / (matrix[i][i]);
	}

	int k1 = (degree + 1) / 2;
	int pk = degree - 1;

	for (int k = 1; k < k1; k++)
	{
		for (int j = k + 1; j < pk; j++)
		{
			double d = 0.0;
			for (int i = k; i < j; i++)
			{
				d = d - matrix[j][i] * inverseMatrix[i][k];
			}
			inverseMatrix[j][k] = d / (inverseMatrix[j][j]);
			inverseMatrix[pk][degree - j] = inverseMatrix[j][k];
		}
		pk = pk - 1;
	}
}


