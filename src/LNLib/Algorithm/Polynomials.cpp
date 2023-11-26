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
#include "MathUtils.h"
#include "ValidationUtils.h"
#include "LNLibExceptions.h"
#include <algorithm>

using namespace LNLib;

double Polynomials::Horner(int degree, const std::vector<double>& coefficients, double paramT)
{
	VALIDATE_ARGUMENT(degree > 0, "degree", "Degree must greater than zero.");
	VALIDATE_ARGUMENT(degree + 1 == static_cast<int>(coefficients.size()), "degree", "Coefficients size equals degree plus one.");
	double result = coefficients[degree];
	for (int i = degree - 1; i >= 0; i--)
	{
		result = result * paramT + coefficients[i];
	}
	return result;
}

double Polynomials::Bernstein(int index, int degree, double paramT)
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

	VALIDATE_ARGUMENT(degree >= 0, "degree", "Degree must greater than or equals zero.");
	VALIDATE_ARGUMENT_RANGE(paramT, 0.0, 1.0);
	
	std::vector<double> temp(degree + 1,0.0);
	temp[degree - index] = 1.0;
	double t1 = 1.0 - paramT;

	for (int k = index; k <= degree; k++)
	{
		for (int j = degree; j >= k; j--)
		{
			temp[j] = t1 * temp[j] + paramT * temp[j - 1];
		}
	}
	return temp[degree];
}

std::vector<double> Polynomials::AllBernstein(int degree, double paramT)
{
	VALIDATE_ARGUMENT(degree > 0, "degree", "Degree must greater than zero.");
	VALIDATE_ARGUMENT_RANGE(paramT, 0.0, 1.0);

	std::vector<double> bernsteinArray(degree + 1);
	bernsteinArray[0] = 1.0;

	double t1 = 1.0 - paramT;
	for (int j = 1; j <= degree; j++)
	{
		double saved = 0.0;
		for (int k = 0; k < j; k++)
		{
			double temp = bernsteinArray[k];
			bernsteinArray[k] = saved + t1 * temp;
			saved = paramT * temp;
		}
		bernsteinArray[j] = saved;
	}
	return bernsteinArray;
}

double Polynomials::Horner(int degreeU, int degreeV, const std::vector<std::vector<double>>& coefficients, UV& uv)
{
	VALIDATE_ARGUMENT(degreeU > 0, "degreeU", "DegreeU must greater than zero.");
	VALIDATE_ARGUMENT(degreeV > 0, "degreeV", "DegreeV must greater than zero.");
	VALIDATE_ARGUMENT(coefficients.size() > 0, "coefficients", "Coefficients size must greater than zero.");
	VALIDATE_ARGUMENT(degreeU + 1 == static_cast<int>(coefficients.size()), "degreeU", "Coefficients row size equals degreeU plus one.");
	VALIDATE_ARGUMENT(degreeV + 1 == static_cast<int>(coefficients[0].size()), "degreeV", "Coefficients column size equals degreeV plus one.");
	VALIDATE_ARGUMENT_RANGE(uv.GetU(), 0.0, 1.0);
	VALIDATE_ARGUMENT_RANGE(uv.GetV(), 0.0, 1.0);

	std::vector<double> temp(degreeU + 1);
	for (int i = 0; i <= degreeU; i++)
	{
		temp[i] = Horner(degreeV, coefficients[i], uv.GetV());
	}
	return Horner(degreeU, temp, uv.GetU());
}

int LNLib::Polynomials::GetKnotMultiplicity(const std::vector<double>& knotVector, double knot)
{
	VALIDATE_ARGUMENT(knotVector.size() > 0, "knotVector", "KnotVector size must greater than zero.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidKnotVector(knotVector), "knotVector", "KnotVector must be a nondecreasing sequence of real numbers.");
	VALIDATE_ARGUMENT_RANGE(knot, knotVector[0], knotVector[knotVector.size() - 1]);
	
	int size = knotVector.size();
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

int LNLib::Polynomials::GetKnotSpanIndex(int degree, const std::vector<double>& knotVector, double paramT)
{
	VALIDATE_ARGUMENT(degree >= 0, "degree", "Degree must greater than or equals zero.");
	VALIDATE_ARGUMENT(knotVector.size() > 0, "knotVector", "KnotVector size must greater than zero.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidKnotVector(knotVector), "knotVector", "KnotVector must be a nondecreasing sequence of real numbers.");
	VALIDATE_ARGUMENT_RANGE(paramT, knotVector[0], knotVector[knotVector.size() - 1]);

	int n = knotVector.size() - degree - 2;
	if (MathUtils::IsGreaterThanOrEqual(paramT, knotVector[n + 1]))
	{
		return n;
	}
	if (MathUtils::IsLessThanOrEqual(paramT, knotVector[degree]))
	{
		return degree;
	}

	int low = 0;
	int high = n + 1;
	int mid = static_cast<int>(floor((low + high) / 2.0));

	while (paramT < knotVector[mid] || 
		   paramT >= knotVector[mid + 1])
	{
		if (paramT < knotVector[mid])
		{
			high = mid;
		}
		else
		{
			low = mid;
		}	
		mid = static_cast<int>(floor((low + high) / 2.0));
	}
	return mid;
}

std::vector<double> LNLib::Polynomials::BasisFunctions(int spanIndex, int degree, const std::vector<double>& knotVector, double paramT)
{
	VALIDATE_ARGUMENT(spanIndex >= 0, "spanIndex", "SpanIndex must greater than or equals zero.");
	VALIDATE_ARGUMENT(degree >= 0, "degree", "Degree must greater than or equals zero.");
	VALIDATE_ARGUMENT(knotVector.size() > 0, "knotVector", "KnotVector size must greater than zero.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidKnotVector(knotVector), "knotVector", "KnotVector must be a nondecreasing sequence of real numbers.");
	VALIDATE_ARGUMENT_RANGE(paramT, knotVector[0], knotVector[knotVector.size() - 1]);

	std::vector<double> basisFunctions(degree + 1);
	basisFunctions[0] = 1.0;

	std::vector<double> left(degree + 1);
	std::vector<double> right(degree + 1);

	for (int j = 1; j <= degree; j++)
	{
		left[j] = paramT - knotVector[spanIndex + 1 - j];
		right[j] = knotVector[spanIndex + j] - paramT;

		double saved = 0.0;

		for (int r = 0; r < j; r++)
		{
			double temp = basisFunctions[r] / (right[r + 1] + left[j - r]);
			basisFunctions[r] = saved + right[r + 1] * temp;
			saved = left[j - r] * temp;
		}
		basisFunctions[j] = saved;
	}
	return basisFunctions;
}

std::vector<std::vector<double>> LNLib::Polynomials::BasisFunctionsDerivatives(int spanIndex, int degree,  int derivative, const std::vector<double>& knotVector, double paramT)
{
	VALIDATE_ARGUMENT(spanIndex >= 0, "spanIndex", "SpanIndex must greater than or equals zero.");
	VALIDATE_ARGUMENT(degree > 0, "degree", "Degree must greater than zero.");
	VALIDATE_ARGUMENT(derivative <= degree, "derivative", "Derivative must not greater than degree.");
	VALIDATE_ARGUMENT(knotVector.size() > 0, "knotVector", "KnotVector size must greater than zero.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidKnotVector(knotVector), "knotVector", "KnotVector must be a nondecreasing sequence of real numbers.");
	VALIDATE_ARGUMENT_RANGE(paramT, knotVector[0], knotVector[knotVector.size() - 1]);

	std::vector<std::vector<double>> derivatives(derivative + 1, std::vector<double>(degree + 1));
	std::vector<std::vector<double>> ndu(degree + 1,std::vector<double>(degree + 1));

	ndu[0][0] = 1.0;

	std::vector<double> left(degree + 1);
	std::vector<double> right(degree + 1);

	double saved = 0.0;
	double temp = 0.0;

	for (int j = 1; j <= degree; j++)
	{
		left[j] = paramT - knotVector[spanIndex + 1 - j];
		right[j] = knotVector[spanIndex + j] - paramT;

		saved = 0.0;
		for (int r = 0; r < j; r++)
		{
			ndu[j][r] = right[r + 1] + left[j - r];
			temp = ndu[r][j - 1] / ndu[j][r];

			ndu[r][j] = saved + right[r + 1] * temp;
			saved = left[j - r] * temp;
		}
		ndu[j][j] = saved;
	}

	for (int j = 0; j <= degree; j++)
	{
		derivatives[0][j] = ndu[j][degree];
	}

	std::vector<std::vector<double>> a(2,std::vector<double>(degree + 1));
	for (int r = 0; r <= degree; r++)
	{
		int s1 = 0; 
		int s2 = 1;
		a[0][0] = 1.0;

		for (int k = 1; k <= derivative; k++)
		{
			double d = 0.0;
			int rk = r - k;
			int pk = degree - k;

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

	int r = degree;
	for (int k = 1; k <= derivative; k++)
	{
		for (int j = 0; j <= degree; j++)
		{
			derivatives[k][j] *= r;
		}
		r *= degree - k;
	}
	return derivatives;
}

double LNLib::Polynomials::OneBasisFunction(int spanIndex, int degree, const std::vector<double>& knotVector, double paramT)
{
	VALIDATE_ARGUMENT(spanIndex >= 0, "spanIndex", "SpanIndex must greater than or equals zero.");
	VALIDATE_ARGUMENT(degree > 0, "degree", "Degree must greater than zero.");
	VALIDATE_ARGUMENT(knotVector.size() > 0, "knotVector", "KnotVector size must greater than zero.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidKnotVector(knotVector), "knotVector", "KnotVector must be a nondecreasing sequence of real numbers.");
	VALIDATE_ARGUMENT_RANGE(paramT, knotVector[0], knotVector[knotVector.size() - 1]);

	int m = knotVector.size() - 1;
	if ((spanIndex == 0 && MathUtils::IsAlmostEqualTo(paramT, knotVector[0])) ||
		(spanIndex == m - degree - 1 && MathUtils::IsAlmostEqualTo(paramT, knotVector[m])))
	{
		return 1.0;
	}
	if (MathUtils::IsLessThan(paramT, knotVector[spanIndex]) || 
		MathUtils::IsGreaterThanOrEqual(paramT, knotVector[spanIndex + degree + 1]))
	{
		return 0.0;
	}

	std::vector<double> N(degree + spanIndex + 1, 0.0);
	for (int j = 0; j <= degree; j++)
	{
		if (MathUtils::IsGreaterThanOrEqual(paramT, knotVector[spanIndex + j]) &&
			MathUtils::IsLessThan(paramT, knotVector[spanIndex + j + 1]))
		{
			N[j] = 1.0;
		}
	}
	
	for (int k = 1; k <= degree; k++)
	{
		double saved = 0.0;
		if (!MathUtils::IsAlmostEqualTo(N[0], 0.0))
		{
			saved = ((paramT - knotVector[spanIndex]) * N[0]) / (knotVector[spanIndex + k] - knotVector[spanIndex]);
		}
						
		for (int j = 0; j < degree - k + 1; j++)
		{
			double knotLeft = knotVector[spanIndex + j + 1];
			double knotRight = knotVector[spanIndex + j + k + 1];
			if (MathUtils::IsAlmostEqualTo(N[j + 1], 0.0))
			{
				N[j] = saved;
				saved = 0.0;
			}
			else
			{
				double temp = N[j + 1] / (knotRight - knotLeft);
				N[j] = saved + (knotRight - paramT) * temp;
				saved = (paramT - knotLeft) * temp;
			}
		}
	}
	return N[0];
}

std::vector<double> LNLib::Polynomials::OneBasisFunctionDerivative(int spanIndex, int degree, int derivative, const std::vector<double>& knotVector, double paramT)
{
	VALIDATE_ARGUMENT(spanIndex >= 0, "spanIndex", "SpanIndex must greater than or equals zero.");
	VALIDATE_ARGUMENT(degree > 0, "degree", "Degree must greater than zero.");
	VALIDATE_ARGUMENT(derivative <= degree, "derivative", "Derivative must not greater than degree.");
	VALIDATE_ARGUMENT(knotVector.size() > 0, "knotVector", "KnotVector size must greater than zero.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidKnotVector(knotVector), "knotVector", "KnotVector must be a nondecreasing sequence of real numbers.");
	VALIDATE_ARGUMENT_RANGE(paramT, knotVector[0], knotVector[knotVector.size() - 1]);

	std::vector<double> derivatives(derivative + 1, 0.0);

	if (MathUtils::IsLessThan(paramT, knotVector[spanIndex]) ||
		MathUtils::IsGreaterThanOrEqual(paramT, knotVector[spanIndex + degree + 1]))
	{
		return derivatives;
	}

	std::vector<std::vector<double>> N(degree + 1, std::vector<double>(degree + 1,0.0));

	for (int j = 0; j <= degree; j++)
	{
		if (MathUtils::IsGreaterThanOrEqual(paramT, knotVector[spanIndex + j]) &&
			MathUtils::IsLessThan(paramT, knotVector[spanIndex + j + 1]))
		{
			N[j][0] = 1.0;
		}
	}
	for (int k = 1; k <= degree; k++)
	{
		double saved = 0.0;
		if (!MathUtils::IsAlmostEqualTo(N[0][k - 1], 0.0))
		{
			saved = ((paramT - knotVector[spanIndex]) * N[0][k - 1]) / (knotVector[spanIndex + k] - knotVector[spanIndex]);
		}
		for (int j = 0; j < degree - k + 1; j++)
		{
			double knotLeft = knotVector[spanIndex + j + 1];
			double knotRight = knotVector[spanIndex + j + k + 1];

			if (MathUtils::IsAlmostEqualTo(N[j + 1][k - 1], 0.0))
			{
				N[j][k] = saved;
				saved = 0.0;
			}
			else
			{
				double temp = N[j + 1][k - 1] / (knotRight - knotLeft);
				N[j][k] = saved + (knotRight - paramT) * temp;
				saved = (paramT - knotLeft) * temp;
			}
		}
	}

	derivatives[0] = N[0][degree];

	for (int k = 1; k <= derivative; k++)
	{
		std::vector<double> ND(k+1,0.0);
		for (int j = 0; j <= k; j++)
		{
			ND[j] = N[j][degree - k];
		}
		for (int jj = 1; jj <= k; jj++)
		{
			double saved = MathUtils::IsAlmostEqualTo(ND[0], 0.0)?
							0.0: ND[0] / (knotVector[spanIndex + degree - k + jj] - knotVector[spanIndex]);
			for (int j = 0; j < k - jj + 1; j++)
			{
				double knotLeft = knotVector[spanIndex + j + 1];
				double knotRight = knotVector[spanIndex + j + degree - k + jj + 1];

				if (MathUtils::IsAlmostEqualTo(ND[j + 1], 0.0))
				{
					ND[j] = (degree - k + jj) * saved;
					saved = 0.0;
				}
				else
				{
					double temp = ND[j + 1] / (knotRight - knotLeft);
					ND[j] = (degree - k + jj) * (saved - temp);
					saved = temp;
				}
			}
		}
		derivatives[k] = ND[0];
	}
	return derivatives;
}

std::vector<std::vector<double>> LNLib::Polynomials::AllBasisFunctions(int spanIndex, int degree, const std::vector<double>& knotVector, double knot)
{
	VALIDATE_ARGUMENT(spanIndex >= 0, "spanIndex", "SpanIndex must greater than or equals zero.");
	VALIDATE_ARGUMENT(degree >= 0, "degree", "Degree must greater than or equals zero.");
	VALIDATE_ARGUMENT(knotVector.size() > 0, "knotVector", "KnotVector size must greater than zero.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidKnotVector(knotVector), "knotVector", "KnotVector must be a nondecreasing sequence of real numbers.");
	VALIDATE_ARGUMENT_RANGE(knot, knotVector[0], knotVector[knotVector.size() - 1]);

	std::vector<std::vector<double>> result(degree + 1, std::vector<double>(degree + 1));
	for (int i = 0; i <= degree; i++)
	{
		std::vector<double> basis = Polynomials::BasisFunctions(spanIndex, i, knotVector, knot);
		for (int j = 0; j <= i; j++)
		{
			result[j][i] = basis[j];
		}
	}
	return result;
}

void LNLib::Polynomials::BezierToPowerMatrix(int degree, std::vector<std::vector<double>>& matrix)
{
	matrix.resize(degree+1, std::vector<double>(degree + 1));
	for (int i = 0; i < degree; i++)
	{
		for (int j = i + 1; j <= degree; j++)
		{
			matrix[i][j] = 0.0;
		}
	}

	matrix[0][0] = matrix[degree][degree] = 1.0;
	matrix[degree][0] = degree % 2 == 0 ? -1.0 : 1.0;

	double sign = -1.0;
	for (int i = 1; i < degree; i++)
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

void LNLib::Polynomials::PowerToBezierMatrix(int degree, const std::vector<std::vector<double>>& matrix, std::vector<std::vector<double>>& inverseMatrix)
{
	inverseMatrix.resize(degree + 1, std::vector<double>(degree + 1));
	for (int i = 0; i < degree; i++)
	{
		for (int j = i + 1; j <= degree; j++)
		{
			inverseMatrix[i][j] = 0.0;
		}
	}

	for (int i = 0; i <= degree; i++)
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
