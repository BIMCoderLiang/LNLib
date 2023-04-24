
#include "Polynomials.h"
#include "ValidationUtils.h"
#include "MathUtils.h"

using namespace LNLib;

double Polynomials::Horner(const std::vector<double>& coefficients, unsigned int degree, double paramT)
{
	double result = coefficients[degree];
	for (int i = degree - 1; i >= 0; i--)
	{
		result = result * paramT + coefficients[i];
	}
	return result;
}

double Polynomials::Bernstein(unsigned int i, unsigned int degree, double paramT)
{
	std::vector<double> temp;
	temp.resize(degree + 1, 0.0);

	temp[degree - i] = 1.0;
	double t1 = 1.0 - paramT;

	for (unsigned int k = i; k <= degree; k++)
	{
		for (unsigned int j = degree; j >= k; j--)
		{
			temp[j] = t1 * temp[j] + paramT * temp[j - 1];
		}
	}
	return temp[degree];
}

void Polynomials::AllBernstein(unsigned int degree, double paramT, std::vector<double>& bernsteinArray)
{
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
}

double Polynomials::Horner(const std::vector<std::vector<double>>& coefficients, unsigned int degreeU, unsigned int degreeV, double paramU, double paramV)
{
	std::vector<double> temp;
	temp.resize(degreeU + 1);

	for (unsigned int i = 0; i <= degreeU; i++)
	{
		temp[i] = Horner(coefficients[i], degreeV, paramV);
	}
	return Horner(temp, degreeU, paramU);
}

int LNLib::Polynomials::GetKnotSpanIndex(unsigned int n, unsigned int degree, double knot, const std::vector<double>& knotVector)
{
	int m = static_cast<int>(knotVector.size()) - degree - 2;
	if (m >= 0)
	{
		if (knot == knotVector[n + 1]) return n;

		int low = degree;
		int high = n + 1;
		int mid = static_cast<int>(floor((low + high) / 2.0));

		while (knot < knotVector[mid] ||
				knot >= knotVector[mid + 1])
		{
			if (knot < knotVector[mid])
				high = mid;
			else
				low = mid;
			mid = static_cast<int>(floor((low + high) / 2.0));
		}
		return mid;
	}
	return -1;
}

void LNLib::Polynomials::BasisFunctions(unsigned int spanIndex, unsigned int degree, double paramT, const std::vector<double>& knotVector, std::vector<double>& basisFunctions)
{
	basisFunctions[0] = 1.0;
	for (unsigned int j = 1; j <= degree; j++)
	{
		std::vector<double> left;
		left.resize(degree+1);
		std::vector<double> right;
		right.resize(degree+1);

		left[j] = paramT - knotVector[spanIndex + 1 - j];
		right[j] = knotVector[spanIndex + j] - paramT;

		double saved = 0.0;

		for (unsigned int r = 0; r < j; r++)
		{
			double temp = basisFunctions[r] / (right[r + 1] + left[j - r]);
			basisFunctions[r] = saved + right[r + 1] * temp;
			saved = left[j - r] * temp;
		}
		basisFunctions[j] = saved;
	}
}

void LNLib::Polynomials::BasisFunctionsDerivatives(unsigned int spanIndex, unsigned int degree, double paramT,  unsigned int derivative, const std::vector<double>& knotVector, std::vector<std::vector<double>>& derivatives)
{
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

	for (unsigned int j = 1; j <= degree; j++)
	{
		double saved = 0.0;

		left[j] = paramT - knotVector[spanIndex + 1 - j];
		right[j] = knotVector[spanIndex + j] - paramT;

		for (unsigned int r = 0; r < j; r++)
		{
			ndu[j][r] = right[r + 1] + left[j - r];
			double temp = ndu[r][j - 1] / ndu[j][r];

			ndu[r][j] = saved + right[r + 1] * temp;
			saved = left[j - r] * temp;
		}
		ndu[j][j] = saved;
	}

	for (unsigned int j = 0; j <= degree; j++)
	{
		derivatives[0][j] = ndu[j][degree];
	}


	std::vector<std::vector<double>> a;
	a.resize(2);
	for (int i = 0; i < 2; i++)
	{
		a[i].resize(degree + 1);
	}

	for (unsigned int r = 0; r <= degree; r++)
	{
		int s1 = 0; 
		int s2 = 1;
		a[0][0] = 1.0;

		for (unsigned int k = 1; k <= derivative; k++)
		{
			double d = 0.0;
			int rk = r - k;
			unsigned int pk = degree - k;

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
	for (unsigned int k = 1; k <= derivative; k++)
	{
		for (unsigned int j = 0; j <= degree; j++)
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

	for (unsigned int j = 0; j <= degree; j++)
	{
		temp[j] = MathUtils::IsGreaterThanOrEqual(paramT, knotVector[index + j]) && 
					MathUtils::IsLessThan(paramT, knotVector[index + j + 1]) ? 
						1.0 : 0.0;
	}
	for (unsigned int k = 1; k <= degree; k++)
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
		for (unsigned int j = 0; j < degree - k + 1; j++)
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
	std::vector<std::vector<double>> basisFunctions;
	basisFunctions.resize(degree + 1);
	for (unsigned int i = 0; i <= degree; i++)
	{
		basisFunctions[i].resize(degree + 1);
	}


	if (MathUtils::IsLessThan(paramT, knotVector[index]) ||
		MathUtils::IsGreaterThanOrEqual(paramT, knotVector[index + degree + 1]))
	{
		for (unsigned int k = 0; k <= derivative; k++)
		{
			derivatives[k] = 0.0;
		}
		return;
	}
	for (unsigned int j = 0; j <= degree; j++)
	{
		basisFunctions[j][0] = MathUtils::IsGreaterThanOrEqual(paramT, knotVector[index + j]) &&
								MathUtils::IsLessThan(paramT, knotVector[index + j + 1]) ?
								  1.0 : 0.0;
	}
	for (unsigned int k = 1; k <= degree; k++)
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
		for (unsigned int j = 0; j < degree - k + 1; j++)
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
	for (unsigned int k = 1; k <= derivative; k++)
	{
		std::vector<double> basisFunctionsDerivative;
		basisFunctionsDerivative.resize(derivative);

		for (unsigned int j = 0; j <= k; j++)
		{
			basisFunctionsDerivative[j] = basisFunctions[j][degree - k];
		}
		for (unsigned int jj = 1; jj <= k; jj++)
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

			for (unsigned int j = 0; j < k - jj + 1; j++)
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
