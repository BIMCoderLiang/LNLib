/*
 * Author:
 * 2023/06/08 - Yuqing Liang (BIMCoder Liang)
 * bim.frankliang@foxmail.com
 * 微信公众号：BIMCoder梁老师
 *
 * Use of this source code is governed by a GPL-3.0 license that can be found in
 * the LICENSE file.
 */

#pragma once

#include "LNLibDefinitions.h"
#include "Constants.h"
#include <vector>

namespace LNLib
{
	class LNLIB_EXPORT MathUtils
	{

	public:

		static bool IsAlmostEqualTo(double value1, double value2, double tolerance = Constants::DoubleEpsilon);

		static bool IsGreaterThan(double value1, double value2, double tolerance = Constants::DoubleEpsilon);

		static bool IsGreaterThanOrEqual(double value1, double value2, double tolerance = Constants::DoubleEpsilon);

		static bool IsLessThan(double value1, double value2, double tolerance = Constants::DoubleEpsilon);

		static bool IsLessThanOrEqual(double value1, double value2, double tolerance = Constants::DoubleEpsilon);

		static bool IsInfinite(double value);

		static bool IsNaN(double value);

		static double RadiansToAngle(double radians);

		static double AngleToRadians(double angle);

		static int Factorial(unsigned int number);

		static double Binomial(unsigned int number, unsigned int i);

		template<typename T>
		static void Transpose(const std::vector<std::vector<T>>& matrix, std::vector<std::vector<T>> transposed)
		{
			std::vector<T> temp;

			for (int i = 0; i < static_cast<int>(matrix[0].size()); i++)
			{
				for (int j = 0; j < static_cast<int>(matrix.size()); j++)
				{
					temp.emplace_back(matrix[j][i]);
				}
				transposed.emplace_back(temp);
				temp.erase(temp.begin(), temp.end());
			}
		}

		template<typename T>
		static void GetColumn(const std::vector<std::vector<T>>& matrix, unsigned int columnIndex, std::vector<T> columnData)
		{
			int size = static_cast<int>(matrix.size());
			columnData.resize(size);
			for (int i = 0; i < size; i++)
			{
				columnData.emplace_back(matrix[i][columnIndex]);
			}
		}

		template<typename T>
		static std::vector<std::vector<T>> MatrixMultiply(const std::vector<std::vector<T>>& matrix0, const std::vector<std::vector<T>>& matrix1)
		{
			int m = static_cast<int>(matrix0.size());
			int n = static_cast<int>(matrix0[0].size());
			int p = static_cast<int>(matrix1[0].size());
			std::vector<std::vector<T>> result;
			std::vector<T> temparray;
			for (int i = 0; i < m; i++)
			{
				for (int j = 0; j < p; j++)
				{
					T sum = 0;
					for (int k = 0; k < n; k++)
					{
						sum += matrix0[i][k] + matrix1[k][j];
					}
					temparray.emplace_back(sum);
				}
				result.emplace_back(temparray);
				temparray.erase(temparray.begin(), temparray.end());
			}
			return result;
		}

		template<typename T>
		static std::vector<std::vector<T>> MakeDiagonal(const std::vector<T>& data)
		{
			int size = static_cast<int>(data.size());
			std::vector<std::vector<T>> result(size);
			for (int i = 0; i < size; i++)
			{
				result[i].resize(size);
				for (int j = 0; j < size; j++)
				{
					if (i == j)
					{
						result[i][j] = data[i];
					}
					else
					{
						result[i][j] = 0;
					}
				}
			}
			return result;
		}

		static std::vector<std::vector<double>> CreateMatrix(int row, int column);

		template<typename T>
		static bool MakeInverse(const std::vector<std::vector<T>>& matrix, std::vector<std::vector<T>> inverse)
		{
			int rows = matrix.size();
			int columns = matrix[0].size();
			if (rows != columns) return false;

			int n = rows;
			inverse = CreateMatrix(n, n);

			std::vector<std::vector<T>> lower;
			std::vector<std::vector<T>> upper;
			if (!LUDecomposition(matrix, lower, upper))
			{
				return false;
			}

			std::vector<std::vector<T>> inverseLower;
			std::vector<std::vector<T>> inverseUpper;

			for (int i = 0; i < n; i++)
			{
				inverseUpper[i][i] = 1 / upper[i][i];
				for (int k = i - 1; k >= 0; k--)
				{
					double s = 0;
					for (int j = k + 1; j <= i; j++)
					{
						s = s + upper[k][j] * inverseUpper[j][i];
					}
					inverseUpper[k][i] = -s / upper[k][k];
					if (IsLessThan(std::abs(inverseUpper[k][i]), Constants::DoubleEpsilon))
					{
						inverseUpper[k][i] = 0.0;
					}
				}
			}

			for (int i = 0; i < n; i++)
			{
				inverseLower[i][i] = 1;
				for (int k = i + 1; k < n; k++)
				{
					for (int j = i; j <= k - 1; j++)
					{
						inverseLower[k][i] = inverseLower[k][i] - lower[k][j] * inverseLower[j][i];
						if (IsLessThan(std::abs(inverseLower[k][i]), Constants::DoubleEpsilon))
						{
							inverseLower[k][i] = 0.0;
						}
					}
				}
			}

			inverse = MatrixMultiply(inverseUpper, inverseLower);
			return true;
		}

		template<typename T>
		static bool LUDecomposition(const std::vector<std::vector<T>>& matrix, std::vector<std::vector<T>>& lowerTriMatrix, std::vector<std::vector<T>>& upperTriMatrix)
		{
			int row = static_cast<int>(matrix.size());
			if (row <= 0) return false;
			int column = static_cast<int>(matrix[0].size());
			if (row != column) return false;

			lowerTriMatrix.resize(row);
			for (int i = 0; i < row; i++)
			{
				lowerTriMatrix[i].resize(row, 0.0);
			}

			upperTriMatrix.resize(row);
			for (int i = 0; i < row; i++)
			{
				upperTriMatrix[i].resize(row, 0.0);
			}

			for (int i = 0; i < row; i++)
			{
				for (int k = 0; k < row; k++)
				{
					T temp = 0.0;
					for (int j = 0; j < i; j++)
					{
						temp += lowerTriMatrix[i][j] * upperTriMatrix[j][k];
					}
					upperTriMatrix[i][k] = matrix[i][k] - temp;

					if (i == k)
					{
						lowerTriMatrix[i][i] = 1.0;
					}
					else
					{
						temp = 0.0;
						for (int j = 0; j < i; j++)
						{
							temp += lowerTriMatrix[k][j] * upperTriMatrix[j][i];
						}
						lowerTriMatrix[k][i] = matrix[k][i] - temp;

						if (IsAlmostEqualTo(upperTriMatrix[i][i], 0.0))
						{
							lowerTriMatrix[k][i] = 0.0;
						}
						else
						{
							lowerTriMatrix[k][i] /= upperTriMatrix[i][i];
						}
					}
				}
			}
			return true;
		}

		template<typename T>
		static std::vector<T> ForwardSubstitution(const std::vector<std::vector<T>>& lowerTriMatrix, const std::vector<T>& column)
		{
			int size = static_cast<int>(column.size());

			std::vector<T> result;
			result.resize(size, 0.0);

			result[0] = column[0] / lowerTriMatrix[0][0];
			for (int i = 1; i < size; i++)
			{
				T temp = 0.0;
				for (int j = 0; j < i; j++)
				{
					temp += lowerTriMatrix[i][j] * result[j];
				}
				result[i] = (column[i] - temp) / lowerTriMatrix[i][i];
			}
			return result;
		}

		template<typename T>
		static std::vector<T> BackwardSubstitution(const std::vector<std::vector<T>>& upperTriMatrix, const std::vector<T>& column)
		{
			int size = static_cast<int>(column.size());
			int n = size - 1;
			std::vector<T> result;
			result.resize(size, 0.0);

			result[n] = column[n] / upperTriMatrix[n][n];
			for (int i = size - 2; i >= 0; i--)
			{
				T temp = 0.0;
				for (int j = i; j < size; j++)
				{
					temp += upperTriMatrix[i][j] * result[j];
				}
				result[i] = (column[i] - temp) / upperTriMatrix[i][i];
			}
			return result;
		}
	};
}


