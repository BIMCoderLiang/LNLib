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

		static std::vector<std::vector<double>> MatrixMultiply(const std::vector<std::vector<double>>& matrix0, const std::vector<std::vector<double>>& matrix1);
	
		static std::vector<std::vector<double>> MakeDiagonal(const std::vector<double>& data);
		
		static std::vector<std::vector<double>> CreateMatrix(int row, int column);

		static bool MakeInverse(const std::vector<std::vector<double>>& matrix, std::vector<std::vector<double>> inverse);
		
		static bool LUDecomposition(const std::vector<std::vector<double>>& matrix, std::vector<std::vector<double>>& lowerTriMatrix, std::vector<std::vector<double>>& upperTriMatrix);
		
		static std::vector<double> ForwardSubstitution(const std::vector<std::vector<double>>& lowerTriMatrix, const std::vector<double>& column);

		static std::vector<double> BackwardSubstitution(const std::vector<std::vector<double>>& upperTriMatrix, const std::vector<double>& column);
	};
}


