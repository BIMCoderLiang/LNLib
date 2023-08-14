/*
 * Author:
 * 2023/06/08 - Yuqing Liang (BIMCoder Liang)
 * bim.frankliang@foxmail.com
 * 微信公众号：BIMCoder梁老师
 *
 * Use of this source code is governed by a GPL-3.0 license that can be found in
 * the LICENSE file.
 */

#include "MathUtils.h"
#include <limits>


bool LNLib::MathUtils::IsAlmostEqualTo(double value1, double value2, double tolerance)
{
    if (IsNaN(value1) || IsNaN(value2))
        return false;
    double eps = (abs(value1) + abs(value2) + 10) * tolerance;
    double delta = value1 - value2;
    return (-eps < delta) && (eps > delta);
}

bool LNLib::MathUtils::IsGreaterThan(double value1, double value2, double tolerance)
{
    if (IsNaN(value1) || IsNaN(value2))
        return false;
    return value1 > value2 && !IsAlmostEqualTo(value1, value2, tolerance);
}

bool LNLib::MathUtils::IsGreaterThanOrEqual(double value1, double value2, double tolerance)
{
    if (IsNaN(value1) || IsNaN(value2))
        return false;
    return (value1 - value2 > tolerance) || IsAlmostEqualTo(value1, value2, tolerance);
}

bool LNLib::MathUtils::IsLessThan(double value1, double value2, double tolerance)
{
    if (IsNaN(value1) || IsNaN(value2))
        return false;
    return value1 < value2 && !IsAlmostEqualTo(value1, value2, tolerance);
}

bool LNLib::MathUtils::IsLessThanOrEqual(double value1, double value2, double tolerance)
{
    if (IsNaN(value1) || IsNaN(value2))
        return false;
    return value1 < value2 || IsAlmostEqualTo(value1, value2, tolerance);
}

bool LNLib::MathUtils::IsInfinite(double value)
{
    constexpr double maxValue = (std::numeric_limits<double>::max)();
    double minValue = -maxValue;
    return !(minValue <= value && value <= maxValue);
}

bool LNLib::MathUtils::IsNaN(double value)
{
    return value != value;
}

double LNLib::MathUtils::RadiansToAngle(double radians)
{
    return radians * 180.0 / Constants::Pi;
}

double LNLib::MathUtils::AngleToRadians(double angle)
{
    return angle * Constants::Pi / 180.0;
}

int LNLib::MathUtils::Factorial(unsigned int number)
{
    if (number == 0)
        return 1;
    else
        return number * Factorial(number - 1);
}

double LNLib::MathUtils::Binomial(unsigned int number, unsigned int i)
{
    return Factorial(number) / (Factorial(i) * Factorial(number - 1));
}

double LNLib::MathUtils::ComputerCubicEquationsWithOneVariable(double cubic, double quadratic, double linear, double constant)
{
    double result;
    double initial = 0.001;
    result = initial - ((cubic * pow(initial,3) + quadratic * pow(initial, 2) + linear * initial + constant) / (3 * cubic * pow(initial, 2) + 2 * quadratic * initial + linear));
    while(MathUtils::IsGreaterThan(abs(result-initial),Constants::DoubleEpsilon))
    {
        initial = result;
        result = initial - ((cubic * pow(initial,3) + quadratic * pow(initial, 2) + linear * initial + constant) / (3 * cubic * pow(initial, 2) + 2 * quadratic * initial + linear));
    }
    return result;
}

std::vector<std::vector<double>> LNLib::MathUtils::MatrixMultiply(const std::vector<std::vector<double>>& matrix0, const std::vector<std::vector<double>>& matrix1)
{
    int m = static_cast<int>(matrix0.size());
    int n = static_cast<int>(matrix0[0].size());
    int p = static_cast<int>(matrix1[0].size());
    std::vector<std::vector<double>> result;
    std::vector<double> temparray;
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < p; j++)
        {
            double sum = 0.0;
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

std::vector<std::vector<double>> LNLib::MathUtils::MakeDiagonal(int size)
{
    std::vector<std::vector<double>> result(size);
    for (int i = 0; i < size; i++)
    {
        result[i].resize(size);
        for (int j = 0; j < size; j++)
        {
            if (i == j)
            {
                result[i][j] = 1;
            }
            else
            {
                result[i][j] = 0;
            }
        }
    }
    return result;
}



std::vector<std::vector<double>> LNLib::MathUtils::CreateMatrix(int row, int column)
{
    std::vector<std::vector<double>> result;
    for (int i = 0; i < row; i++)
    {
        std::vector<double> v(column, 0);
        result.emplace_back(v);
    }
    return result;
}

bool LNLib::MathUtils::MakeInverse(const std::vector<std::vector<double>>& matrix, std::vector<std::vector<double>>& inverse)
{
    int rows = matrix.size();
    int columns = matrix[0].size();
    if (rows != columns) return false;

    int n = rows;
    inverse = CreateMatrix(n, n);

    std::vector<std::vector<double>> lower;
    std::vector<std::vector<double>> upper;
    if (!LUDecomposition(matrix, lower, upper))
    {
        return false;
    }

    std::vector<std::vector<double>> inverseLower(n);
    std::vector<std::vector<double>> inverseUpper(n);
    for (int i = 0; i < n; i++)
    {
        inverseLower[i].resize(n);
        inverseUpper[i].resize(n);
    }

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

    for (int i = 0; i < n; i++)           
    {
        for (int j = 0; j < n; j++)
        {
            for (int k = 0; k < n; k++)
            {
                inverse[i][j] += inverseUpper[i][k] * inverseLower[k][j];
            }
        }
    }
    return true;
}

bool LNLib::MathUtils::LUDecomposition(const std::vector<std::vector<double>>& matrix, std::vector<std::vector<double>>& lowerTriMatrix, std::vector<std::vector<double>>& upperTriMatrix)
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
            double temp = 0.0;
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

std::vector<double> LNLib::MathUtils::ForwardSubstitution(const std::vector<std::vector<double>>& lowerTriMatrix, const std::vector<double>& column)
{
    int size = static_cast<int>(column.size());

    std::vector<double> result;
    result.resize(size, 0.0);

    result[0] = column[0] / lowerTriMatrix[0][0];
    for (int i = 1; i < size; i++)
    {
        double temp = 0.0;
        for (int j = 0; j < i; j++)
        {
            temp += lowerTriMatrix[i][j] * result[j];
        }
        result[i] = (column[i] - temp) / lowerTriMatrix[i][i];
    }
    return result;
}

std::vector<double> LNLib::MathUtils::BackwardSubstitution(const std::vector<std::vector<double>>& upperTriMatrix, const std::vector<double>& column)
{
    int size = static_cast<int>(column.size());
    int n = size - 1;
    std::vector<double> result;
    result.resize(size, 0.0);

    result[n] = column[n] / upperTriMatrix[n][n];
    for (int i = size - 2; i >= 0; i--)
    {
        double temp = 0.0;
        for (int j = i; j < size; j++)
        {
            temp += upperTriMatrix[i][j] * result[j];
        }
        result[i] = (column[i] - temp) / upperTriMatrix[i][i];
    }
    return result;
}
