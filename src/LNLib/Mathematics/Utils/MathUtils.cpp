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

std::vector<std::vector<double>> LNLib::MathUtils::MatrixMultiply(const std::vector<std::vector<double>>& left, const std::vector<std::vector<double>>& right)
{
    int m = static_cast<int>(left.size());
    int n = static_cast<int>(left[0].size());
    int p = static_cast<int>(right[0].size());

    std::vector<std::vector<double>> result(m, std::vector<double>(p,0.0));
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < p; j++)
        {
            for (int k = 0; k < n; k++)
            {
                result[i][j] += left[i][k] * right[k][j];
            }
        }
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

bool LNLib::MathUtils::IsSquareMatrix(const std::vector<std::vector<double>>& matrix)
{
    int row = static_cast<int>(matrix.size());
    int column = static_cast<int>(matrix[0].size());
    return row == column;
}

double LNLib::MathUtils::GetDeterminant(const std::vector<std::vector<double>>& matrix)
{
    if (!IsSquareMatrix(matrix))
    {
        return 0.0;
    }
       
    int n = static_cast<int>(matrix.size());
    std::vector<std::vector<double>> temp = matrix;
    double det = 1.0;
    for (int i = 0; i < n; i++)
    {
        int pivot = i;
        for (int j = i + 1; j < n; j++)
        {
            if (abs(temp[j][i]) > abs(temp[pivot][i]))
            {
                pivot = j;
            }
        }
        if (pivot != i)
        {
            std::swap(temp[i], temp[pivot]);
            det *= -1;
        }
        if (temp[i][i] == 0)
        {
            return 0.0;
        }
        det *= temp[i][i];
        for (int j = i + 1; j < n; j++)
        {
            double factor = temp[j][i] / temp[i][i];
            for (int k = i + 1; k < n; k++)
            {
                temp[j][k] -= factor * temp[i][k];
            }
        }
    }
    return det;
}

bool LNLib::MathUtils::MakeInverse(const std::vector<std::vector<double>>& matrix, std::vector<std::vector<double>>& inverse)
{
    if (!IsSquareMatrix(matrix))
    {
        return false;
    }
    double det = GetDeterminant(matrix);
    if (IsAlmostEqualTo(det, 0.0))
    {
        return false;
    }

    int n = static_cast<int>(matrix.size());
    std::vector<std::vector<double>> lower;
    std::vector<std::vector<double>> upper;
    if (LUDecomposition(matrix, lower, upper))
    {
        std::vector<std::vector<double>> inverseLower(n, std::vector<double>(n));
        std::vector<std::vector<double>> inverseUpper(n, std::vector<double>(n));

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
                if (IsAlmostEqualTo(abs(s), 0.0))
                {
                    inverseUpper[k][i] = 0.0;
                }
                else
                {
                    inverseUpper[k][i] = -s / upper[k][k];
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
                    double temp = inverseLower[k][i] - lower[k][j] * inverseLower[j][i];
                    if (IsAlmostEqualTo(temp, 0.0))
                    {
                        inverseLower[k][i] = 0.0;
                    }
                    else
                    {
                        inverseLower[k][i] = temp;
                    }
                }
            }
        }

        inverse = MatrixMultiply(inverseUpper, inverseLower);
        return true;
    }
    else
    {
        bool rs = true;
        std::vector<std::vector<double>> tempInverse;
        for (int i = 0; i < n; i++)
        {
            std::vector<double> b(n, 0.0);
            b[i] = 1;
            std::vector<double> pivot;
            if (!LUPDecomposition(matrix, lower, upper, pivot))
            {
                rs = false;
                break;
            }
            std::vector<double> x(n);
            std::vector<double> y(n);

            for (int i = 0; i < n; i++)
            {
                y[i] = b[pivot[i]];
                for (int j = 0; j < i; j++)
                {
                    y[i] = y[i] - lower[i][j] * y[j];
                }
            }
            for (int i = n - 1; i >= 0; i--)
            {
                x[i] = y[i];
                for (int j = n - 1; j > i; j--)
                {
                    x[i] = x[i] - upper[i][j] * x[j];
                }
                x[i] /= upper[i][i];
            }
            tempInverse.emplace_back(x);
        }
        if (rs)
        {
            Transpose(tempInverse, inverse);
            return true;
        }
        else
        {
            //to do...
            return true;
        }
    }
}

bool LNLib::MathUtils::LUDecomposition(const std::vector<std::vector<double>>& matrix, std::vector<std::vector<double>>& lowerTriMatrix, std::vector<std::vector<double>>& upperTriMatrix)
{
    if (!IsSquareMatrix(matrix))
    {
        return false;
    }

    int n = static_cast<int>(matrix.size());
    lowerTriMatrix.resize(n);
    upperTriMatrix.resize(n);
    for (int i = 0; i < n; i++)
    {
        lowerTriMatrix[i].resize(n);
        upperTriMatrix[i].resize(n);
        for (int j = 0; j < n; j++)
        {
            upperTriMatrix[i][j] = 0;
            if (i == j)
            {
                lowerTriMatrix[i][j] = 1.0;
            }
        }
    }
    
    for (int i = 0; i < n; i++)
    {
        double sum = 0;
        for (int j = i; j < n; j++)
        {
            for (int k = 0; k <= i - 1; k++)
            {
                sum += lowerTriMatrix[i][k] * upperTriMatrix[k][j];
            }
            double temp = matrix[i][j] - sum;
            if (IsAlmostEqualTo(temp, 0.0))
            {
                if (i == j) return false;
                upperTriMatrix[i][j] = 0.0;
            }
            else
            {
                upperTriMatrix[i][j] = temp;
            }
            sum = 0.0;
        }

        for (int x = i + 1; x < n; x++)
        {
            for (int k = 0; k <= i - 1; k++)
            {
                sum += lowerTriMatrix[x][k] * upperTriMatrix[k][i];
            }
            double temp = matrix[x][i] - sum;
            if (IsAlmostEqualTo(temp, 0.0))
            {
                lowerTriMatrix[x][i] = 0.0;
            }
            else
            {
                lowerTriMatrix[x][i] = temp / upperTriMatrix[i][i];
            }
            sum = 0.0;
        }
    }
    return true;
}

bool LNLib::MathUtils::LUPDecomposition(const std::vector<std::vector<double>>& matrix, std::vector<std::vector<double>>& lowerTriMatrix, std::vector<std::vector<double>>& upperTriMatrix, std::vector<double>& pivot)
{
    if (!IsSquareMatrix(matrix))
    {
        return false;
    }

    int n = static_cast<int>(matrix.size());
    std::vector<std::vector<double>> copy(n, std::vector<double>(n));
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            copy[i][j] = matrix[i][j];
            lowerTriMatrix[i][j] = 0.0;
            upperTriMatrix[i][j] = 0.0;
        }
    }
    pivot.resize(n);
    for (int i = 0; i < n; i++)
    {
        pivot[i] = i;
    }

    int row = 0;
    for (int i = 0; i < n - 1; i++)
    {
        double p = 0.0;
        for (int j = i; j < n; j++)
        {
            if (IsGreaterThan(abs(copy[j][i]),p))
            {
                p = abs(copy[j][i]);
                row = j;
            }
        }
        if (IsAlmostEqualTo(p,0.0))
        {
            return false;
        }

        int tmp = pivot[i];
        pivot[i] = pivot[row];
        pivot[row] = tmp;

        double tmp2 = 0.0;
        for (int j = 0; j < n; j++)
        {
            tmp2 = copy[i][j];
            copy[i][j] = copy[row][j];
            copy[row][j] = tmp2;
        }

        double u = copy[i][i];
        double l = 0.0;
        for (int j = i + 1; j < n; j++)
        {
            l = copy[j][i] / u;
            copy[j][i] = l;
            for (int k = i + 1; k < n; k++)
            {
                copy[j][k] = copy[j][k] - copy[i][k] * l;
            }
        }
    }

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j <= i; j++)
        {
            if (i != j)
            {
                lowerTriMatrix[i][j] = copy[i][j];
            }
            else
            {
                lowerTriMatrix[i][j] = 1;
            }
        }
        for (int k = i; k < n; k++)
        {
            upperTriMatrix[i][k] = copy[i][k];
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
