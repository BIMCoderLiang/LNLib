/*
 * Author:
 * 2023/06/08 - Yuqing Liang (BIMCoder Liang)
 * bim.frankliang@foxmail.com
 * 微信公众号：BIMCoder梁老师
 *
 * Use of this source code is governed by a GPL-3.0 license that can be found in
 * the LICENSE file.
 */


/*
 * 
 * NoTICE: The Eigen C++ Library is Only used in MathUitls.cpp 
 * since we should provide uniform LNLib API Style to the public.
 * 
 */

#include "MathUtils.h"

#include <Eigen/Dense>

#include <cmath>
#include <limits>

bool LNLib::MathUtils::IsAlmostEqualTo(double value1, double value2, double tolerance)
{
    if (IsNaN(value1) || IsNaN(value2))
        return false;
    double eps = (std::abs(value1) + std::abs(value2) + 10) * tolerance;
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

int LNLib::MathUtils::Factorial(int number)
{
    if (number == 0)
        return 1;
    else
        return number * Factorial(number - 1);
}

double LNLib::MathUtils::Binomial(int number, int i)
{
    return Factorial(number) / (Factorial(i) * Factorial(number - i));
}

double LNLib::MathUtils::ComputerCubicEquationsWithOneVariable(double cubic, double quadratic, double linear, double constant)
{
    double result;
    double initial = 0.001;
    result = initial - ((cubic * std::pow(initial,3) + quadratic * std::pow(initial, 2) + linear * initial + constant) / (3 * cubic * std::pow(initial, 2) + 2 * quadratic * initial + linear));
    while(MathUtils::IsGreaterThan(std::abs(result-initial),Constants::DoubleEpsilon))
    {
        initial = result;
        result = initial - ((cubic * std::pow(initial,3) + quadratic * std::pow(initial, 2) + linear * initial + constant) / (3 * cubic * std::pow(initial, 2) + 2 * quadratic * initial + linear));
    }
    return result;
}

std::vector<std::vector<double>> LNLib::MathUtils::MatrixMultiply(const std::vector<std::vector<double>>& left, const std::vector<std::vector<double>>& right)
{
    int m = left.size();
    int n = left[0].size();
    int p = right[0].size();

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

double LNLib::MathUtils::GetDeterminant(const std::vector<std::vector<double>>& matrix)
{
    Eigen::MatrixXd m(matrix.size(), matrix[0].size());
    for (int i = 0; i < matrix.size(); i++)
    {
        m.row(i) = Eigen::VectorXd::Map(&matrix[i][0], matrix[i].size());
    }
    return m.determinant();
}

bool LNLib::MathUtils::MakeInverse(const std::vector<std::vector<double>>& matrix, std::vector<std::vector<double>>& inverse)
{
    int row = matrix.size();
    int column = matrix[0].size();

    if (!(row == column))
    {
        return false;
    }
    Eigen::MatrixXd m(matrix.size(), matrix[0].size());
    for (int i = 0; i < row; i++)
    {
        m.row(i) = Eigen::VectorXd::Map(&matrix[i][0], matrix[i].size());
    }

    auto result = m.inverse();
    inverse.resize(row);
    for (int k = 0; k < row; k++)
    {
        inverse[k].resize(column);
    }

    for (int col = 0; col < result.cols(); ++col)
    {
        for (int row = 0; row < result.rows(); ++row)
        {
            inverse[row][col] = result(row, col);
        }
    }
    return true;
}

std::vector<std::vector<double>> LNLib::MathUtils::SolveLinearSystem(const std::vector<std::vector<double>>& matrix, const std::vector<std::vector<double>>& right)
{
    std::vector<std::vector<double>> result(matrix.size(), std::vector<double>(right[0].size()));

    Eigen::MatrixXd m(matrix.size(), matrix[0].size());
    for (int i = 0; i < matrix.size(); i++)
    {
        m.row(i) = Eigen::VectorXd::Map(&matrix[i][0], matrix[i].size());
    }
        
    Eigen::MatrixXd r(right.size(),  right[0].size());
    for (int i = 0; i < right.size(); i++)
    {
        r.row(i) = Eigen::VectorXd::Map(&right[i][0], right[i].size());
    }

    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> solve = m.lu().solve(r);
    for (int col = 0; col < solve.cols(); ++col)
    {
        for (int row = 0; row < solve.rows(); ++row)
        {
            result[row][col] = solve(row, col);
        }
    }
    return result;
}


