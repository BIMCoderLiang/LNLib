/*
 * Author:
 * 2023/06/08 - Yuqing Liang (BIMCoder Liang)
 * bim.frankliang@foxmail.com
 * 
 *
 * Use of this source code is governed by a LGPL-2.1 license that can be found in
 * the LICENSE file.
 */


/*
 * 
 * NOTICE: The Eigen C++ Library is Only used in MathUitls.cpp 
 * since we should provide uniform style API to the public.
 * 
 */

#include "MathUtils.h"

#include <Eigen/Dense>

#include <cmath>
#include <limits>

namespace LNLib {

    static const std::vector<std::vector<double>> pascal = [] {
        std::vector<std::vector<double>> table(21, std::vector<double>(21, 0.0));
        for (int i = 0; i <= 20; ++i) {
            table[i][0] = table[i][i] = 1.0;
            for (int j = 1; j < i; ++j) {
                table[i][j] = table[i - 1][j - 1] + table[i - 1][j];
            }
        }
        return table;
        }();
}

bool LNLib::MathUtils::IsAlmostEqualTo(double a, double b, double tolerance)
{
    if (IsNaN(a) || IsNaN(b)) return false;
    if (IsInfinite(a) || IsInfinite(b)) return a == b;

    double diff = std::abs(a - b);
    double absMax = std::max(std::abs(a), std::abs(b));
    return diff <= tolerance * (1.0 + absMax);
}

bool LNLib::MathUtils::IsGreaterThan(double a, double b, double tolerance)
{
    if (IsNaN(a) || IsNaN(b)) return false;
    return (a - b) > tolerance * (1.0 + std::max(std::abs(a), std::abs(b)));
}

bool LNLib::MathUtils::IsGreaterThanOrEqual(double a, double b, double tolerance)
{
    if (IsNaN(a) || IsNaN(b)) return false;
    double eps = tolerance * (1.0 + std::max(std::abs(a), std::abs(b)));
    return (a - b) >= -eps;
}

bool LNLib::MathUtils::IsLessThan(double a, double b, double tolerance)
{
    if (IsNaN(a) || IsNaN(b)) return false;
    return (b - a) > tolerance * (1.0 + std::max(std::abs(a), std::abs(b)));
}

bool LNLib::MathUtils::IsLessThanOrEqual(double a, double b, double tolerance)
{
    if (IsNaN(a) || IsNaN(b)) return false;
    double eps = tolerance * (1.0 + std::max(std::abs(a), std::abs(b)));
    return (a - b) <= eps;
}

bool LNLib::MathUtils::IsInfinite(double value)
{
    return std::isinf(value);
}

bool LNLib::MathUtils::IsNaN(double value)
{
    return std::isnan(value);
}

double LNLib::MathUtils::RadiansToAngle(double radians)
{
    return radians * 180.0 / Constants::Pi;
}

double LNLib::MathUtils::AngleToRadians(double angle)
{
    return angle * Constants::Pi / 180.0;
}

long long LNLib::MathUtils::Factorial(int n)
{
    if (n < 0) return -1;
    if (n > 20) {
        throw std::invalid_argument("Factorial overflow for n > 20");
    }
    long long result = 1;
    for (int i = 2; i <= n; ++i) {
        result *= i;
    }
    return result;
}

double LNLib::MathUtils::Binomial(int number, int i)
{
    if (number < 0 || i < 0 || i > number) return 0.0;
    if (number <= 20) return pascal[number][i];

    i = std::min(i, number - i);
    double res = 1.0;
    for (int j = 1; j <= i; ++j)
        res = res * (number - i + j) / j;
    return res;
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

std::vector<std::vector<double>> LNLib::MathUtils::SolveLinearSystem(
    const std::vector<std::vector<double>>& matrix,
    const std::vector<std::vector<double>>& right)
{
    if (matrix.empty() || right.empty() || matrix[0].empty() || right[0].empty()) {
        return {};
    }

    int rows = static_cast<int>(matrix.size());
    int cols = static_cast<int>(matrix[0].size());
    int rhsCols = static_cast<int>(right[0].size());

    Eigen::MatrixXd A(rows, cols);
    for (int i = 0; i < rows; ++i) {
        A.row(i) = Eigen::VectorXd::Map(&matrix[i][0], matrix[i].size());
    }

    Eigen::MatrixXd B(right.size(), rhsCols);
    for (int i = 0; i < right.size(); ++i) {
        B.row(i) = Eigen::VectorXd::Map(&right[i][0], right[i].size());
    }

    Eigen::MatrixXd X = A.colPivHouseholderQr().solve(B);

    std::vector<std::vector<double>> result(cols, std::vector<double>(rhsCols));
    for (int j = 0; j < rhsCols; ++j) {
        for (int i = 0; i < cols; ++i) {
            result[i][j] = X(i, j);
        }
    }
    return result;
}
bool LNLib::MathUtils::SolveLinearSystemBanded(int matrixDimension, const std::vector<std::vector<double>>& matrix, int bandwidth, const std::vector<std::vector<double>>& right, std::vector<std::vector<double>>& result)
{
    int sbw = bandwidth / 2;
    int n = matrixDimension;

    Eigen::SparseMatrix<double> A(n, n);
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(n * bandwidth);

    for (int i = 0; i < n; i++) 
    {
        for (int k = 0; k < bandwidth; k++) 
        {
            int j = i - sbw + k;
            if (j >= 0 && j < n) 
            {
                double value = matrix[i][k];
                if (value != 0.0)
                {
                    triplets.push_back(Eigen::Triplet<double>(i, j, value));
                }
            }
        }
    }
    A.setFromTriplets(triplets.begin(), triplets.end());

    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.analyzePattern(A);
    solver.factorize(A);

    if (solver.info() != Eigen::Success) {
        return false;
    }

    for (int k = 0; k < right[0].size(); k++)
    {
        Eigen::VectorXd rhs(n), solution(n);
        for (int i = 0; i < n; i++)
        {
            rhs(i) = right[i][k];
        }
        solution = solver.solve(rhs);
        if (solver.info() != Eigen::Success)
        {
            return false;
            break;
        }
        for (int i = 0; i < n; i++)
        {
            result[i][k] = solution(i);
        }
    }
    return true;
}


