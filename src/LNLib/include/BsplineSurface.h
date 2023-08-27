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
#include "Polynomials.h"
#include "UV.h"
#include "ValidationUtils.h"
#include "LNLibExceptions.h"
#include <vector>

namespace LNLib
{
	class XYZ;
	class XYZW;
	class LNLIB_EXPORT BsplineSurface
	{
	public:

		/// <summary>
		/// The NURBS Book 2nd Edition Page103
		/// Algorithm A3.5
		/// Compute surface point.
		/// </summary>
		template <typename T>
		static T GetPointOnSurface(int degreeU, int degreeV, const std::vector<double>& knotVectorU, const std::vector<double>& knotVectorV, UV uv, const std::vector<std::vector<T>>& controlPoints)
		{
			VALIDATE_ARGUMENT(degreeU > 0, "degreeU", "Degree must greater than zero.");
			VALIDATE_ARGUMENT(degreeV > 0, "degreeV", "Degree must greater than zero.");
			VALIDATE_ARGUMENT(knotVectorU.size() > 0, "knotVectorU", "KnotVector size must greater than zero.");
			VALIDATE_ARGUMENT(ValidationUtils::IsValidKnotVector(knotVectorU), "knotVectorU", "KnotVector must be a nondecreasing sequence of real numbers.");
			VALIDATE_ARGUMENT_RANGE(uv.GetU(), knotVectorU[0], knotVectorU[knotVectorU.size() - 1]);
			VALIDATE_ARGUMENT(knotVectorV.size() > 0, "knotVectorV", "KnotVector size must greater than zero.");
			VALIDATE_ARGUMENT(ValidationUtils::IsValidKnotVector(knotVectorV), "knotVectorV", "KnotVector must be a nondecreasing sequence of real numbers.");
			VALIDATE_ARGUMENT_RANGE(uv.GetV(), knotVectorV[0], knotVectorV[knotVectorV.size() - 1]);
			VALIDATE_ARGUMENT(controlPoints.size() > 0, "controlPoints", "ControlPoints must contains one point at least.");
			VALIDATE_ARGUMENT(ValidationUtils::IsValidBspline(degreeU, knotVectorU.size(), controlPoints.size()), "controlPoints", "Arguments must fit: m = n + p + 1");
			VALIDATE_ARGUMENT(ValidationUtils::IsValidBspline(degreeV, knotVectorV.size(), controlPoints[0].size()), "controlPoints", "Arguments must fit: m = n + p + 1");

			int uSpanIndex = Polynomials::GetKnotSpanIndex(degreeU, knotVectorU, uv.GetU());
			std::vector<double> basisFunctionsU = Polynomials::BasisFunctions(uSpanIndex, degreeU, knotVectorU, uv.GetU());

			int vSpanIndex = Polynomials::GetKnotSpanIndex(degreeV, knotVectorV, uv.GetV());
			std::vector<double> basisFunctionsV = Polynomials::BasisFunctions(vSpanIndex, degreeV, knotVectorV, uv.GetV());

			int uind = uSpanIndex - degreeU;
			T point;
			for (int l = 0; l <= degreeV; l++)
			{
				T temp = T();
				int vind = vSpanIndex - degreeV + 1;
				for (int k = 0; k <= degreeU; k++)
				{
					temp += basisFunctionsU[k] * controlPoints[uind + k][vind];
				}
				point += basisFunctionsV[l] * temp;
			}
			return point;
		}

		/// <summary>
		/// The NURBS Book 2nd Edition Page111
		/// Algorithm A3.6
		/// Compute surface derivatives. (Usually Use)
		/// </summary>
		template <typename T>
		static std::vector<std::vector<T>> ComputeDerivatives(int degreeU, int degreeV, int derivative, const std::vector<double>& knotVectorU, const std::vector<double>& knotVectorV, UV uv, const std::vector<std::vector<T>>& controlPoints)
		{
			VALIDATE_ARGUMENT(degreeU > 0, "degreeU", "Degree must greater than zero.");
			VALIDATE_ARGUMENT(degreeV > 0, "degreeU", "Degree must greater than zero.");
			VALIDATE_ARGUMENT(derivative > 0, "derivative", "derivative must greater than zero.");
			VALIDATE_ARGUMENT(knotVectorU.size() > 0, "knotVectorU", "KnotVector size must greater than zero.");
			VALIDATE_ARGUMENT(ValidationUtils::IsValidKnotVector(knotVectorU), "knotVectorU", "KnotVector must be a nondecreasing sequence of real numbers.");
			VALIDATE_ARGUMENT_RANGE(uv.GetU(), knotVectorU[0], knotVectorU[knotVectorU.size() - 1]);
			VALIDATE_ARGUMENT(knotVectorV.size() > 0, "knotVectorV", "KnotVector size must greater than zero.");
			VALIDATE_ARGUMENT(ValidationUtils::IsValidKnotVector(knotVectorV), "knotVectorV", "KnotVector must be a nondecreasing sequence of real numbers.");
			VALIDATE_ARGUMENT_RANGE(uv.GetV(), knotVectorV[0], knotVectorV[knotVectorV.size() - 1]);
			VALIDATE_ARGUMENT(controlPoints.size() > 0, "controlPoints", "ControlPoints must contains one point at least.");
			VALIDATE_ARGUMENT(ValidationUtils::IsValidBspline(degreeU, knotVectorU.size(), controlPoints.size()), "controlPoints", "Arguments must fit: m = n + p + 1");
			VALIDATE_ARGUMENT(ValidationUtils::IsValidBspline(degreeV, knotVectorV.size(), controlPoints[0].size()), "controlPoints", "Arguments must fit: m = n + p + 1");

			std::vector<std::vector<T>> derivatives(derivative + 1, std::vector<T>(derivative + 1));

			int du = std::min(derivative, degreeU);
			int dv = std::min(derivative, degreeV);

			int uSpanIndex = Polynomials::GetKnotSpanIndex(degreeU, knotVectorU, uv.GetU());
			std::vector<std::vector<double>> derivativeBasisFunctionsU = Polynomials::BasisFunctionsDerivatives(uSpanIndex, degreeU, du, knotVectorU, uv.GetU());

			int vSpanIndex = Polynomials::GetKnotSpanIndex(degreeV, knotVectorV, uv.GetV());
			std::vector<std::vector<double>> derivativeBasisFunctionsV = Polynomials::BasisFunctionsDerivatives(vSpanIndex, degreeV, dv, knotVectorV, uv.GetV());

			std::vector<T> temp(degreeV + 1);

			for (int k = 0; k <= du; k++)
			{
				for (int s = 0; s <= degreeV; s++)
				{
					for (int r = 0; r <= degreeU; r++)
					{
						temp[s] += derivativeBasisFunctionsU[k][r] * controlPoints[uSpanIndex - degreeU + r][vSpanIndex - degreeV + s];
					}
					int dd = std::min(derivative - k, dv);
					for (int l = 0; l <= dd; l++)
					{
						for (int s = 0; s <= degreeV; s++)
						{
							derivatives[k][l] += derivativeBasisFunctionsV[l][s] * temp[s];
						}
					}
				}
			}
			return derivatives;
		}

		/// <summary>
		/// The NURBS Book 2nd Edition Page114.
		/// Algorithm A3.7
		/// Compute control points of derivative surfaces.
		/// </summary>
		template <typename T>
		static std::vector<std::vector<std::vector<std::vector<T>>>> ComputeControlPointsOfDerivatives(int degreeU, int degreeV, int derivative, int minSpanIndexU, int maxSpanIndexU, int minSpanIndexV, int maxSpanIndexV, const std::vector<double>& knotVectorU, const std::vector<double>& knotVectorV, UV uv, const std::vector<std::vector<T>>& controlPoints)
		{
			VALIDATE_ARGUMENT(degreeU > 0, "degreeU", "Degree must greater than zero.");
			VALIDATE_ARGUMENT(degreeV > 0, "degreeU", "Degree must greater than zero.");
			VALIDATE_ARGUMENT(derivative > 0, "derivative", "derivative must greater than zero.");
			VALIDATE_ARGUMENT_RANGE(minSpanIndexU, 0, maxSpanIndexU);
			VALIDATE_ARGUMENT_RANGE(minSpanIndexV, 0, maxSpanIndexV);
			VALIDATE_ARGUMENT(knotVectorU.size() > 0, "knotVectorU", "KnotVector size must greater than zero.");
			VALIDATE_ARGUMENT(ValidationUtils::IsValidKnotVector(knotVectorU), "knotVectorU", "KnotVector must be a nondecreasing sequence of real numbers.");
			VALIDATE_ARGUMENT_RANGE(uv.GetU(), knotVectorU[0], knotVectorU[knotVectorU.size() - 1]);
			VALIDATE_ARGUMENT(knotVectorV.size() > 0, "knotVectorV", "KnotVector size must greater than zero.");
			VALIDATE_ARGUMENT(ValidationUtils::IsValidKnotVector(knotVectorV), "knotVectorV", "KnotVector must be a nondecreasing sequence of real numbers.");
			VALIDATE_ARGUMENT_RANGE(uv.GetV(), knotVectorV[0], knotVectorV[knotVectorV.size() - 1]);
			VALIDATE_ARGUMENT(controlPoints.size() > 0, "controlPoints", "ControlPoints must contains one point at least.");
			VALIDATE_ARGUMENT(ValidationUtils::IsValidBspline(degreeU, knotVectorU.size(), controlPoints.size()), "controlPoints", "Arguments must fit: m = n + p + 1");
			VALIDATE_ARGUMENT(ValidationUtils::IsValidBspline(degreeV, knotVectorV.size(), controlPoints[0].size()), "controlPoints", "Arguments must fit: m = n + p + 1");

			std::vector<std::vector<std::vector<std::vector<T>>>> controlPointsOfDerivative(derivative + 1,
				std::vector<std::vector<std::vector<T>>>(maxSpanIndexU - minSpanIndexU + 1,
					std::vector<std::vector<T>>(maxSpanIndexV - minSpanIndexV + 1)));

			int du = std::min(derivative, degreeU);
			int dv = std::min(derivative, degreeV);
			int rangeU = maxSpanIndexU - minSpanIndexU;
			int rangeV = maxSpanIndexV - minSpanIndexV;

			for (int j = minSpanIndexV; j <= maxSpanIndexV; j++)
			{
				std::vector<T> points;
				for (int i = 0; i < controlPoints.size(); i++)
				{
					points.emplace_back(controlPoints[i][j]);
				}

				std::vector<std::vector<T>> temp = BsplineCurve::ComputeControlPointsOfDerivatives(degreeU, du, minSpanIndexU, maxSpanIndexU, knotVectorU, points);
				for (int k = 0; k <= du; k++)
				{
					for (int i = 0; i <= rangeU - k; i++)
					{
						controlPointsOfDerivative[k][0][i][j - minSpanIndexV] = temp[k][i];
					}
				}
			}
			for (int k = 0; k < du; k++)
			{
				for (int i = 0; i <= rangeU - k; i++)
				{
					std::vector<T> points = controlPointsOfDerivative[k][0][i];
					int dd = std::min(derivative - k, dv);
					std::vector<std::vector<T>> temp = BsplineCurve::ComputeControlPointsOfDerivatives(degreeV, dd, 0, rangeV, knotVectorV, points);
					for (int l = 1; l <= dd; l++)
					{
						for (int j = 0; j < rangeV - l; j++)
						{
							controlPointsOfDerivative[k][l][i][j] = temp[l][j];
						}
					}
				}
			}
			return controlPointsOfDerivative;
		}

		/// <summary>
		/// The NURBS Book 2nd Edition Page115.
		/// Algorithm A3.8
		/// Compute surface derivatives.
		/// </summary>
		template <typename T>
		static std::vector<std::vector<T>> ComputeDerivativesByAllBasisFunctions(int degreeU, int degreeV, int derivative, const std::vector<double>& knotVectorU, const std::vector<double>& knotVectorV, UV uv, const std::vector<std::vector<T>>& controlPoints)
		{
			VALIDATE_ARGUMENT(degreeU > 0, "degreeU", "Degree must greater than zero.");
			VALIDATE_ARGUMENT(degreeV > 0, "degreeU", "Degree must greater than zero.");
			VALIDATE_ARGUMENT(derivative > 0, "derivative", "derivative must greater than zero.");
			VALIDATE_ARGUMENT(knotVectorU.size() > 0, "knotVectorU", "KnotVector size must greater than zero.");
			VALIDATE_ARGUMENT(ValidationUtils::IsValidKnotVector(knotVectorU), "knotVectorU", "KnotVector must be a nondecreasing sequence of real numbers.");
			VALIDATE_ARGUMENT_RANGE(uv.GetU(), knotVectorU[0], knotVectorU[knotVectorU.size() - 1]);
			VALIDATE_ARGUMENT(knotVectorV.size() > 0, "knotVectorV", "KnotVector size must greater than zero.");
			VALIDATE_ARGUMENT(ValidationUtils::IsValidKnotVector(knotVectorV), "knotVectorV", "KnotVector must be a nondecreasing sequence of real numbers.");
			VALIDATE_ARGUMENT_RANGE(uv.GetV(), knotVectorV[0], knotVectorV[knotVectorV.size() - 1]);
			VALIDATE_ARGUMENT(controlPoints.size() > 0, "controlPoints", "ControlPoints must contains one point at least.");
			VALIDATE_ARGUMENT(ValidationUtils::IsValidBspline(degreeU, knotVectorU.size(), controlPoints.size()), "controlPoints", "Arguments must fit: m = n + p + 1");
			VALIDATE_ARGUMENT(ValidationUtils::IsValidBspline(degreeV, knotVectorV.size(), controlPoints[0].size()), "controlPoints", "Arguments must fit: m = n + p + 1");

			std::vector<std::vector<T>> derivatives(derivative + 1, std::vector<T>(derivative + 1));

			int du = std::min(derivative, degreeU);
			int dv = std::min(derivative, degreeV);

			int uSpanIndex = Polynomials::GetKnotSpanIndex(degreeU, uv.GetU(), knotVectorU);
			std::vector<std::vector<double>> allBasisFunctionsU(degreeU + 1);

			for (int i = 0; i <= degreeU; i++)
			{
				std::vector<double> basisFunctions = Polynomials::BasisFunctions(uSpanIndex, i, knotVectorU, uv.GetU());
				allBasisFunctionsU[i] = basisFunctions;
			}

			int vSpanIndex = Polynomials::GetKnotSpanIndex(degreeV, uv.GetV(), knotVectorV);
			std::vector<std::vector<double>> allBasisFunctionsV(degreeV + 1);
			for (int i = 0; i <= degreeV; i++)
			{
				allBasisFunctionsV[i] = Polynomials::BasisFunctions(vSpanIndex, i, knotVectorV, uv.GetV());
			}

			std::vector<std::vector<std::vector<std::vector<T>>>> controlPointsOfDerivative = ComputeControlPointsOfDerivatives(controlPoints, knotVectorU, knotVectorV, degreeU, degreeV, uv, uSpanIndex - degreeU, uSpanIndex, vSpanIndex - degreeV, vSpanIndex, derivative);

			for (int k = 0; k <= du; k++)
			{
				int dd = std::min(derivative - k, dv);
				for (int l = 0; l <= dd; l++)
				{
					for (int i = 0; i <= degreeV - l; i++)
					{
						T temp = T();
						for (int j = 0; j <= degreeU - k; j++)
						{
							temp += allBasisFunctionsU[j][degreeU - k] * controlPointsOfDerivative[k][l][j][i];
						}
						derivatives[k][l] += allBasisFunctionsV[i][degreeV - l] * temp;
					}
				}
			}
			return derivatives;
		}
	};
}


