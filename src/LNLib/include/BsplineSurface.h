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
		static T GetPointOnSurface(const std::vector<std::vector<T>>& controlPoints, const std::vector<double>& knotVectorU, const std::vector<double>& knotVectorV, int degreeU, int degreeV, UV uv)
		{
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
		static std::vector<std::vector<T>> ComputeDerivatives(const std::vector<std::vector<T>>& controlPoints, const std::vector<double>& knotVectorU, const std::vector<double>& knotVectorV, int degreeU, int degreeV, UV uv, int derivative)
		{
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
		static std::vector<std::vector<std::vector<std::vector<XYZ>>>> ComputeControlPointsOfDerivatives(const std::vector<std::vector<XYZ>>& controlPoints, const std::vector<double>& knotVectorU, const std::vector<double>& knotVectorV, int degreeU, int degreeV, UV uv, int minSpanIndexU, int maxSpanIndexU, int minSpanIndexV, int maxSpanIndexV, int derivative);

		/// <summary>
		/// The NURBS Book 2nd Edition Page115.
		/// Algorithm A3.8
		/// Compute surface derivatives.
		/// </summary>
		template <typename T>
		static std::vector<std::vector<T>> ComputeDerivativesByAllBasisFunctions(const std::vector<std::vector<T>>& controlPoints, const std::vector<double>& knotVectorU, const std::vector<double>& knotVectorV, unsigned int degreeU, unsigned int degreeV, UV uv, unsigned int derivative)
		{
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
				std::vector<double> basisFunctions = Polynomials::BasisFunctions(vSpanIndex, i, knotVectorV, uv.GetV());
				allBasisFunctionsV[i] = basisFunctions;
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


