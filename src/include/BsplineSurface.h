/**
 * Author:
 * 2023/06/08 - Yuqing Liang (BIMCoder Liang)
 * bim.frankliang@foxmail.com
 * 微信公众号：BIMCoder梁老师
 *
 * Use of this source code is governed by a GPL-3.0 license license that can be found in
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
		static void GetPointOnSurface(const std::vector<std::vector<T>>& controlPoints, const std::vector<double>& knotVectorU, const std::vector<double>& knotVectorV, unsigned int degreeU, unsigned int degreeV, UV uv, T& point)
		{
			int n = static_cast<int>(knotVectorU.size() - degreeU - 2);
			int uSpanIndex = Polynomials::GetKnotSpanIndex(n, degreeU, uv.GetU(), knotVectorU);

			std::vector<double> basisFunctionsU;
			Polynomials::BasisFunctions(uSpanIndex, degreeU, uv.GetU(), knotVectorU, basisFunctionsU);

			int m = static_cast<int>(knotVectorV.size() - degreeV - 2);
			int vSpanIndex = Polynomials::GetKnotSpanIndex(m, degreeV, uv.GetV(), knotVectorV);

			std::vector<double> basisFunctionsV;
			Polynomials::BasisFunctions(vSpanIndex, degreeV, uv.GetV(), knotVectorV, basisFunctionsV);

			int uind = uSpanIndex - degreeU;
			for (int l = 0; l <= static_cast<int>(degreeV); l++)
			{
				T temp = T();
				int vind = vSpanIndex - degreeV + 1;
				for (int k = 0; k <= static_cast<int>(degreeU); k++)
				{
					temp += basisFunctionsU[k] * controlPoints[uind + k][vind];
				}
				point += basisFunctionsV[l] * temp;
			}
		}

		/// <summary>
		/// The NURBS Book 2nd Edition Page111
		/// Algorithm A3.6
		/// Compute surface derivatives. (Usually Use)
		/// </summary>
		template <typename T>
		static void ComputeDerivatives(const std::vector<std::vector<T>>& controlPoints, const std::vector<double>& knotVectorU, const std::vector<double>& knotVectorV, unsigned int degreeU, unsigned int degreeV, UV uv, unsigned int derivative, std::vector<std::vector<T>>& derivatives)		
		{
			derivatives.resize(derivative + 1);
			for (int i = 0; i <= static_cast<int>(derivative); i++)
			{
				derivatives[i].resize(derivative + 1);
			}

			int du = std::min(derivative, degreeU);
			int dv = std::min(derivative, degreeV);

			int n = static_cast<int>(knotVectorU.size() - degreeU - 2);
			int uSpanIndex = Polynomials::GetKnotSpanIndex(n, degreeU, uv.GetU(), knotVectorU);
			std::vector<std::vector<double>> derivativeBasisFunctionsU;
			Polynomials::BasisFunctionsDerivatives(uSpanIndex, degreeU, uv.GetU(), du, knotVectorU, derivativeBasisFunctionsU);

			int m = static_cast<int>(knotVectorV.size() - degreeV - 2);
			int vSpanIndex = Polynomials::GetKnotSpanIndex(m, degreeV, uv.GetV(), knotVectorV);
			std::vector<std::vector<double>> derivativeBasisFunctionsV;
			Polynomials::BasisFunctionsDerivatives(vSpanIndex, degreeV, uv.GetV(), dv, knotVectorV, derivativeBasisFunctionsV);

			std::vector<T> temp;
			temp.resize(degreeV + 1);

			for (int k = 0; k <= du; k++)
			{
				for (int s = 0; s <= static_cast<int>(degreeV); s++)
				{
					for (int r = 0; r <= static_cast<int>(degreeU); r++)
					{
						temp[s] += derivativeBasisFunctionsU[k][r] * controlPoints[uSpanIndex - degreeU + r][vSpanIndex - degreeV + s];
					}
					int dd = std::min(static_cast<int>(derivative) - k, dv);
					for (int l = 0; l <= dd; l++)
					{
						for (int s = 0; s <= static_cast<int>(degreeV); s++)
						{
							derivatives[k][l] += derivativeBasisFunctionsV[l][s] * temp[s];
						}
					}
				}
			}
		}


		/// <summary>
		/// The NURBS Book 2nd Edition Page114.
		/// Algorithm A3.7
		/// Compute control points of derivative surfaces.
		/// </summary>
		static void ComputeControlPointsOfDerivatives(const std::vector<std::vector<XYZ>>& controlPoints, const std::vector<double>& knotVectorU, const std::vector<double>& knotVectorV, unsigned int degreeU, unsigned int degreeV, UV uv, unsigned int derMinU, unsigned int derMaxU, unsigned int derMinV, unsigned int derMaxV, unsigned int derivative, std::vector<std::vector<std::vector<std::vector<XYZ>>>>& controlPointsOfDerivative);

		/// <summary>
		/// The NURBS Book 2nd Edition Page115.
		/// Algorithm A3.8
		/// Compute surface derivatives.
		/// </summary>
		template <typename T>
		static void ComputeDerivativesByAllBasisFunctions(const std::vector<std::vector<T>>& controlPoints, const std::vector<double>& knotVectorU, const std::vector<double>& knotVectorV, unsigned int degreeU, unsigned int degreeV, UV uv, unsigned int derivative, std::vector<std::vector<T>>& derivatives)
		{
			derivatives.resize(derivative + 1);
			for (int i = 0; i <= derivative; i++)
			{
				derivatives[i].resize(derivative + 1);
			}

			int du = std::min(derivative, degreeU);
			int dv = std::min(derivative, degreeV);

			int n = static_cast<int>(knotVectorU.size() - degreeU) - 2;
			int uSpanIndex = Polynomials::GetKnotSpanIndex(n, degreeU, uv.GetU(), knotVectorU);
			std::vector<std::vector<double>> allBasisFunctionsU;
			allBasisFunctionsU.resize(degreeU + 1);
			for (int i = 0; i <= static_cast<int>(degreeU); i++)
			{
				std::vector<double> basisFunctions;
				Polynomials::BasisFunctions(uSpanIndex, i, uv.GetU(), knotVectorU, basisFunctions);
				allBasisFunctionsU[i] = basisFunctions;
			}

			int m = static_cast<int>(knotVectorV.size() - degreeV) - 2;
			int vSpanIndex = Polynomials::GetKnotSpanIndex(m, degreeV, uv.GetV(), knotVectorV);
			std::vector<std::vector<double>> allBasisFunctionsV;
			allBasisFunctionsV.resize(degreeV + 1);
			for (int i = 0; i <= static_cast<int>(degreeV); i++)
			{
				std::vector<double> basisFunctions;
				Polynomials::BasisFunctions(vSpanIndex, i, uv.GetV(), knotVectorV, basisFunctions);
				allBasisFunctionsV[i] = basisFunctions;
			}

			std::vector<std::vector<std::vector<std::vector<T>>>> controlPointsOfDerivative;
			ComputeControlPointsOfDerivatives(controlPoints, knotVectorU, knotVectorV, degreeU, degreeV, uv, uSpanIndex - degreeU, uSpanIndex, vSpanIndex - degreeV, vSpanIndex, derivative, controlPointsOfDerivative);

			for (int k = 0; k <= du; k++)
			{
				int dd = std::min(static_cast<int>(derivative - k), dv);
				for (int l = 0; l <= dd; l++)
				{
					for (int i = 0; i <= static_cast<int>(degreeV) - l; i++)
					{
						T temp = T();
						for (int j = 0; j <= static_cast<int>(degreeU) - k; j++)
						{
							temp += allBasisFunctionsU[j][degreeU - k] * controlPointsOfDerivative[k][l][j][i];
						}
						derivatives[k][l] += allBasisFunctionsV[i][degreeV - l] * temp;
					}
				}
			}
		}
	};
}


