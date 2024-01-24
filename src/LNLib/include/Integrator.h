/*
 * Author:
 * 2024/01/08 - Yuqing Liang (BIMCoder Liang)
 * bim.frankliang@foxmail.com
 * 微信公众号：BIMCoder梁老师
 *
 * Use of this source code is governed by a GPL-3.0 license that can be found in
 * the LICENSE file.
 */

#pragma once
#include "LNLibDefinitions.h"
#include "LNObject.h"
#include <vector>

namespace LNLib
{
	class LNLIB_EXPORT IntegrationFunction
	{
	public:
		virtual double operator()(double parameter, void* customData) = 0;
	};

	class LNLIB_EXPORT Integrator
	{

	public:
		template<typename Function>
		static double Simpson(Function function, void* customData, double start, double end)
		{
			double st = (*function)(start, customData);
			double mt = (*function)((start + end) / 2.0, customData);
			double et = (*function)((end) / 2.0, customData);
			double result = ((end - start) / 6.0) * (st + 4 * mt + et);
			return result;
		}

		static double Simpson(double start, double end, std::vector<double> odds, std::vector<double> evens, double delta);

		/// <summary>
		/// According to https://github.com/Pomax/bezierjs
		/// Order is set 24.
		/// </summary>
		static const std::vector<double> GaussLegendreAbscissae;
		static const std::vector<double> GaussLegendreWeights;

		/// <summary>
		/// According to https://github.com/chrisidefix/nurbs
		/// </summary>
		static std::vector<double> ChebyshevSeries(int size = 100);
		template<typename Function>
		static double ClenshawCurtisQuadrature(Function function, void* customData, double start, double end, std::vector<double>& series, double epsilon = Constants::DistanceEpsilon)
		{
			double integration;
			int j, k, l;
			double err, esf, eref, erefh, hh, ir, iback, irback, ba, ss, x, y, fx, errir;
			int lenw = series.size() - 1;
			esf = 10;
			ba = 0.5 * (end - start);
			ss = 2 * series[lenw];
			x = ba * series[lenw];
			series[0] = 0.5 * (*function)(start, customData);
			series[3] = 0.5 * (*function)(end, customData);
			series[2] = (*function)(start + x, customData);
			series[4] = (*function)(end - x, customData);
			series[1] = (*function)(start + ba, customData);
			eref = 0.5 * (fabs(series[0]) + fabs(series[1]) + fabs(series[2]) + fabs(series[3]) + fabs(series[4]));
			series[0] += series[3];
			series[2] += series[4];
			ir = series[0] + series[1] + series[2];
			integration = series[0] * series[lenw - 1] + series[1] * series[lenw - 2] + series[2] * series[lenw - 3];
			erefh = eref * sqrt(epsilon);
			eref *= epsilon;
			hh = 0.25;
			l = 2;
			k = lenw - 5;
			do {
				iback = integration;
				irback = ir;
				x = ba * series[k + 1];
				y = 0;
				integration = series[0] * series[k];
				for (j = 1; j <= l; j++) {
					x += y;
					y += ss * (ba - x);
					fx = (*function)(start + x, customData) + (*function)(end - x, customData);
					ir += fx;
					integration += series[j] * series[k - j] + fx * series[k - j - l];
					series[j + l] = fx;
				}
				ss = 2 * series[k + 1];
				err = esf * l * fabs(integration - iback);
				hh *= 0.25;
				errir = hh * fabs(ir - 2 * irback);
				l *= 2;
				k -= l + 2;
			} while ((err > erefh || errir > eref) && k > 4 * l);
			integration *= end - start;
			if (err > erefh || errir > eref)
			{
				err *= -fabs(end - start);
			}
			else
			{
				err = eref * fabs(end - start);
			}
			return integration;
		}
	};
}


