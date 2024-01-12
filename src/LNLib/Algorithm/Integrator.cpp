/*
 * Author:
 * 2024/01/08 - Yuqing Liang (BIMCoder Liang)
 * bim.frankliang@foxmail.com
 * 微信公众号：BIMCoder梁老师
 *
 * Use of this source code is governed by a GPL-3.0 license that can be found in
 * the LICENSE file.
 */

#include "Integrator.h"

namespace LNLib
{
	double Integrator::Simpson(double start, double end, double startTangentLength, double middleTangentLength, double endTangentLength)
	{
		double result = ((end - start) / 6.0) * (startTangentLength + 4 * middleTangentLength + endTangentLength);
		return result;
	}

	double Integrator::Simpson(double start, double end, std::vector<double> odds, std::vector<double> evens, double delta)
	{
		double oddsSum = 0.0;
		double evensSum = 0.0;
		for (int i = 0; i < odds.size(); i++)
		{
			oddsSum += 4 * odds[i];
		}
		for (int i = 0; i < evens.size(); i++)
		{
			evensSum += 2 * evens[i];
		}
		double result = (delta / 3.0) * (start + oddsSum + evensSum + end);
		return result;
	}
}


