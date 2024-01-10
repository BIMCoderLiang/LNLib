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
}


