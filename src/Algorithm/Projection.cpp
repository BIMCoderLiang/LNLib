/*
 * Author:
 * 2023/06/29 - Yuqing Liang (BIMCoder Liang)
 * bim.frankliang@foxmail.com
 * 微信公众号：BIMCoder梁老师
 *
 * Use of this source code is governed by a GPL-3.0 license license that can be found in
 * the LICENSE file.
 */

#include "Projection.h"
#include "XYZ.h"

void LNLib::Projection::PointToLine(const XYZ& origin, const XYZ& vector, const XYZ& Point, XYZ& ProjectPoint)
{
    XYZ temp = vector;
    XYZ nTemp = temp.Normalize();

    ProjectPoint = origin + ((Point - origin).DotProduct(nTemp)) * nTemp;
}
