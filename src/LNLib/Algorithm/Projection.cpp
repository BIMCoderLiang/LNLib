/*
 * Author:
 * 2023/06/29 - Yuqing Liang (BIMCoder Liang)
 * bim.frankliang@foxmail.com
 * 
 *
 * Use of this source code is governed by a GPL-3.0 license that can be found in
 * the LICENSE file.
 */

#include "Projection.h"
#include "XYZ.h"
#include "MathUtils.h"

LNLib::XYZ LNLib::Projection::PointToRay(const XYZ& origin, const XYZ& vector, const XYZ& Point)
{
    XYZ pt = Point - origin;
    double param = pt.DotProduct(vector);
    return origin + param * vector;
}

bool LNLib::Projection::PointToLine(const XYZ& start, const XYZ& end, const XYZ& point, XYZ& projectPoint)
{
    double length = start.Distance(end);
    if (MathUtils::IsAlmostEqualTo(length, 0.0))
    {
        return false;
    }
    double param = (point - start).DotProduct(end - start) / (length * length);
    if (MathUtils::IsGreaterThan(param, 1.0) || MathUtils::IsLessThan(param,0.0))
    {
        return false;
    }
    param = std::max(0.0, std::min(1.0, param));
    projectPoint = start + param * (end - start);
    return true;
}

LNLib::XYZ LNLib::Projection::Stereographic(const XYZ& pointOnSphere, double radius)
{
    double x = pointOnSphere.GetX();
    double y = pointOnSphere.GetY();
    double z = pointOnSphere.GetZ();

    double alpha = -2 * radius / (z - 2 * radius);
    double u = alpha * x;
    double v = alpha * y;

    return XYZ(u, v, 0);
}
