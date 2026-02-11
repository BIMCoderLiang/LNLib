/*
 * Author:
 * 2023/06/29 - Yuqing Liang (BIMCoder Liang)
 * bim.frankliang@foxmail.com
 * 
 *
 * Use of this source code is governed by a LGPL-2.1 license that can be found in
 * the LICENSE file.
 */

#include "Projection.h"
#include "XYZ.h"
#include "MathUtils.h"

LNLib::XYZ LNLib::Projection::PointToRay(const XYZ& origin, const XYZ& vector, const XYZ& point)
{
    XYZ diff = point - origin;
    double dot = diff.DotProduct(vector);
    double lenSqr = vector.DotProduct(vector);

    double t = dot / lenSqr;
    return origin + t * vector;
}

bool LNLib::Projection::PointToLine(const XYZ& start, const XYZ& end, const XYZ& point, XYZ& projectPoint)
{
    XYZ direction = end - start;
    double lenSqr = direction.DotProduct(direction);

    if (MathUtils::IsAlmostEqualTo(lenSqr, 0.0)) {
        return false;
    }

    double t = (point - start).DotProduct(direction) / lenSqr;

    if (MathUtils::IsLessThan(t, 0.0) || 
        MathUtils::IsGreaterThan(t, 1.0)) {
        return false;
    }

    projectPoint = start + t * direction;
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
