/*
 * Author:
 * 2025/11/23 - Yuqing Liang (BIMCoder Liang)
 * bim.frankliang@foxmail.com
 * 
 *
 * Use of this source code is governed by a LGPL-2.1 license that can be found in
 * the LICENSE file.
 */

#include <pybind11/pybind11.h>
#include "XYZ.h"
#include "Intersection.h"

namespace py = pybind11;
void cstrIntersection(py::module_&m)
{
    py::class_<LNLib::Intersection>(m, "Intersection")
        .def_static("ComputeRays", [](const LNLib::XYZ& point0, const LNLib::XYZ& vector0, const LNLib::XYZ& point1, const LNLib::XYZ& vector1) {
        double param0, param1;
        LNLib::XYZ intersectPoint;
        LNLib::CurveCurveIntersectionType type = LNLib::Intersection::ComputeRays(point0, vector0, point1, vector1, param0, param1, intersectPoint);
        return py::make_tuple(type, param0, param1, intersectPoint);
            })
        .def_static("ComputeLineAndPlane", [](const LNLib::XYZ& normal, const LNLib::XYZ& pointOnPlane, const LNLib::XYZ& pointOnLine, const LNLib::XYZ& lineDirection) {
        LNLib::XYZ intersectPoint;
        LNLib::LinePlaneIntersectionType type = LNLib::Intersection::ComputeLineAndPlane(normal, pointOnPlane, pointOnLine, lineDirection, intersectPoint);
        return py::make_tuple(type, intersectPoint);
            });
}
