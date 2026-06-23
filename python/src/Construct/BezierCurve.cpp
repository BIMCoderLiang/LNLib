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
#include "LNObject.h"
#include "BezierCurve.h"

namespace py = pybind11;
void cstrBezierCurve(py::module_&m)
{
    py::class_<LNLib::BezierCurve>(m, "BezierCurve")
        .def_static("GetPointOnCurveByBernstein_XYZ", &LNLib::BezierCurve::GetPointOnCurveByBernstein<LNLib::XYZ>)
        .def_static("GetPointOnCurveByBernstein_XYZW", &LNLib::BezierCurve::GetPointOnCurveByBernstein<LNLib::XYZW>)
        .def_static("GetPointOnCurveByDeCasteljau_XYZ", &LNLib::BezierCurve::GetPointOnCurveByDeCasteljau<LNLib::XYZ>)
        .def_static("GetPointOnCurveByDeCasteljau_XYZW", &LNLib::BezierCurve::GetPointOnCurveByDeCasteljau<LNLib::XYZW>);
}
