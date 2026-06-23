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
#include "BezierSurface.h"

namespace py = pybind11;
void cstrBezierSurface(py::module_&m)
{
    py::class_<LNLib::BezierSurface>(m, "BezierSurface")
        .def_static("GetPointOnSurfaceByDeCasteljau_XYZ", &LNLib::BezierSurface::GetPointOnSurfaceByDeCasteljau<LNLib::XYZ>)
        .def_static("GetPointOnSurfaceByDeCasteljau_XYZW", &LNLib::BezierSurface::GetPointOnSurfaceByDeCasteljau<LNLib::XYZW>);

}
