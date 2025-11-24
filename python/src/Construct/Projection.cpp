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
#include "Projection.h"

namespace py = pybind11;
void cstrProjection(py::module_&m)
{
    py::class_<LNLib::Projection>(m, "Projection")
        .def_static("PointToRay", &LNLib::Projection::PointToRay)
        .def_static("PointToLine", &LNLib::Projection::PointToLine)
        .def_static("Stereographic", &LNLib::Projection::Stereographic);
}
