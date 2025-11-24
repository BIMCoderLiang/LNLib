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
        .def_static("ComputeRays", &LNLib::Intersection::ComputeRays)
        .def_static("ComputeLineAndPlane", &LNLib::Intersection::ComputeLineAndPlane);
}
