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
#include "KnotVectorUtils.h"

namespace py = pybind11;
void cstrKnotVectorUtils(py::module_&m)
{
    py::class_<LNLib::KnotVectorUtils>(m, "KnotVectorUtils")
        .def_static("GetContinuity", &LNLib::KnotVectorUtils::GetContinuity)
        .def_static("Rescale", &LNLib::KnotVectorUtils::Rescale)
        .def_static("GetInsertedKnotElement", py::overload_cast<int, const std::vector<double>&, double, double>(&LNLib::KnotVectorUtils::GetInsertedKnotElement))
        .def_static("GetKnotMultiplicityMap", &LNLib::KnotVectorUtils::GetKnotMultiplicityMap)
        .def_static("GetInternalKnotMultiplicityMap", &LNLib::KnotVectorUtils::GetInternalKnotMultiplicityMap)
        .def_static("GetInsertedKnotElement", [](const std::vector<double>& knotVector0, const std::vector<double>& knotVector1) {
        std::vector<double> insertElements0, insertElements1;
        LNLib::KnotVectorUtils::GetInsertedKnotElement(knotVector0, knotVector1, insertElements0, insertElements1);
        return py::make_tuple(insertElements0, insertElements1);
            })
        .def_static("GetInsertedKnotElements", py::overload_cast<const std::vector<std::vector<double>>&>(&LNLib::KnotVectorUtils::GetInsertedKnotElements))
        .def_static("GetInsertedKnotElements", py::overload_cast<int, const std::vector<double>&>(&LNLib::KnotVectorUtils::GetInsertedKnotElements))
        .def_static("IsUniform", &LNLib::KnotVectorUtils::IsUniform);
}
