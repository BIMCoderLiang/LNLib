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
#include "ValidationUtils.h"

namespace py = pybind11;
void cstrValidationUtils(py::module_&m)
{
    py::class_<LNLib::ValidationUtils>(m, "ValidationUtils")
        .def_static("IsValidBezier", &LNLib::ValidationUtils::IsValidBezier)
        .def_static("IsValidKnotVector", &LNLib::ValidationUtils::IsValidKnotVector)
        .def_static("IsValidBspline", &LNLib::ValidationUtils::IsValidBspline)
        .def_static("IsValidNurbs", &LNLib::ValidationUtils::IsValidNurbs)
        .def_static("IsValidDegreeReduction", &LNLib::ValidationUtils::IsValidDegreeReduction)
        .def_static("ComputeCurveModifyTolerance", &LNLib::ValidationUtils::ComputeCurveModifyTolerance);
}
