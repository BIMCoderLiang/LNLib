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
#include "LNEnums.h"

namespace py = pybind11;
void cstrEnums(py::module_&m)
{
    py::enum_<LNLib::CurveCurveIntersectionType>(m, "CurveCurveIntersectionType")
        .value("Intersecting", LNLib::CurveCurveIntersectionType::Intersecting)
        .value("Parallel", LNLib::CurveCurveIntersectionType::Parallel)
        .value("Coincident", LNLib::CurveCurveIntersectionType::Coincident)
        .value("Skew", LNLib::CurveCurveIntersectionType::Skew);

    py::enum_<LNLib::LinePlaneIntersectionType>(m, "LinePlaneIntersectionType")
        .value("Intersecting", LNLib::LinePlaneIntersectionType::Intersecting)
        .value("Parallel", LNLib::LinePlaneIntersectionType::Parallel)
        .value("On", LNLib::LinePlaneIntersectionType::On);

    py::enum_<LNLib::CurveNormal>(m, "CurveNormal")
        .value("Normal", LNLib::CurveNormal::Normal)
        .value("Binormal", LNLib::CurveNormal::Binormal);

    py::enum_<LNLib::SurfaceDirection>(m, "SurfaceDirection")
        .value("All", LNLib::SurfaceDirection::All)
        .value("UDirection", LNLib::SurfaceDirection::UDirection)
        .value("VDirection", LNLib::SurfaceDirection::VDirection);

    py::enum_<LNLib::SurfaceCurvature>(m, "SurfaceCurvature")
        .value("Maximum", LNLib::SurfaceCurvature::Maximum)
        .value("Minimum", LNLib::SurfaceCurvature::Minimum)
        .value("Gauss", LNLib::SurfaceCurvature::Gauss)
        .value("Mean", LNLib::SurfaceCurvature::Mean)
        .value("Abs", LNLib::SurfaceCurvature::Abs)
        .value("Rms", LNLib::SurfaceCurvature::Rms);

    py::enum_<LNLib::IntegratorType>(m, "IntegratorType")
        .value("Simpson", LNLib::IntegratorType::Simpson)
        .value("GaussLegendre", LNLib::IntegratorType::GaussLegendre)
        .value("Chebyshev", LNLib::IntegratorType::Chebyshev);

    py::enum_<LNLib::OffsetType>(m, "OffsetType")
        .value("TillerAndHanson", LNLib::OffsetType::TillerAndHanson)
        .value("PieglAndTiller", LNLib::OffsetType::PieglAndTiller);
}
