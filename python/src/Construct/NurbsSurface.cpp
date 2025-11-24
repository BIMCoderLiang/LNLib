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
#include "NurbsSurface.h"

namespace py = pybind11;
void cstrNurbsSurface(py::module_&m)
{
    py::class_<LNLib::NurbsSurface>(m, "NurbsSurface")
        .def_static("GetPointOnSurface", &LNLib::NurbsSurface::GetPointOnSurface)
        .def_static("ComputeRationalSurfaceDerivatives", &LNLib::NurbsSurface::ComputeRationalSurfaceDerivatives)
        .def_static("ComputeRationalSurfaceFirstOrderDerivative", &LNLib::NurbsSurface::ComputeRationalSurfaceFirstOrderDerivative)
        .def_static("Curvature", &LNLib::NurbsSurface::Curvature)
        .def_static("Normal", &LNLib::NurbsSurface::Normal)
        .def_static("Swap", &LNLib::NurbsSurface::Swap)
        .def_static("Reverse", &LNLib::NurbsSurface::Reverse)
        .def_static("InsertKnot", &LNLib::NurbsSurface::InsertKnot)
        .def_static("RefineKnotVector", &LNLib::NurbsSurface::RefineKnotVector)
        .def_static("DecomposeToBeziers", &LNLib::NurbsSurface::DecomposeToBeziers)
        .def_static("RemoveKnot", &LNLib::NurbsSurface::RemoveKnot)
        .def_static("ElevateDegree", &LNLib::NurbsSurface::ElevateDegree)
        .def_static("ReduceDegree", &LNLib::NurbsSurface::ReduceDegree)
        .def_static("EquallyTessellate", &LNLib::NurbsSurface::EquallyTessellate)
        .def_static("IsClosed", &LNLib::NurbsSurface::IsClosed)
        .def_static("GetParamOnSurface", &LNLib::NurbsSurface::GetParamOnSurface)
        .def_static("GetParamOnSurfaceByGSA", &LNLib::NurbsSurface::GetParamOnSurfaceByGSA)
        .def_static("Reparametrize", &LNLib::NurbsSurface::Reparametrize)
        .def_static("GetUVTangent", &LNLib::NurbsSurface::GetUVTangent)
        .def_static("CreateBilinearSurface", &LNLib::NurbsSurface::CreateBilinearSurface)
        .def_static("CreateCylindricalSurface", &LNLib::NurbsSurface::CreateCylindricalSurface)
        .def_static("CreateRuledSurface", &LNLib::NurbsSurface::CreateRuledSurface)
        .def_static("CreateRevolvedSurface", &LNLib::NurbsSurface::CreateRevolvedSurface)
        .def_static("NonUniformScaling", &LNLib::NurbsSurface::NonUniformScaling)
        .def_static("MakeCornerFilletSurface", &LNLib::NurbsSurface::MakeCornerFilletSurface)
        .def_static("GlobalInterpolation", &LNLib::NurbsSurface::GlobalInterpolation)
        .def_static("BicubicLocalInterpolation", &LNLib::NurbsSurface::BicubicLocalInterpolation)
        .def_static("GlobalApproximation", &LNLib::NurbsSurface::GlobalApproximation)
        .def_static("CreateSwungSurface", &LNLib::NurbsSurface::CreateSwungSurface)
        .def_static("CreateLoftSurface", &LNLib::NurbsSurface::CreateLoftSurface)
        .def_static("CreateGeneralizedTranslationalSweepSurface", &LNLib::NurbsSurface::CreateGeneralizedTranslationalSweepSurface)
        .def_static("CreateSweepSurface", py::overload_cast<const LNLib::LN_NurbsCurve&, const LNLib::LN_NurbsCurve&, int, LNLib::LN_NurbsSurface&>(&LNLib::NurbsSurface::CreateSweepSurface))
        .def_static("CreateSweepSurface", py::overload_cast<const LNLib::LN_NurbsCurve&, const LNLib::LN_NurbsCurve&, int, int, LNLib::LN_NurbsSurface&>(&LNLib::NurbsSurface::CreateSweepSurface))
        .def_static("CreateGordonSurface", &LNLib::NurbsSurface::CreateGordonSurface)
        .def_static("CreateCoonsSurface", &LNLib::NurbsSurface::CreateCoonsSurface)
        .def_static("ApproximateArea", &LNLib::NurbsSurface::ApproximateArea)
        .def_static("Triangulate", &LNLib::NurbsSurface::Triangulate);
}
