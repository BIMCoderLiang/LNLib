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
#include "XYZ.h"
#include "Matrix4d.h"
#include "NurbsCurve.h"

namespace py = pybind11;
void cstrNurbsCurve(py::module_&m)
{
    py::class_<LNLib::NurbsCurve>(m, "NurbsCurve")
        .def_static("GetPointOnCurve", &LNLib::NurbsCurve::GetPointOnCurve)
        .def_static("ComputeRationalCurveDerivatives", &LNLib::NurbsCurve::ComputeRationalCurveDerivatives)
        .def_static("CanComputerDerivative", &LNLib::NurbsCurve::CanComputerDerivative)
        .def_static("Curvature", &LNLib::NurbsCurve::Curvature)
        .def_static("Torsion", &LNLib::NurbsCurve::Torsion)
        .def_static("InsertKnot", &LNLib::NurbsCurve::InsertKnot)
        .def_static("GetPointOnCurveByCornerCut", &LNLib::NurbsCurve::GetPointOnCurveByCornerCut)
        .def_static("RefineKnotVector", &LNLib::NurbsCurve::RefineKnotVector)
        .def_static("DecomposeToBeziers", &LNLib::NurbsCurve::DecomposeToBeziers)
        .def_static("RemoveKnot", &LNLib::NurbsCurve::RemoveKnot)
        .def_static("RemoveExcessiveKnots", &LNLib::NurbsCurve::RemoveExcessiveKnots)
        .def_static("ElevateDegree", &LNLib::NurbsCurve::ElevateDegree)
        .def_static("ReduceDegree", &LNLib::NurbsCurve::ReduceDegree)
        .def_static("EquallyTessellate", &LNLib::NurbsCurve::EquallyTessellate)
        .def_static("IsClosed", &LNLib::NurbsCurve::IsClosed)
        .def_static("GetParamOnCurve", py::overload_cast<const LNLib::LN_NurbsCurve&, const LNLib::XYZ&>(&LNLib::NurbsCurve::GetParamOnCurve))
        .def_static("CreateTransformed", &LNLib::NurbsCurve::CreateTransformed)
        .def_static("Reparametrize", py::overload_cast<const LNLib::LN_NurbsCurve&, double, double, LNLib::LN_NurbsCurve&>(&LNLib::NurbsCurve::Reparametrize))
        .def_static("Reparametrize", py::overload_cast<const LNLib::LN_NurbsCurve&, double, double, double, double, LNLib::LN_NurbsCurve&>(&LNLib::NurbsCurve::Reparametrize))
        .def_static("Reverse", &LNLib::NurbsCurve::Reverse)
        .def_static("SplitAt", &LNLib::NurbsCurve::SplitAt)
        .def_static("Segment", &LNLib::NurbsCurve::Segment)
        .def_static("Merge", &LNLib::NurbsCurve::Merge)
        .def_static("Offset", &LNLib::NurbsCurve::Offset)
        .def_static("CreateLine", &LNLib::NurbsCurve::CreateLine)
        .def_static("CreateCubicHermite", &LNLib::NurbsCurve::CreateCubicHermite)
        .def_static("CreateArc", &LNLib::NurbsCurve::CreateArc)
        .def_static("CreateOneConicArc", &LNLib::NurbsCurve::CreateOneConicArc)
        .def_static("SplitArc", &LNLib::NurbsCurve::SplitArc)
        .def_static("CreateOpenConic", &LNLib::NurbsCurve::CreateOpenConic)
        .def_static("GlobalInterpolation", py::overload_cast<int, const std::vector<LNLib::XYZ>&, LNLib::LN_NurbsCurve&, const std::vector<double>&>(&LNLib::NurbsCurve::GlobalInterpolation))
        .def_static("GlobalInterpolation", py::overload_cast<int, const std::vector<LNLib::XYZ>&, const std::vector<LNLib::XYZ>&, double, LNLib::LN_NurbsCurve&>(&LNLib::NurbsCurve::GlobalInterpolation))
        .def_static("CubicLocalInterpolation", &LNLib::NurbsCurve::CubicLocalInterpolation)
        .def_static("LeastSquaresApproximation", &LNLib::NurbsCurve::LeastSquaresApproximation)
        .def_static("WeightedAndContrainedLeastSquaresApproximation", &LNLib::NurbsCurve::WeightedAndContrainedLeastSquaresApproximation)
        .def_static("ComputerRemoveKnotErrorBound", &LNLib::NurbsCurve::ComputerRemoveKnotErrorBound)
        .def_static("RemoveKnotsByGivenBound", &LNLib::NurbsCurve::RemoveKnotsByGivenBound)
        .def_static("GlobalApproximationByErrorBound", &LNLib::NurbsCurve::GlobalApproximationByErrorBound)
        .def_static("FitWithConic", &LNLib::NurbsCurve::FitWithConic)
        .def_static("FitWithCubic", &LNLib::NurbsCurve::FitWithCubic)
        .def_static("Normal", &LNLib::NurbsCurve::Normal)
        .def_static("ProjectNormal", &LNLib::NurbsCurve::ProjectNormal)
        .def_static("ControlPointReposition", &LNLib::NurbsCurve::ControlPointReposition)
        .def_static("WeightModification", &LNLib::NurbsCurve::WeightModification)
        .def_static("NeighborWeightsModification", &LNLib::NurbsCurve::NeighborWeightsModification)
        .def_static("Warping", &LNLib::NurbsCurve::Warping)
        .def_static("Flattening", &LNLib::NurbsCurve::Flattening)
        .def_static("Bending", &LNLib::NurbsCurve::Bending)
        .def_static("ConstraintBasedModification", &LNLib::NurbsCurve::ConstraintBasedModification)
        .def_static("IsClamp", &LNLib::NurbsCurve::IsClamp)
        .def_static("ToClampCurve", &LNLib::NurbsCurve::ToClampCurve)
        .def_static("IsPeriodic", &LNLib::NurbsCurve::IsPeriodic)
        .def_static("ToUnclampCurve", &LNLib::NurbsCurve::ToUnclampCurve)
        .def_static("IsLinear", &LNLib::NurbsCurve::IsLinear)
        .def_static("IsArc", &LNLib::NurbsCurve::IsArc)
        .def_static("ApproximateLength", &LNLib::NurbsCurve::ApproximateLength)
        .def_static("GetParamOnCurve", py::overload_cast<const LNLib::LN_NurbsCurve&, double, LNLib::IntegratorType>(&LNLib::NurbsCurve::GetParamOnCurve))
        .def_static("GetParamsOnCurve", &LNLib::NurbsCurve::GetParamsOnCurve)
        .def_static("Tessellate", &LNLib::NurbsCurve::Tessellate);
}
