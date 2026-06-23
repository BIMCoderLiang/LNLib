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
#include <pybind11/stl.h> 
#include "LNObject.h"
#include "XYZ.h"
#include "Matrix4d.h"
#include "NurbsCurve.h"

namespace py = pybind11;
using namespace LNLib;

void cstrNurbsCurve(py::module_&m)
{
    py::class_<LNLib::NurbsCurve>(m, "NurbsCurve")
        .def_static("Check", &LNLib::NurbsCurve::Check)
        .def_static("GetPointOnCurve", &LNLib::NurbsCurve::GetPointOnCurve)
        .def_static("ComputeRationalCurveDerivatives", &LNLib::NurbsCurve::ComputeRationalCurveDerivatives)
        .def_static("CanComputeDerivative", &LNLib::NurbsCurve::CanComputeDerivative)
        .def_static("Curvature", &LNLib::NurbsCurve::Curvature)
        .def_static("Torsion", &LNLib::NurbsCurve::Torsion)
        .def_static("InsertKnot", [](const LNLib::LN_NurbsCurve& curve, double insertKnot, int times) {
        LNLib::LN_NurbsCurve result;
        int ret = LNLib::NurbsCurve::InsertKnot(curve, insertKnot, times, result);
        return py::make_tuple(ret, result);
            })
        .def_static("GetPointOnCurveByCornerCut", &LNLib::NurbsCurve::GetPointOnCurveByCornerCut)
        .def_static("RefineKnotVector", [](const LNLib::LN_NurbsCurve& curve, const std::vector<double>& insertKnotElements) {
        LNLib::LN_NurbsCurve result;
        LNLib::NurbsCurve::RefineKnotVector(curve, insertKnotElements, result);
        return result;
            })
        .def_static("DecomposeToBeziers", &LNLib::NurbsCurve::DecomposeToBeziers)
        .def_static("RemoveKnot", [](const LNLib::LN_NurbsCurve& curve, double removeKnot, int times) {
        LNLib::LN_NurbsCurve result;
        bool ret = LNLib::NurbsCurve::RemoveKnot(curve, removeKnot, times, result);
        return py::make_tuple(ret, result);
            })
        .def_static("RemoveExcessiveKnots", [](const LNLib::LN_NurbsCurve& curve) {
        LNLib::LN_NurbsCurve result;
        LNLib::NurbsCurve::RemoveExcessiveKnots(curve, result);
        return result;
            })
        .def_static("ElevateDegree", [](const LNLib::LN_NurbsCurve& curve, int times) {
        LNLib::LN_NurbsCurve result;
        LNLib::NurbsCurve::ElevateDegree(curve, times, result);
        return result;
            })
        .def_static("ReduceDegree", [](const LNLib::LN_NurbsCurve& curve) {
        LNLib::LN_NurbsCurve result;
        bool ret = LNLib::NurbsCurve::ReduceDegree(curve, result);
        return py::make_tuple(ret, result);
            })
        .def_static("EquallyTessellate", [](const LNLib::LN_NurbsCurve& curve) {
        std::vector<LNLib::XYZ> tessellatedPoints;
        std::vector<double> correspondingKnots;
        LNLib::NurbsCurve::EquallyTessellate(curve, tessellatedPoints, correspondingKnots);
        return py::make_tuple(tessellatedPoints, correspondingKnots);
            })
        .def_static("IsClosed", &LNLib::NurbsCurve::IsClosed)
        .def_static("GetParamOnCurve", py::overload_cast<const LNLib::LN_NurbsCurve&, const LNLib::XYZ&>(&LNLib::NurbsCurve::GetParamOnCurve))
        .def_static("CreateTransformed", [](const LNLib::LN_NurbsCurve& curve, const LNLib::Matrix4d& matrix) {
        LNLib::LN_NurbsCurve result;
        LNLib::NurbsCurve::CreateTransformed(curve, matrix, result);
        return result;
            })
        .def_static("Reparametrize", [](const LNLib::LN_NurbsCurve& curve, double min, double max) {
        LNLib::LN_NurbsCurve result;
        LNLib::NurbsCurve::Reparametrize(curve, min, max, result);
        return result;
            })
        .def_static("Reparametrize", [](const LNLib::LN_NurbsCurve& curve, double alpha, double beta, double gamma, double delta) {
        LNLib::LN_NurbsCurve result;
        LNLib::NurbsCurve::Reparametrize(curve, alpha, beta, gamma, delta, result);
        return result;
            })
        .def_static("Reverse", [](const LNLib::LN_NurbsCurve& curve) {
        LNLib::LN_NurbsCurve result;
        LNLib::NurbsCurve::Reverse(curve, result);
        return result;
            })
        .def_static("SplitAt", [](const LNLib::LN_NurbsCurve& curve, double parameter) {
        LNLib::LN_NurbsCurve left, right;
        bool ret = LNLib::NurbsCurve::SplitAt(curve, parameter, left, right);
        return py::make_tuple(ret, left, right);
            })
        .def_static("Extract", [](const LNLib::LN_NurbsCurve& curve, double startParameter, double endParameter) {
        LNLib::LN_NurbsCurve result;
        LNLib::NurbsCurve::Extract(curve, startParameter, endParameter, result);
        return result;
            })
        .def_static("Merge", [](const LNLib::LN_NurbsCurve& left, const LNLib::LN_NurbsCurve& right) {
        LNLib::LN_NurbsCurve result;
        bool ret = LNLib::NurbsCurve::Merge(left, right, result);
        return py::make_tuple(ret, result);
            })
        .def_static("Offset", [](const LNLib::LN_NurbsCurve& curve, double offset, LNLib::OffsetType type) {
        LNLib::LN_NurbsCurve result;
        LNLib::NurbsCurve::Offset(curve, offset, type, result);
        return result;
            })
        .def_static("CreateLine", [](const LNLib::XYZ& start, const LNLib::XYZ& end) {
        LNLib::LN_NurbsCurve result;
        LNLib::NurbsCurve::CreateLine(start, end, result);
        return result;
            })
        .def_static("CreateCubicHermite", [](const std::vector<LNLib::XYZ>& throughPoints, const std::vector<LNLib::XYZ>& tangents) {
        LNLib::LN_NurbsCurve result;
        LNLib::NurbsCurve::CreateCubicHermite(throughPoints, tangents, result);
        return result;
            })
        .def_static("CreateArc", [](const LNLib::XYZ& center, const LNLib::XYZ& xAxis, const LNLib::XYZ& yAxis, double startRad, double endRad, double xRadius, double yRadius) {
        LNLib::LN_NurbsCurve result;
        bool ret = LNLib::NurbsCurve::CreateArc(center, xAxis, yAxis, startRad, endRad, xRadius, yRadius, result);
        return py::make_tuple(ret, result);
            })
        .def_static("CreateOneConicArc", [](const LNLib::XYZ& start, const LNLib::XYZ& startTangent, const LNLib::XYZ& end, const LNLib::XYZ& endTangent, const LNLib::XYZ& pointOnConic) {
        LNLib::XYZ projectPoint;
        double projectPointWeight;
        bool ret = LNLib::NurbsCurve::CreateOneConicArc(start, startTangent, end, endTangent, pointOnConic, projectPoint, projectPointWeight);
        return py::make_tuple(ret, projectPoint, projectPointWeight);
            })
        .def_static("SplitArc", [](const LNLib::XYZ& start, const LNLib::XYZ& projectPoint, double projectPointWeight, const LNLib::XYZ& end, double insertWeight) {
        LNLib::XYZ insertPointAtStartSide, splitPoint, insertPointAtEndSide;
        LNLib::NurbsCurve::SplitArc(start, projectPoint, projectPointWeight, end, insertPointAtStartSide, splitPoint, insertPointAtEndSide, insertWeight);
        return py::make_tuple(insertPointAtStartSide, splitPoint, insertPointAtEndSide);
            })
        .def_static("CreateOpenConic", [](const LNLib::XYZ& start, const LNLib::XYZ& startTangent, const LNLib::XYZ& end, const LNLib::XYZ& endTangent, const LNLib::XYZ& pointOnConic) {
        LNLib::LN_NurbsCurve result;
        bool ret = LNLib::NurbsCurve::CreateOpenConic(start, startTangent, end, endTangent, pointOnConic, result);
        return py::make_tuple(ret, result);
            })
        .def_static("GlobalInterpolation", [](int degree, const std::vector<LNLib::XYZ>& throughPoints, const std::vector<double>& params) {
        LNLib::LN_NurbsCurve result;
        LNLib::NurbsCurve::GlobalInterpolation(degree, throughPoints, result, params);
        return result;
            }, py::arg("degree"), py::arg("throughPoints"), py::arg("params") = std::vector<double>{})
        .def_static("GlobalInterpolation", [](int degree, const std::vector<LNLib::XYZ>& throughPoints, const std::vector<LNLib::XYZ>& tangents, double tangentFactor) {
        LNLib::LN_NurbsCurve result;
        LNLib::NurbsCurve::GlobalInterpolation(degree, throughPoints, tangents, tangentFactor, result);
        return result;
            })
        .def_static("CubicLocalInterpolation", [](const std::vector<LNLib::XYZ>& throughPoints) {
        LNLib::LN_NurbsCurve result;
        bool ret = LNLib::NurbsCurve::CubicLocalInterpolation(throughPoints, result);
        return py::make_tuple(ret, result);
            })
        .def_static("LeastSquaresApproximation", [](int degree, const std::vector<LNLib::XYZ>& throughPoints, int controlPointsCount) {
        LNLib::LN_NurbsCurve result;
        bool ret = LNLib::NurbsCurve::LeastSquaresApproximation(degree, throughPoints, controlPointsCount, result);
        return py::make_tuple(ret, result);
            })
        .def_static("WeightedAndContrainedLeastSquaresApproximation", [](int degree, const std::vector<LNLib::XYZ>& throughPoints, const std::vector<double>& throughPointWeights, const std::vector<LNLib::XYZ>& tangents, const std::vector<int>& tangentIndices, const std::vector<double>& tangentWeights, int controlPointsCount) {
        LNLib::LN_NurbsCurve result;
        bool ret = LNLib::NurbsCurve::WeightedAndContrainedLeastSquaresApproximation(degree, throughPoints, throughPointWeights, tangents, tangentIndices, tangentWeights, controlPointsCount, result);
        return py::make_tuple(ret, result);
            })
        .def_static("ComputerRemoveKnotErrorBound", &LNLib::NurbsCurve::ComputerRemoveKnotErrorBound)
        .def_static("RemoveKnotsByGivenBound", [](const LNLib::LN_NurbsCurve& curve, const std::vector<double>& params, double maxError) {
        LNLib::LN_NurbsCurve result;
        std::vector<double> errors;
        LNLib::NurbsCurve::RemoveKnotsByGivenBound(curve, params, errors, maxError, result);
        return py::make_tuple(errors, result);
            })
        .def_static("GlobalApproximationByErrorBound", [](int degree, const std::vector<LNLib::XYZ>& throughPoints, double maxError) {
        LNLib::LN_NurbsCurve result;
        LNLib::NurbsCurve::GlobalApproximationByErrorBound(degree, throughPoints, maxError, result);
        return result;
            })
        .def_static("FitWithConic", [](const std::vector<LNLib::XYZ>& throughPoints, int startPointIndex, int endPointIndex, const LNLib::XYZ& startTangent, const LNLib::XYZ& endTangent, double maxError) {
        std::vector<LNLib::XYZW> middleControlPoints;
        bool ret = LNLib::NurbsCurve::FitWithConic(throughPoints, startPointIndex, endPointIndex, startTangent, endTangent, maxError, middleControlPoints);
        return py::make_tuple(ret, middleControlPoints);
            })
        .def_static("FitWithCubic", [](const std::vector<LNLib::XYZ>& throughPoints, int startPointIndex, int endPointIndex, const LNLib::XYZ& startTangent, const LNLib::XYZ& endTangent, double maxError) {
        std::vector<LNLib::XYZW> middleControlPoints;
        bool ret = LNLib::NurbsCurve::FitWithCubic(throughPoints, startPointIndex, endPointIndex, startTangent, endTangent, maxError, middleControlPoints);
        return py::make_tuple(ret, middleControlPoints);
            })
        .def_static("Normal", &LNLib::NurbsCurve::Normal)
        .def_static("ProjectNormal", &LNLib::NurbsCurve::ProjectNormal)
        .def_static("ControlPointReposition", [](const LNLib::LN_NurbsCurve& curve, double parameter, int moveIndex, LNLib::XYZ moveDirection, double moveDistance) {
        LNLib::LN_NurbsCurve result;
        bool ret = LNLib::NurbsCurve::ControlPointReposition(curve, parameter, moveIndex, moveDirection, moveDistance, result);
        return py::make_tuple(ret, result);
            })
        .def_static("WeightModification", [](const LNLib::LN_NurbsCurve& curve, double parameter, int moveIndex, double moveDistance) {
        LNLib::LN_NurbsCurve result;
        LNLib::NurbsCurve::WeightModification(curve, parameter, moveIndex, moveDistance, result);
        return result;
            })
        .def_static("NeighborWeightsModification", [](const LNLib::LN_NurbsCurve& curve, double parameter, int moveIndex, double moveDistance, double scale) {
        LNLib::LN_NurbsCurve result;
        bool ret = LNLib::NurbsCurve::NeighborWeightsModification(curve, parameter, moveIndex, moveDistance, scale, result);
        return py::make_tuple(ret, result);
            })
        .def_static("Warping", [](const LNLib::LN_NurbsCurve& curve, const std::vector<double>& warpShape, double warpDistance, const LNLib::XYZ& planeNormal, double startParameter, double endParameter) {
        LNLib::LN_NurbsCurve result;
        LNLib::NurbsCurve::Warping(curve, warpShape, warpDistance, planeNormal, startParameter, endParameter, result);
        return result;
            })
        .def_static("Flattening", [](const LNLib::LN_NurbsCurve& curve, LNLib::XYZ lineStartPoint, LNLib::XYZ lineEndPoint, double startParameter, double endParameter) {
        LNLib::LN_NurbsCurve result;
        bool ret = LNLib::NurbsCurve::Flattening(curve, lineStartPoint, lineEndPoint, startParameter, endParameter, result);
        return py::make_tuple(ret, result);
            })
        .def_static("Bending", [](const LNLib::LN_NurbsCurve& curve, double startParameter, double endParameter, LNLib::XYZ bendCenter, double radius, double crossRatio) {
        LNLib::LN_NurbsCurve result;
        LNLib::NurbsCurve::Bending(curve, startParameter, endParameter, bendCenter, radius, crossRatio, result);
        return result;
            })
        .def_static("ConstraintBasedModification", [](const LNLib::LN_NurbsCurve& curve, const std::vector<double>& constraintParams, const std::vector<LNLib::XYZ>& derivativeConstraints, const std::vector<int>& appliedIndices, const std::vector<int>& appliedDegree, const std::vector<int>& fixedControlPointIndices) {
        LNLib::LN_NurbsCurve result;
        LNLib::NurbsCurve::ConstraintBasedModification(curve, constraintParams, derivativeConstraints, appliedIndices, appliedDegree, fixedControlPointIndices, result);
        return result;
            })
        .def_static("IsClamp", &LNLib::NurbsCurve::IsClamp)
        .def_static("ToClampCurve", [](const LNLib::LN_NurbsCurve& curve) {
        LNLib::LN_NurbsCurve result;
        LNLib::NurbsCurve::ToClampCurve(curve, result);
        return result;
            })
        .def_static("IsPeriodic", &LNLib::NurbsCurve::IsPeriodic)
        .def_static("ToUnclampCurve", [](const LNLib::LN_NurbsCurve& curve) {
        LNLib::LN_NurbsCurve result;
        LNLib::NurbsCurve::ToUnclampCurve(curve, result);
        return result;
            })
        .def_static("IsLinear", &LNLib::NurbsCurve::IsLinear)
        .def_static("IsArc", [](const LNLib::LN_NurbsCurve& curve) {
        LNLib::LN_ArcInfo arcInfo;
        bool ret = LNLib::NurbsCurve::IsArc(curve, arcInfo);
        return py::make_tuple(ret, arcInfo);
            })
        .def_static("ApproximateLength", py::overload_cast<const LNLib::LN_NurbsCurve&, LNLib::IntegratorType>(&LNLib::NurbsCurve::ApproximateLength), py::arg("curve"), py::arg("type") = LNLib::IntegratorType::GaussLegendre)
        .def_static("ApproximateLength", py::overload_cast<const LNLib::LN_NurbsCurve&, double, double, LNLib::IntegratorType>(&LNLib::NurbsCurve::ApproximateLength), py::arg("curve"), py::arg("startParameter"), py::arg("endParameter"), py::arg("type") = LNLib::IntegratorType::GaussLegendre)
        .def_static("Extend", [](const LNLib::LN_NurbsCurve& curve, double delta, bool isFromStart, LNLib::ExtensionType type) {
        LNLib::LN_NurbsCurve result;
        LNLib::NurbsCurve::Extend(curve, delta, isFromStart, type, result);
        return result;
            })
        .def_static("GetParamOnCurve", py::overload_cast<const LNLib::LN_NurbsCurve&, double, LNLib::IntegratorType>(&LNLib::NurbsCurve::GetParamOnCurve), py::arg("curve"), py::arg("givenLength"), py::arg("type") = LNLib::IntegratorType::GaussLegendre)
        .def_static("GetParamsOnCurve", &LNLib::NurbsCurve::GetParamsOnCurve, py::arg("curve"), py::arg("givenLength"), py::arg("type") = LNLib::IntegratorType::GaussLegendre)
        .def_static("Tessellate", &LNLib::NurbsCurve::Tessellate);
}
