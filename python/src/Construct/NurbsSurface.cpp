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
        .def_static("Check", &LNLib::NurbsSurface::Check)
        .def_static("GetIsoCurve", &LNLib::NurbsSurface::GetIsoCurve)
        .def_static("GetPointOnSurface", &LNLib::NurbsSurface::GetPointOnSurface)
        .def_static("ComputeRationalSurfaceDerivatives", &LNLib::NurbsSurface::ComputeRationalSurfaceDerivatives)
        .def_static("ComputeRationalSurfaceFirstOrderDerivative", [](const LNLib::LN_NurbsSurface& surface, LNLib::UV uv) {
        LNLib::XYZ S, Su, Sv;
        LNLib::NurbsSurface::ComputeRationalSurfaceFirstOrderDerivative(surface, uv, S, Su, Sv);
        return py::make_tuple(S, Su, Sv);
            })
        .def_static("Curvature", &LNLib::NurbsSurface::Curvature)
        .def_static("GetBoundingBox", &LNLib::NurbsSurface::GetBoundingBox)
        .def_static("GetOrientedBoundingBox", &LNLib::NurbsSurface::GetOrientedBoundingBox)
        .def_static("Normal", &LNLib::NurbsSurface::Normal)
        .def_static("Swap", [](const LNLib::LN_NurbsSurface& surface) {
        LNLib::LN_NurbsSurface result;
        LNLib::NurbsSurface::Swap(surface, result);
        return result;
            })
        .def_static("Reverse", [](const LNLib::LN_NurbsSurface& surface, LNLib::SurfaceDirection direction) {
        LNLib::LN_NurbsSurface result;
        LNLib::NurbsSurface::Reverse(surface, direction, result);
        return result;
            })
        .def_static("InsertKnot", [](const LNLib::LN_NurbsSurface& surface, double insertKnot, int times, bool isUDirection) {
        LNLib::LN_NurbsSurface result;
        LNLib::NurbsSurface::InsertKnot(surface, insertKnot, times, isUDirection, result);
        return result;
            })
        .def_static("RefineKnotVector", [](const LNLib::LN_NurbsSurface& surface, const std::vector<double>& insertKnotElements, bool isUDirection) {
        LNLib::LN_NurbsSurface result;
        LNLib::NurbsSurface::RefineKnotVector(surface, insertKnotElements, isUDirection, result);
        return result;
            })
        .def_static("DecomposeToBeziers", &LNLib::NurbsSurface::DecomposeToBeziers)
        .def_static("RemoveKnot", [](const LNLib::LN_NurbsSurface& surface, double removeKnot, int times, bool isUDirection) {
        LNLib::LN_NurbsSurface result;
        LNLib::NurbsSurface::RemoveKnot(surface, removeKnot, times, isUDirection, result);
        return result;
            })
        .def_static("ElevateDegree", [](const LNLib::LN_NurbsSurface& surface, int times, bool isUDirection) {
        LNLib::LN_NurbsSurface result;
        LNLib::NurbsSurface::ElevateDegree(surface, times, isUDirection, result);
        return result;
            })
        .def_static("ReduceDegree", [](const LNLib::LN_NurbsSurface& surface, bool isUDirection) {
        LNLib::LN_NurbsSurface result;
        bool ret = LNLib::NurbsSurface::ReduceDegree(surface, isUDirection, result);
        return py::make_tuple(ret, result);
            })
        .def_static("EquallyTessellate", [](const LNLib::LN_NurbsSurface& surface) {
        std::vector<LNLib::XYZ> tessellatedPoints;
        std::vector<LNLib::UV> correspondingKnots;
        LNLib::NurbsSurface::EquallyTessellate(surface, tessellatedPoints, correspondingKnots);
        return py::make_tuple(tessellatedPoints, correspondingKnots);
            })
        .def_static("IsClosed", &LNLib::NurbsSurface::IsClosed)
        .def_static("GetParamOnSurface", &LNLib::NurbsSurface::GetParamOnSurface)
        .def_static("GetParamOnSurfaceByGSA", &LNLib::NurbsSurface::GetParamOnSurfaceByGSA)
        .def_static("Reparametrize", [](const LNLib::LN_NurbsSurface& surface, double minU, double maxU, double minV, double maxV) {
        LNLib::LN_NurbsSurface result;
        LNLib::NurbsSurface::Reparametrize(surface, minU, maxU, minV, maxV, result);
        return result;
            })
        .def_static("GetUVTangent", [](const LNLib::LN_NurbsSurface& surface, const LNLib::UV param, const LNLib::XYZ& tangent) {
        LNLib::UV uvTangent;
        bool ret = LNLib::NurbsSurface::GetUVTangent(surface, param, tangent, uvTangent);
        return py::make_tuple(ret, uvTangent);
            })
        .def_static("CreateBilinearSurface", [](const LNLib::XYZ& topLeftPoint, const LNLib::XYZ& topRightPoint, const LNLib::XYZ& bottomLeftPoint, const LNLib::XYZ& bottomRightPoint) {
        LNLib::LN_NurbsSurface result;
        LNLib::NurbsSurface::CreateBilinearSurface(topLeftPoint, topRightPoint, bottomLeftPoint, bottomRightPoint, result);
        return result;
            })
        .def_static("CreateCylindricalSurface", [](const LNLib::XYZ& origin, const LNLib::XYZ& xAxis, const LNLib::XYZ& yAxis, double startRad, double endRad, double radius, double height) {
        LNLib::LN_NurbsSurface result;
        bool ret = LNLib::NurbsSurface::CreateCylindricalSurface(origin, xAxis, yAxis, startRad, endRad, radius, height, result);
        return py::make_tuple(ret, result);
            })
        .def_static("CreateRuledSurface", [](const LNLib::LN_NurbsCurve& curve0, const LNLib::LN_NurbsCurve& curve1) {
        LNLib::LN_NurbsSurface result;
        LNLib::NurbsSurface::CreateRuledSurface(curve0, curve1, result);
        return result;
            })
        .def_static("CreateRevolvedSurface", [](const LNLib::XYZ& origin, const LNLib::XYZ& axis, double rad, const LNLib::LN_NurbsCurve& profile) {
        LNLib::LN_NurbsSurface result;
        bool ret = LNLib::NurbsSurface::CreateRevolvedSurface(origin, axis, rad, profile, result);
        return py::make_tuple(ret, result);
            })
        .def_static("NonUniformScaling", &LNLib::NurbsSurface::NonUniformScaling)
        .def_static("MakeCornerFilletSurface", [](const LNLib::LN_NurbsCurve& arc) {
        LNLib::LN_NurbsSurface result;
        LNLib::NurbsSurface::MakeCornerFilletSurface(arc, result);
        return result;
            })
        .def_static("GlobalInterpolation", [](const std::vector<std::vector<LNLib::XYZ>>& throughPoints, int degreeU, int degreeV) {
        LNLib::LN_NurbsSurface result;
        LNLib::NurbsSurface::GlobalInterpolation(throughPoints, degreeU, degreeV, result);
        return result;
            })
        .def_static("BicubicLocalInterpolation", [](const std::vector<std::vector<LNLib::XYZ>>& throughPoints) {
        LNLib::LN_NurbsSurface result;
        bool ret = LNLib::NurbsSurface::BicubicLocalInterpolation(throughPoints, result);
        return py::make_tuple(ret, result);
            })
        .def_static("GlobalApproximation", [](const std::vector<std::vector<LNLib::XYZ>>& throughPoints, int degreeU, int degreeV, int controlPointsRows, int controlPointsColumns) {
        LNLib::LN_NurbsSurface result;
        bool ret = LNLib::NurbsSurface::GlobalApproximation(throughPoints, degreeU, degreeV, controlPointsRows, controlPointsColumns, result);
        return py::make_tuple(ret, result);
            })
        .def_static("CreateSwungSurface", [](const LNLib::LN_NurbsCurve& profile, const LNLib::LN_NurbsCurve& trajectory, double scale) {
        LNLib::LN_NurbsSurface result;
        bool ret = LNLib::NurbsSurface::CreateSwungSurface(profile, trajectory, scale, result);
        return py::make_tuple(ret, result);
            })
        .def_static("CreateLoftSurface", [](const std::vector<LNLib::LN_NurbsCurve>& sections, int customTrajectoryDegree, const std::vector<double>& customTrajectoryKnotVector) {
        LNLib::LN_NurbsSurface result;
        LNLib::NurbsSurface::CreateLoftSurface(sections, result, customTrajectoryDegree, customTrajectoryKnotVector);
        return result;
            }, py::arg("sections"), py::arg("customTrajectoryDegree") = 0, py::arg("customTrajectoryKnotVector") = std::vector<double>{})
        .def_static("CreateGeneralizedTranslationalSweepSurface", [](const LNLib::LN_NurbsCurve& profile, const LNLib::LN_NurbsCurve& trajectory) {
        LNLib::LN_NurbsSurface result;
        LNLib::NurbsSurface::CreateGeneralizedTranslationalSweepSurface(profile, trajectory, result);
        return result;
            })
        .def_static("CreateSweepSurface", [](const LNLib::LN_NurbsCurve& profile, const LNLib::LN_NurbsCurve& trajectory, int minimumProfiles) {
        LNLib::LN_NurbsSurface result;
        LNLib::NurbsSurface::CreateSweepSurface(profile, trajectory, minimumProfiles, result);
        return result;
            })
        .def_static("CreateSweepSurface", [](const LNLib::LN_NurbsCurve& profile, const LNLib::LN_NurbsCurve& trajectory, int minimumProfiles, int customTrajectoryDegree) {
        LNLib::LN_NurbsSurface result;
        LNLib::NurbsSurface::CreateSweepSurface(profile, trajectory, minimumProfiles, customTrajectoryDegree, result);
        return result;
            })
        .def_static("CreateGordonSurface", [](const std::vector<LNLib::LN_NurbsCurve>& uCurves, const std::vector<LNLib::LN_NurbsCurve>& vCurves, const std::vector<std::vector<LNLib::XYZ>>& intersectionPoints) {
        LNLib::LN_NurbsSurface result;
        LNLib::NurbsSurface::CreateGordonSurface(uCurves, vCurves, intersectionPoints, result);
        return result;
            })
        .def_static("CreateCoonsSurface", [](const LNLib::LN_NurbsCurve& leftCurve, const LNLib::LN_NurbsCurve& bottomCurve, const LNLib::LN_NurbsCurve& rightCurve, const LNLib::LN_NurbsCurve& topCurve) {
        LNLib::LN_NurbsSurface result;
        LNLib::NurbsSurface::CreateCoonsSurface(leftCurve, bottomCurve, rightCurve, topCurve, result);
        return result;
            })
        .def_static("ApproximateArea", &LNLib::NurbsSurface::ApproximateArea, py::arg("surface"), py::arg("type") = LNLib::IntegratorType::GaussLegendre)
        .def_static("Triangulate", &LNLib::NurbsSurface::Triangulate);
}
