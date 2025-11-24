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
#include "XYZW.h"
#include "LNObject.h"

namespace py = pybind11;
void cstrObject(py::module_&m)
{
    py::class_<LNLib::LN_BezierCurve<LNLib::XYZW>>(m, "BezierCurve")
        .def(py::init<>())
        .def_readwrite("Degree", &LNLib::LN_BezierCurve<LNLib::XYZW>::Degree)
        .def_readwrite("ControlPoints", &LNLib::LN_BezierCurve<LNLib::XYZW>::ControlPoints);

    py::class_<LNLib::LN_BezierSurface<LNLib::XYZW>>(m, "BezierSurface")
        .def(py::init<>())
        .def_readwrite("DegreeU", &LNLib::LN_BezierSurface<LNLib::XYZW>::DegreeU)
        .def_readwrite("DegreeV", &LNLib::LN_BezierSurface<LNLib::XYZW>::DegreeV)
        .def_property("ControlPoints",
            [](const LNLib::LN_BsplineSurface<LNLib::XYZW>& surface) {
                return surface.ControlPoints;
            },
            [](LNLib::LN_BsplineSurface<LNLib::XYZW>& surface, const std::vector<std::vector<LNLib::XYZW>>& controlPoints) {
                surface.ControlPoints = controlPoints;
            });

    py::class_<LNLib::LN_NurbsCurve>(m, "NurbsCurve")
        .def(py::init<>())
        .def_readwrite("Degree", &LNLib::LN_NurbsCurve::Degree)
        .def_readwrite("KnotVector", &LNLib::LN_NurbsCurve::KnotVector)
        .def_readwrite("ControlPoints", &LNLib::LN_NurbsCurve::ControlPoints);

    py::class_<LNLib::LN_NurbsSurface>(m, "NurbsSurface")
        .def(py::init<>())
        .def_readwrite("DegreeU", &LNLib::LN_NurbsSurface::DegreeU)
        .def_readwrite("DegreeV", &LNLib::LN_NurbsSurface::DegreeV)
        .def_readwrite("KnotVectorU", &LNLib::LN_NurbsSurface::KnotVectorU)
        .def_readwrite("KnotVectorV", &LNLib::LN_NurbsSurface::KnotVectorV)
        .def_property("ControlPoints",
            [](const LNLib::LN_NurbsSurface& surface) {
                return surface.ControlPoints;
            },
            [](LNLib::LN_NurbsSurface& surface, const std::vector<std::vector<LNLib::XYZW>>& controlPoints) {
                surface.ControlPoints = controlPoints;
            });

    py::class_<LNLib::LN_Mesh>(m, "Mesh")
        .def(py::init<>())
        .def_readwrite("Vertices", &LNLib::LN_Mesh::Vertices)
        .def_readwrite("Faces", &LNLib::LN_Mesh::Faces)
        .def_readwrite("UVs", &LNLib::LN_Mesh::UVs)
        .def_readwrite("UVIndices", &LNLib::LN_Mesh::UVIndices)
        .def_readwrite("Normals", &LNLib::LN_Mesh::Normals)
        .def_readwrite("NormalIndices", &LNLib::LN_Mesh::NormalIndices);

    py::class_<LNLib::LN_ArcInfo>(m, "ArcInfo")
        .def(py::init<>())
        .def_readwrite("Radius", &LNLib::LN_ArcInfo::Radius)
        .def_readwrite("Center", &LNLib::LN_ArcInfo::Center);
}
