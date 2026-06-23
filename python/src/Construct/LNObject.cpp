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
#include "XYZW.h"
#include "UV.h"

namespace py = pybind11;
void cstrObject(py::module_&m)
{
    py::class_<LNLib::LN_BezierCurve<LNLib::XYZ>>(m, "LN_BezierCurve_XYZ")
        .def(py::init<>())
        .def_readwrite("Degree", &LNLib::LN_BezierCurve<LNLib::XYZ>::Degree)
        .def_readwrite("ControlPoints", &LNLib::LN_BezierCurve<LNLib::XYZ>::ControlPoints);

    py::class_<LNLib::LN_BezierCurve<LNLib::XYZW>>(m, "LN_BezierCurve_XYZW")
        .def(py::init<>())
        .def_readwrite("Degree", &LNLib::LN_BezierCurve<LNLib::XYZW>::Degree)
        .def_readwrite("ControlPoints", &LNLib::LN_BezierCurve<LNLib::XYZW>::ControlPoints);

    py::class_<LNLib::LN_BsplineCurve<LNLib::XYZ>>(m, "LN_BsplineCurve_XYZ")
        .def(py::init<>())
        .def_readwrite("Degree", &LNLib::LN_BsplineCurve<LNLib::XYZ>::Degree)
        .def_readwrite("KnotVector", &LNLib::LN_BsplineCurve<LNLib::XYZ>::KnotVector)
        .def_readwrite("ControlPoints", &LNLib::LN_BsplineCurve<LNLib::XYZ>::ControlPoints);

    py::class_<LNLib::LN_BsplineCurve<LNLib::XYZW>>(m, "LN_NurbsCurve", py::module_local())
        .def(py::init<>())
        .def_readwrite("Degree", &LNLib::LN_BsplineCurve<LNLib::XYZW>::Degree)
        .def_readwrite("KnotVector", &LNLib::LN_BsplineCurve<LNLib::XYZW>::KnotVector)
        .def_readwrite("ControlPoints", &LNLib::LN_BsplineCurve<LNLib::XYZW>::ControlPoints);

    py::class_<LNLib::LN_Mesh>(m, "LN_Mesh")
        .def(py::init<>())
        .def_readwrite("Vertices", &LNLib::LN_Mesh::Vertices)
        .def_readwrite("Faces", &LNLib::LN_Mesh::Faces)
        .def_readwrite("UVs", &LNLib::LN_Mesh::UVs)
        .def_readwrite("UVIndices", &LNLib::LN_Mesh::UVIndices)
        .def_readwrite("Normals", &LNLib::LN_Mesh::Normals)
        .def_readwrite("NormalIndices", &LNLib::LN_Mesh::NormalIndices);

    py::class_<LNLib::LN_ArcInfo>(m, "LN_ArcInfo")
        .def(py::init<>())
        .def_readwrite("Radius", &LNLib::LN_ArcInfo::Radius)
        .def_readwrite("Center", &LNLib::LN_ArcInfo::Center);

    py::class_<LNLib::LN_BoundingBox3d>(m, "LN_BoundingBox3d")
        .def(py::init<>())
        .def_readwrite("MinPoint", &LNLib::LN_BoundingBox3d::MinPoint)
        .def_readwrite("MaxPoint", &LNLib::LN_BoundingBox3d::MaxPoint)
        .def("Intersects", &LNLib::LN_BoundingBox3d::Intersects);

    py::class_<LNLib::LN_OrientedBoundingBox3d>(m, "LN_OrientedBoundingBox3d")
        .def(py::init<>())
        .def_readwrite("Center", &LNLib::LN_OrientedBoundingBox3d::Center)
        .def_property("Axes",
            [](const LNLib::LN_OrientedBoundingBox3d& self) {
                return std::vector<LNLib::XYZ>{self.Axes[0], self.Axes[1], self.Axes[2]};
            },
            [](LNLib::LN_OrientedBoundingBox3d& self, const std::vector<LNLib::XYZ>& val) {
                if (val.size() != 3) throw std::runtime_error("Axes size must be 3");
                self.Axes[0] = val[0];
                self.Axes[1] = val[1];
                self.Axes[2] = val[2];
            })
        .def_property("HalfExtents",
            [](const LNLib::LN_OrientedBoundingBox3d& self) {
                return std::vector<double>{self.HalfExtents[0], self.HalfExtents[1], self.HalfExtents[2]};
            },
            [](LNLib::LN_OrientedBoundingBox3d& self, const std::vector<double>& val) {
                if (val.size() != 3) throw std::runtime_error("HalfExtents size must be 3");
                self.HalfExtents[0] = val[0];
                self.HalfExtents[1] = val[1];
                self.HalfExtents[2] = val[2];
            });
}
