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
#include "XYZW.h"
#include "Matrix4d.h"

namespace py = pybind11;
void cstrMatrix4d(py::module_&m)
{
    py::class_<LNLib::Matrix4d>(m, "Matrix4d")
        .def(py::init<>())
        .def(py::init<LNLib::XYZ, LNLib::XYZ, LNLib::XYZ, LNLib::XYZ>())
        .def(py::init<double, double, double, double,
            double, double, double, double,
            double, double, double, double,
            double, double, double, double>())
        .def_static("CreateReflection", py::overload_cast<const LNLib::XYZ&, double>(&LNLib::Matrix4d::CreateReflection))
        .def_static("CreateReflection", py::overload_cast<const LNLib::XYZ&>(&LNLib::Matrix4d::CreateReflection))
        .def_static("CreateRotation", &LNLib::Matrix4d::CreateRotation)
        .def_static("CreateRotationAtPoint", &LNLib::Matrix4d::CreateRotationAtPoint)
        .def_static("CreateTranslation", &LNLib::Matrix4d::CreateTranslation)
        .def_static("CreateScale", &LNLib::Matrix4d::CreateScale)
        .def("SetBasisX", &LNLib::Matrix4d::SetBasisX)
        .def("GetBasisX", &LNLib::Matrix4d::GetBasisX)
        .def("SetBasisY", &LNLib::Matrix4d::SetBasisY)
        .def("GetBasisY", &LNLib::Matrix4d::GetBasisY)
        .def("SetBasisZ", &LNLib::Matrix4d::SetBasisZ)
        .def("GetBasisZ", &LNLib::Matrix4d::GetBasisZ)
        .def("SetBasisW", &LNLib::Matrix4d::SetBasisW)
        .def("GetBasisW", &LNLib::Matrix4d::GetBasisW)
        .def("GetElement", &LNLib::Matrix4d::GetElement)
        .def("SetElement", &LNLib::Matrix4d::SetElement)
        .def("Multiply", &LNLib::Matrix4d::Multiply)
        .def("OfPoint", &LNLib::Matrix4d::OfPoint)
        .def("OfWeightedPoint", &LNLib::Matrix4d::OfWeightedPoint)
        .def("OfVector", &LNLib::Matrix4d::OfVector)
        .def("GetInverse", &LNLib::Matrix4d::GetInverse)
        .def("GetTranspose", &LNLib::Matrix4d::GetTranspose)
        .def("GetScale", &LNLib::Matrix4d::GetScale)
        .def("GetDeterminant", &LNLib::Matrix4d::GetDeterminant)
        .def("IsIdentity", &LNLib::Matrix4d::IsIdentity)
        .def("HasReflection", &LNLib::Matrix4d::HasReflection)
        .def("IsTranslation", &LNLib::Matrix4d::IsTranslation)
        .def("__mul__", [](const LNLib::Matrix4d& m1, const LNLib::Matrix4d& m2) {
        return m1 * m2;
            })
        .def("__add__", [](const LNLib::Matrix4d& m1, const LNLib::Matrix4d& m2) {
        return m1 + m2;
            })
        .def("__sub__", [](const LNLib::Matrix4d& m1, const LNLib::Matrix4d& m2) {
        return m1 - m2;
            });
}
