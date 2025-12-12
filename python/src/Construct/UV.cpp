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
#include "Constants.h"
#include "UV.h"

namespace py = pybind11;

void cstrUV(py::module_&m)
{
    py::class_<LNLib::UV>(m, "UV")
        .def(py::init<>())
        .def(py::init<double, double>())
        .def("SetU", &LNLib::UV::SetU)
        .def("GetU", &LNLib::UV::GetU)
        .def("SetV", &LNLib::UV::SetV)
        .def("GetV", &LNLib::UV::GetV)
        .def("U", py::overload_cast<>(&LNLib::UV::U, py::const_))
        .def("U", py::overload_cast<>(&LNLib::UV::U))
        .def("V", py::overload_cast<>(&LNLib::UV::V, py::const_))
        .def("V", py::overload_cast<>(&LNLib::UV::V))
        .def("IsZero", &LNLib::UV::IsZero, py::arg("epsilon") = LNLib::Constants::DoubleEpsilon)
        .def("IsUnit", &LNLib::UV::IsUnit, py::arg("epsilon") = LNLib::Constants::DoubleEpsilon)
        .def("IsAlmostEqualTo", &LNLib::UV::IsAlmostEqualTo)
        .def("Length", &LNLib::UV::Length)
        .def("SqrLength", &LNLib::UV::SqrLength)
        .def("Normalize", &LNLib::UV::Normalize)
        .def("Add", &LNLib::UV::Add)
        .def("Substract", &LNLib::UV::Substract)
        .def("Negative", &LNLib::UV::Negative)
        .def("DotProduct", &LNLib::UV::DotProduct)
        .def("CrossProduct", &LNLib::UV::CrossProduct)
        .def("Distance", &LNLib::UV::Distance)
        .def("__getitem__", [](const LNLib::UV& uv, int index) {
        if (index < 0 || index >= 2) throw py::index_error();
        return uv[index];
            })
        .def("__setitem__", [](LNLib::UV& uv, int index, double value) {
        if (index < 0 || index >= 2) throw py::index_error();
        uv[index] = value;
            })
        .def("__add__", [](const LNLib::UV& uv1, const LNLib::UV& uv2) {
        return uv1 + uv2;
            })
        .def("__sub__", [](const LNLib::UV& uv1, const LNLib::UV& uv2) {
        return uv1 - uv2;
            })
        .def("__mul__", [](const LNLib::UV& uv1, const LNLib::UV& uv2) {
        return uv1 * uv2;
            })
        .def("__mul__", [](const LNLib::UV& uv, double d) {
        return uv * d;
            })
        .def("__rmul__", [](const LNLib::UV& uv, double d) {
        return uv * d;
            })
        .def("__truediv__", [](const LNLib::UV& uv, double d) {
        return uv / d;
            })
        .def("__neg__", [](const LNLib::UV& uv) {
        return -uv;
            })
        .def("__repr__", [](const LNLib::UV& uv) {
        return "UV(" + std::to_string(uv.GetU()) + ", " + std::to_string(uv.GetV()) + ")";
            });	
}



