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

namespace py = pybind11;
void cstrXYZW(py::module_&m)
{
    py::class_<LNLib::XYZW>(m, "XYZW")
        .def(py::init<>())
        .def(py::init<LNLib::XYZ, double>())
        .def(py::init<double, double, double, double>())
        .def("GetWX", &LNLib::XYZW::GetWX)
        .def("GetWY", &LNLib::XYZW::GetWY)
        .def("GetWZ", &LNLib::XYZW::GetWZ)
        .def("SetW", &LNLib::XYZW::SetW)
        .def("GetW", &LNLib::XYZW::GetW)
        .def("WX", py::overload_cast<>(&LNLib::XYZW::WX, py::const_))
        .def("WX", py::overload_cast<>(&LNLib::XYZW::WX))
        .def("WY", py::overload_cast<>(&LNLib::XYZW::WY, py::const_))
        .def("WY", py::overload_cast<>(&LNLib::XYZW::WY))
        .def("WZ", py::overload_cast<>(&LNLib::XYZW::WZ, py::const_))
        .def("WZ", py::overload_cast<>(&LNLib::XYZW::WZ))
        .def("W", py::overload_cast<>(&LNLib::XYZW::W, py::const_))
        .def("W", py::overload_cast<>(&LNLib::XYZW::W))
        .def("ToXYZ", &LNLib::XYZW::ToXYZ)
        .def("IsAlmostEqualTo", &LNLib::XYZW::IsAlmostEqualTo)
        .def("Distance", &LNLib::XYZW::Distance)
        .def("__getitem__", [](const LNLib::XYZW& xyzw, int index) {
        if (index < 0 || index >= 4) throw py::index_error();
        return xyzw[index];
            })
        .def("__setitem__", [](LNLib::XYZW& xyzw, int index, double value) {
        if (index < 0 || index >= 4) throw py::index_error();
        xyzw[index] = value;
            })
        .def("__add__", [](const LNLib::XYZW& xyzw1, const LNLib::XYZW& xyzw2) {
        return xyzw1 + xyzw2;
            })
        .def("__sub__", [](const LNLib::XYZW& xyzw1, const LNLib::XYZW& xyzw2) {
        return xyzw1 - xyzw2;
            })
        .def("__mul__", [](const LNLib::XYZW& xyzw, double d) {
        return xyzw * d;
            })
        .def("__rmul__", [](const LNLib::XYZW& xyzw, double d) {
        return xyzw * d;
            })
        .def("__truediv__", [](const LNLib::XYZW& xyzw, double d) {
        return xyzw / d;
            })
        .def("__repr__", [](const LNLib::XYZW& xyzw) {
        return "XYZW(" + std::to_string(xyzw.GetWX()) + ", " + std::to_string(xyzw.GetWY()) + ", " + std::to_string(xyzw.GetWZ()) + ", " + std::to_string(xyzw.GetW()) + ")";
            });
}
