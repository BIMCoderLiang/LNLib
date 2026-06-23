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
#include <pybind11/operators.h>
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
        .def_property("WX",
            [](const LNLib::XYZW& self) { return self.WX(); },
            [](LNLib::XYZW& self, double val) { self.WX() = val; })
        .def_property("WY",
            [](const LNLib::XYZW& self) { return self.WY(); },
            [](LNLib::XYZW& self, double val) { self.WY() = val; })
        .def_property("WZ",
            [](const LNLib::XYZW& self) { return self.WZ(); },
            [](LNLib::XYZW& self, double val) { self.WZ() = val; })
        .def_property("W",
            [](const LNLib::XYZW& self) { return self.W(); },
            [](LNLib::XYZW& self, double val) { self.W() = val; })
        .def("ToXYZ", &LNLib::XYZW::ToXYZ, py::arg("divideWeight"))
        .def("IsAlmostEqualTo", &LNLib::XYZW::IsAlmostEqualTo)
        .def("Distance", &LNLib::XYZW::Distance)
        .def("__getitem__", [](const LNLib::XYZW& self, int index) { return self[index]; })
        .def("__setitem__", [](LNLib::XYZW& self, int index, double val) { self[index] = val; })
        .def(py::self + py::self)
        .def(py::self - py::self)
        .def(py::self * double())
        .def(double() * py::self)
        .def(py::self / double());
}
