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

namespace py = pybind11;

void cstrObject(py::module_& m);
void cstrEnums(py::module_& m);
void cstrUV(py::module_& m);
void cstrXYZ(py::module_& m);
void cstrXYZW(py::module_& m);
void cstrMatrix4d(py::module_& m);
void cstrPolynomials(py::module_& m);
void cstrValidationUtils(py::module_& m);
void cstrProjection(py::module_& m);
void cstrIntersection(py::module_& m);
void cstrKnotVectorUtils(py::module_& m);
void cstrBezierCurve(py::module_& m);
void cstrBezierSurface(py::module_& m);
void cstrNurbsCurve(py::module_& m);
void cstrNurbsSurface(py::module_& m);



