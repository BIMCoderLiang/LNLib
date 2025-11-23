/*
 * Author:
 * 2025/11/21 - Yuqing Liang (BIMCoder Liang)
 * bim.frankliang@foxmail.com
 * 
 *
 * Use of this source code is governed by a LGPL-2.1 license that can be found in
 * the LICENSE file.
 */

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "PyLNLib.h"

namespace py = pybind11;

void construct(py::module_ &);

PYBIND11_MODULE(PyLNLib, m)
{
	construct(m);
}

void construct(py::module_&m)
{
	cstrEnums( m);
	cstrUV(m);
	cstrXYZ(m);
	cstrXYZW(m);
	cstrMatrix4d(m);
	cstrObject(m);
	cstrPolynomials(m);
	cstrValidationUtils(m);
	cstrProjection(m);
	cstrIntersection(m);
	cstrKnotVectorUtils(m);
	cstrBezierCurve(m);
	cstrBezierSurface(m);
	cstrNurbsCurve(m);
	cstrNurbsSurface(m);
}



