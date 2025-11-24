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
#include "UV.h"
#include "Polynomials.h"

namespace py = pybind11;
void cstrPolynomials(py::module_&m)
{
    py::class_<LNLib::Polynomials>(m, "Polynomials")
        .def_static("Horner", py::overload_cast<int, const std::vector<double>&, double>(&LNLib::Polynomials::Horner))
        .def_static("Bernstein", &LNLib::Polynomials::Bernstein)
        .def_static("AllBernstein", &LNLib::Polynomials::AllBernstein)
        .def_static("Horner", py::overload_cast<int, int, const std::vector<std::vector<double>>&, LNLib::UV&>(&LNLib::Polynomials::Horner))
        .def_static("GetKnotMultiplicity", &LNLib::Polynomials::GetKnotMultiplicity)
        .def_static("GetKnotSpanIndex", &LNLib::Polynomials::GetKnotSpanIndex)
        .def_static("BasisFunctions", [](int spanIndex, int degree, const std::vector<double>& knotVector, double paramT) {
        double basisFunctions[LNLib::Constants::NURBSMaxDegree + 1];
        LNLib::Polynomials::BasisFunctions(spanIndex, degree, knotVector, paramT, basisFunctions);
        std::vector<double> result(basisFunctions, basisFunctions + degree + 1);
        return result;
            })
        .def_static("BasisFunctionsDerivatives", &LNLib::Polynomials::BasisFunctionsDerivatives)
        .def_static("BasisFunctionsFirstOrderDerivative", [](int spanIndex, int degree, const std::vector<double>& knotVector, double paramT) {
        double derivatives[2][LNLib::Constants::NURBSMaxDegree + 1];
        LNLib::Polynomials::BasisFunctionsFirstOrderDerivative(spanIndex, degree, knotVector, paramT, derivatives);

        py::list result;
        for (int i = 0; i < 2; i++) {
            py::list row;
            for (int j = 0; j <= degree; j++) {
                row.append(derivatives[i][j]);
            }
            result.append(row);
        }
        return result;
            })
        .def_static("OneBasisFunction", &LNLib::Polynomials::OneBasisFunction)
        .def_static("OneBasisFunctionDerivative", &LNLib::Polynomials::OneBasisFunctionDerivative)
        .def_static("AllBasisFunctions", &LNLib::Polynomials::AllBasisFunctions)
        .def_static("BezierToPowerMatrix", &LNLib::Polynomials::BezierToPowerMatrix)
        .def_static("PowerToBezierMatrix", &LNLib::Polynomials::PowerToBezierMatrix);
}
