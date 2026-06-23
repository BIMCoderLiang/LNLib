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
        .def_static("Horner", [](int degreeU, int degreeV, const std::vector<std::vector<double>>& coefficients, LNLib::UV uv) {
        return LNLib::Polynomials::Horner(degreeU, degreeV, coefficients, uv);
            })
        .def_static("Bernstein", &LNLib::Polynomials::Bernstein)
        .def_static("AllBernstein", &LNLib::Polynomials::AllBernstein)
        .def_static("GetKnotMultiplicity", &LNLib::Polynomials::GetKnotMultiplicity)
        .def_static("GetKnotSpanIndex", &LNLib::Polynomials::GetKnotSpanIndex)
        .def_static("BasisFunctions", [](int spanIndex, int degree, const std::vector<double>& knotVector, double paramT) {
        std::vector<double> basis(LNLib::Constants::NURBSMaxDegree + 1);
        LNLib::Polynomials::BasisFunctions(spanIndex, degree, knotVector, paramT, basis.data());
        return basis;
            })
        .def_static("BasisFunctionsDerivatives", &LNLib::Polynomials::BasisFunctionsDerivatives)
        .def_static("BasisFunctionsFirstOrderDerivative", [](int spanIndex, int degree, const std::vector<double>& knotVector, double paramT) {
        double derivatives[2][LNLib::Constants::NURBSMaxDegree + 1];
        LNLib::Polynomials::BasisFunctionsFirstOrderDerivative(spanIndex, degree, knotVector, paramT, derivatives);
        return std::vector<std::vector<double>>{
            std::vector<double>(derivatives[0], derivatives[0] + LNLib::Constants::NURBSMaxDegree + 1),
                std::vector<double>(derivatives[1], derivatives[1] + LNLib::Constants::NURBSMaxDegree + 1)
        };
            })
        .def_static("OneBasisFunction", &LNLib::Polynomials::OneBasisFunction)
        .def_static("OneBasisFunctionDerivative", &LNLib::Polynomials::OneBasisFunctionDerivative)
        .def_static("AllBasisFunctions", &LNLib::Polynomials::AllBasisFunctions)
        .def_static("BezierToPowerMatrix", &LNLib::Polynomials::BezierToPowerMatrix)
        .def_static("PowerToBezierMatrix", &LNLib::Polynomials::PowerToBezierMatrix);
}
