/*
 * Author:
 * 2026/05/24 - Yuqing Liang (BIMCoder Liang)
 * bim.frankliang@foxmail.com
 * 
 *
 * Use of this source code is governed by a LGPL-2.1 license that can be found in
 * the LICENSE file.
 */

#include "NurbsVolume.h"
#include "UVW.h"
#include "XYZ.h"
#include "XYZW.h"
#include "Polynomials.h"
#include "MathUtils.h"
#include "NurbsCurve.h"

LNLib::XYZ LNLib::NurbsVolume::GetPointOnVolume(const LN_NurbsVolume& volume, UVW uvw)
{
	int degreeU = volume.DegreeU;
	int degreeV = volume.DegreeV;
	int degreeW = volume.DegreeW;

	std::vector<double> kvU = volume.KnotVectorU;
	std::vector<double> kvV = volume.KnotVectorV;
	std::vector<double> kvW = volume.KnotVectorW;

	std::vector<std::vector<std::vector<LNLib::XYZW>>> controlPoints = volume.ControlPoints;

	double u = uvw.GetU();
	double v = uvw.GetV();
	double w = uvw.GetW();

	int uspan = LNLib::Polynomials::GetKnotSpanIndex(degreeU, kvU, u);
	double basisU[LNLib::Constants::NURBSMaxDegree + 1];
	LNLib::Polynomials::BasisFunctions(uspan, degreeU, kvU, u, basisU);

	int vspan = LNLib::Polynomials::GetKnotSpanIndex(degreeV, kvV, v);
	double basisV[LNLib::Constants::NURBSMaxDegree + 1];
	LNLib::Polynomials::BasisFunctions(vspan, degreeV, kvV, v, basisV);

	int wspan = LNLib::Polynomials::GetKnotSpanIndex(degreeW, kvW, w);
	double basisW[LNLib::Constants::NURBSMaxDegree + 1];
	LNLib::Polynomials::BasisFunctions(wspan, degreeW, kvW, w, basisW);

	LNLib::XYZW Sw;

	std::vector<LNLib::XYZW> temp1(degreeV + 1);
	std::vector<LNLib::XYZW> temp2(degreeW + 1);

	for (int j = 0; j <= degreeW; ++j)
	{
		std::fill(temp1.begin(), temp1.end(), LNLib::XYZW());

		for (int l = 0; l <= degreeV; ++l)
		{
			for (int k = 0; k <= degreeU; ++k)
			{
				temp1[l] += (basisU[k] * controlPoints[uspan - degreeU + k][vspan - degreeV + l][wspan - degreeW + j]);
			}
		}
		for (int l = 0; l <= degreeV; ++l)
		{
			temp2[j] += basisV[l] * temp1[l];
		}
	}
	for (int j = 0; j <= degreeW; ++j)
	{
		Sw += basisW[j] * temp2[j];
	}

	return Sw.ToXYZ(true);
}

std::vector<std::vector<std::vector<LNLib::XYZ>>> LNLib::NurbsVolume::ComputeVolumeRationalDerivatives(const LN_NurbsVolume& volume, int derivative, UVW uvw)
{
	int degreeU = volume.DegreeU;
	int degreeV = volume.DegreeV;
	int degreeW = volume.DegreeW;

	std::vector<double> kvU = volume.KnotVectorU;
	std::vector<double> kvV = volume.KnotVectorV;
	std::vector<double> kvW = volume.KnotVectorW;

	std::vector<std::vector<std::vector<LNLib::XYZW>>> controlPoints = volume.ControlPoints;

	std::vector<std::vector<std::vector<LNLib::XYZW>>> ders(
		derivative + 1,
		std::vector<std::vector<LNLib::XYZW>>(
			derivative + 1,
			std::vector<LNLib::XYZW>(derivative + 1)
		)
	);

	std::vector<LNLib::XYZW> temp1(degreeV + 1);
	std::vector<LNLib::XYZW> temp2(degreeW + 1);

	int du = std::min(derivative, degreeU);
	int dv = std::min(derivative, degreeV);
	int dw = std::min(derivative, degreeW);

	double u = uvw.GetU();
	double v = uvw.GetV();
	double w = uvw.GetW();

	int uspan = LNLib::Polynomials::GetKnotSpanIndex(degreeU, kvU, u);
	std::vector<std::vector<double>> Nu = LNLib::Polynomials::BasisFunctionsDerivatives(uspan, degreeU, derivative, kvU, u);

	int vspan = LNLib::Polynomials::GetKnotSpanIndex(degreeV, kvV, v);
	std::vector<std::vector<double>> Nv = LNLib::Polynomials::BasisFunctionsDerivatives(vspan, degreeV, derivative, kvV, v);

	int wspan = LNLib::Polynomials::GetKnotSpanIndex(degreeW, kvW, w);
	std::vector<std::vector<double>> Nw = LNLib::Polynomials::BasisFunctionsDerivatives(wspan, degreeW, derivative, kvW, w);

	for (int a = 0; a <= du; ++a)
	{
		int ddb = std::min(derivative - a, dv);
		for (int b = 0; b <= ddb; ++b)
		{
			std::fill(temp2.begin(), temp2.end(), LNLib::XYZW());
			for (int z = 0; z <= degreeW; ++z)
			{
				std::fill(temp1.begin(), temp1.end(), LNLib::XYZW());
				for (int y = 0; y <= degreeV; ++y)
				{
					for (int x = 0; x <= degreeU; ++x)
					{
						temp1[y] += Nu[a][x] * controlPoints[uspan - degreeU + x][vspan - degreeV + y][wspan - degreeW + z];
					}
				}
				for (int y = 0; y <= degreeV; ++y)
				{
					temp2[z] += Nv[b][y] * temp2[y];
				}
			}
			int ddc = std::min(derivative - a - b, dw);
			for (int c = 0; c <= ddc; ++c)
			{
				for (int z = 0; z <= degreeW; ++z)
				{
					ders[a][b][c] += Nw[c][z] * temp2[z];
				}
			}
		}
	}

	std::vector<std::vector<std::vector<LNLib::XYZ>>> Aders(
		derivative + 1,
		std::vector<std::vector<LNLib::XYZ>>(
			derivative + 1,
			std::vector<LNLib::XYZ>(derivative + 1)
		)
	);

	std::vector<std::vector<std::vector<double>>> wders(
		derivative + 1,
		std::vector<std::vector<double>>(
			derivative + 1,
			std::vector<double>(derivative + 1)
		)
	);

	std::vector<std::vector<std::vector<LNLib::XYZ>>> derivatives(
		derivative + 1,
		std::vector<std::vector<LNLib::XYZ>>(
			derivative + 1,
			std::vector<LNLib::XYZ>(derivative + 1)
		)
	);

	for (int i = 0; i < ders.size(); i++)
	{
		for (int j = 0; j < ders[0].size(); j++)
		{
			for (int k = 0; k < ders[0][0].size(); k++)
			{
				Aders[i][j][k] = ders[i][j][k].ToXYZ(false);
				wders[i][j][k] = ders[i][j][k].GetW();
			}
		}
	}

	for (int a = 0; a <= derivative; a++)
	{
		for (int b = 0; b <= derivative - a; b++)
		{
			for (int c = 0; c <= derivative - a - b; c++)
			{
				LNLib::XYZ g = Aders[a][b][c];
				for (int s = 1; s <= c; s++)
				{
					g -= LNLib::MathUtils::Binomial(c, s) * wders[0][0][s] * derivatives[a][b][c - s];
				}
				for (int j = 1; j <= b; j++)
				{
					for (int s = 0; s <= c; ++s)
					{
						g -= LNLib::MathUtils::Binomial(b, j) * LNLib::MathUtils::Binomial(c, s) *
							wders[0][j][s] * derivatives[a][b - j][c - s];
					}
				}
				for (int i = 1; i <= a; i++)
				{
					for (int j = 0; j <= b; j++)
					{
						for (int s = 0; s <= c; s++)
						{
							g -= LNLib::MathUtils::Binomial(a, i) * LNLib::MathUtils::Binomial(b, j) * LNLib::MathUtils::Binomial(c, s) *
								wders[i][j][s] * derivatives[a - i][b - j][c - s];
						}
					}
				}
				derivatives[a][b][c] = g / wders[0][0][0];
			}
		}
	}
	return derivatives;
}

void LNLib::NurbsVolume::InsertKnot(const LN_NurbsVolume& volume, double insertKnot, int times, VolumeDirection direction, LN_NurbsVolume& result)
{
	int degreeU = volume.DegreeU;
	int degreeV = volume.DegreeV;
	int degreeW = volume.DegreeW;

	std::vector<double> kvU = volume.KnotVectorU;
	std::vector<double> kvV = volume.KnotVectorV;
	std::vector<double> kvW = volume.KnotVectorW;

	std::vector<std::vector<std::vector<LNLib::XYZW>>> controlPoints = volume.ControlPoints;

	int n = controlPoints.size() - 1;
	int m = controlPoints[0].size() - 1;
	int l = controlPoints[0][0].size() - 1;

	result.DegreeU = degreeU;
	result.DegreeV = degreeV;
	result.DegreeW = degreeW;

	std::vector<std::vector<std::vector<LNLib::XYZW>>> newControlPoints;
	if (direction == VolumeDirection::U)
	{
		newControlPoints.resize(n + times + 1, std::vector<std::vector<LNLib::XYZW>>(m + 1, std::vector<LNLib::XYZW>(l + 1)));
	}
	else if (direction == VolumeDirection::V)
	{
		newControlPoints.resize(n + 1, std::vector<std::vector<LNLib::XYZW>>(m + times + 1, std::vector<LNLib::XYZW>(l + 1)));
	}
	else if (direction == VolumeDirection::W)
	{
		newControlPoints.resize(n + 1, std::vector<std::vector<LNLib::XYZW>>(m + 1, std::vector<LNLib::XYZW>(l + times + 1)));
	}

	if (direction == VolumeDirection::U)
	{
		result.KnotVectorV = kvV;
		result.KnotVectorW = kvW;

		std::vector<double> newKvU;

		for (int sep = 0; sep <= l; sep++)
		{
			for (int col = 0; col <= m; col++)
			{
				std::vector<LNLib::XYZW> uControlPoints(n + 1);
				for (int i = 0; i <= n; ++i) {
					uControlPoints[i] = controlPoints[i][col][sep];
				}

				LNLib::LN_NurbsCurve curve;
				curve.Degree = degreeU;
				curve.KnotVector = kvU;
				curve.ControlPoints = uControlPoints;

				LNLib::LN_NurbsCurve resultCurve;
				LNLib::NurbsCurve::InsertKnot(curve, insertKnot, times, resultCurve);

				newKvU = resultCurve.KnotVector;
				std::vector<LNLib::XYZW> iControlPoints = resultCurve.ControlPoints;

				for (int i = 0; i < iControlPoints.size(); i++)
				{
					newControlPoints[i][col][sep] = iControlPoints[i];
				}
			}
		}

		result.KnotVectorU = newKvU;
		result.ControlPoints = newControlPoints;
	}
	else if (direction == VolumeDirection::V)
	{
		result.KnotVectorU = kvU;
		result.KnotVectorW = kvW;

		std::vector<double> newKvV;

		for (int sep = 0; sep <= l; sep++)
		{
			for (int row = 0; row <= n; row++)
			{
				std::vector<LNLib::XYZW> vControlPoints(m + 1);
				for (int j = 0; j <= m; ++j) {
					vControlPoints[j] = controlPoints[row][j][sep];
				}

				LNLib::LN_NurbsCurve curve;
				curve.Degree = degreeV;
				curve.KnotVector = kvV;
				curve.ControlPoints = vControlPoints;

				LNLib::LN_NurbsCurve resultCurve;
				LNLib::NurbsCurve::InsertKnot(curve, insertKnot, times, resultCurve);

				newKvV = resultCurve.KnotVector;
				std::vector<LNLib::XYZW> iControlPoints = resultCurve.ControlPoints;

				for (int j = 0; j < iControlPoints.size(); j++)
				{
					newControlPoints[row][j][sep] = iControlPoints[j];
				}
			}
		}

		result.KnotVectorV = newKvV;
		result.ControlPoints = newControlPoints;
	}
	else if (direction == VolumeDirection::W)
	{
		result.KnotVectorU = kvU;
		result.KnotVectorV = kvV;

		std::vector<double> newKvW;

		for (int col = 0; col <= m; col++)
		{
			for (int row = 0; row <= n; row++)
			{
				std::vector<LNLib::XYZW> wControlPoints(l + 1);
				for (int k = 0; k <= l; ++k) {
					wControlPoints[k] = controlPoints[row][col][k];
				}

				LNLib::LN_NurbsCurve curve;
				curve.Degree = degreeW;
				curve.KnotVector = kvW;
				curve.ControlPoints = wControlPoints;

				LNLib::LN_NurbsCurve resultCurve;
				LNLib::NurbsCurve::InsertKnot(curve, insertKnot, times, resultCurve);

				newKvW = resultCurve.KnotVector;
				std::vector<LNLib::XYZW> iControlPoints = resultCurve.ControlPoints;

				for (int k = 0; k < iControlPoints.size(); k++)
				{
					newControlPoints[row][col][k] = iControlPoints[k];
				}
			}
		}

		result.KnotVectorW = newKvW;
		result.ControlPoints = newControlPoints;
	}
}

bool LNLib::NurbsVolume::SplitAt(const LN_NurbsVolume& volume, double parameter, VolumeDirection direction, LN_NurbsVolume& first, LN_NurbsVolume& second)
{
	int degreeU = volume.DegreeU;
	int degreeV = volume.DegreeV;
	int degreeW = volume.DegreeW;

	std::vector<double> kvU = volume.KnotVectorU;
	std::vector<double> kvV = volume.KnotVectorV;
	std::vector<double> kvW = volume.KnotVectorW;

	std::vector<std::vector<std::vector<LNLib::XYZW>>> controlPoints = volume.ControlPoints;

	LN_NurbsVolume insertVolume;
	int k, s, rr;

	if (direction == VolumeDirection::U)
	{
		if (LNLib::MathUtils::IsAlmostEqualTo(parameter, kvU[0]) ||
			LNLib::MathUtils::IsAlmostEqualTo(parameter, kvU[kvU.size() - 1]))
		{
			first = volume;
			second = volume;
			return true;
		}

		k = LNLib::Polynomials::GetKnotSpanIndex(degreeU, kvU, parameter);
		s = LNLib::Polynomials::GetKnotMultiplicity(kvU, parameter);
		rr = degreeU - s;

		if (rr > 0)
		{
			InsertKnot(volume, parameter, rr, direction, insertVolume);

			kvU = insertVolume.KnotVectorU;
			kvV = insertVolume.KnotVectorV;
			kvW = insertVolume.KnotVectorW;
			controlPoints = insertVolume.ControlPoints;

			k = LNLib::Polynomials::GetKnotSpanIndex(degreeU, kvU, parameter);
			s = LNLib::Polynomials::GetKnotMultiplicity(kvU, parameter);
		}

		std::vector<double> kU_L = kvU;
		kU_L.resize(k + 2);
		kU_L.back() = parameter;

		std::vector<double> kU_R;
		kU_R.push_back(parameter);
		kU_R.insert(kU_R.end(), kvU.begin() + (k - s + 1), kvU.end());

		std::vector<std::vector<std::vector<LNLib::XYZW>>> cp_L(controlPoints.begin(), controlPoints.begin() + (k - s + 1));
		std::vector<std::vector<std::vector<LNLib::XYZW>>> cp_R(controlPoints.begin() + (k - s), controlPoints.end());

		first.DegreeU = degreeU; first.DegreeV = degreeV; first.DegreeW = degreeW;
		first.KnotVectorU = kU_L; first.KnotVectorV = kvV; first.KnotVectorW = kvW;
		first.ControlPoints = cp_L;

		second.DegreeU = degreeU; second.DegreeV = degreeV; second.DegreeW = degreeW;
		second.KnotVectorU = kU_R; second.KnotVectorV = kvV; second.KnotVectorW = kvW;
		second.ControlPoints = cp_R;
	}
	else if (direction == VolumeDirection::V)
	{
		if (LNLib::MathUtils::IsAlmostEqualTo(parameter, kvV[0]) ||
			LNLib::MathUtils::IsAlmostEqualTo(parameter, kvV[kvV.size() - 1]))
		{
			first = volume;
			second = volume;
			return true;
		}

		k = LNLib::Polynomials::GetKnotSpanIndex(degreeV, kvV, parameter);
		s = LNLib::Polynomials::GetKnotMultiplicity(kvV, parameter);
		rr = degreeV - s;

		if (rr > 0)
		{
			InsertKnot(volume, parameter, rr, direction, insertVolume);

			kvU = insertVolume.KnotVectorU;
			kvV = insertVolume.KnotVectorV;
			kvW = insertVolume.KnotVectorW;
			controlPoints = insertVolume.ControlPoints;

			k = LNLib::Polynomials::GetKnotSpanIndex(degreeV, kvV, parameter);
			s = LNLib::Polynomials::GetKnotMultiplicity(kvV, parameter);
		}

		std::vector<double> kV_L = kvV;
		kV_L.resize(k + 2);
		kV_L.back() = parameter;

		std::vector<double> kV_R;
		kV_R.push_back(parameter);
		kV_R.insert(kV_R.end(), kvV.begin() + (k - s + 1), kvV.end());

		int splitIdx = k - s;
		std::vector<std::vector<std::vector<LNLib::XYZW>>> cp_L;
		std::vector<std::vector<std::vector<LNLib::XYZW>>> cp_R;

		cp_L.resize(controlPoints.size());
		cp_R.resize(controlPoints.size());

		for (size_t i = 0; i < controlPoints.size(); ++i)
		{
			cp_L[i] = std::vector<std::vector<LNLib::XYZW>>(controlPoints[i].begin(), controlPoints[i].begin() + splitIdx + 1);
			cp_R[i] = std::vector<std::vector<LNLib::XYZW>>(controlPoints[i].begin() + splitIdx, controlPoints[i].end());
		}

		first.DegreeU = degreeU; first.DegreeV = degreeV; first.DegreeW = degreeW;
		first.KnotVectorU = kvU; first.KnotVectorV = kV_L; first.KnotVectorW = kvW;
		first.ControlPoints = cp_L;

		second.DegreeU = degreeU; second.DegreeV = degreeV; second.DegreeW = degreeW;
		second.KnotVectorU = kvU; second.KnotVectorV = kV_R; second.KnotVectorW = kvW;
		second.ControlPoints = cp_R;
	}
	else if (direction == VolumeDirection::W)
	{
		if (LNLib::MathUtils::IsAlmostEqualTo(parameter, kvW[0]) ||
			LNLib::MathUtils::IsAlmostEqualTo(parameter, kvW[kvW.size() - 1]))
		{
			first = volume;
			second = volume;
			return true;
		}

		k = LNLib::Polynomials::GetKnotSpanIndex(degreeW, kvW, parameter);
		s = LNLib::Polynomials::GetKnotMultiplicity(kvW, parameter);
		rr = degreeW - s;

		if (rr > 0)
		{
			InsertKnot(volume, parameter, rr, direction, insertVolume);

			kvU = insertVolume.KnotVectorU;
			kvV = insertVolume.KnotVectorV;
			kvW = insertVolume.KnotVectorW;
			controlPoints = insertVolume.ControlPoints;

			k = LNLib::Polynomials::GetKnotSpanIndex(degreeW, kvW, parameter);
			s = LNLib::Polynomials::GetKnotMultiplicity(kvW, parameter);
		}

		std::vector<double> kW_L = kvW;
		kW_L.resize(k + 2);
		kW_L.back() = parameter;

		std::vector<double> kW_R;
		kW_R.push_back(parameter);
		kW_R.insert(kW_R.end(), kvW.begin() + (k - s + 1), kvW.end());

		int splitIdx = k - s;
		std::vector<std::vector<std::vector<LNLib::XYZW>>> cp_L;
		std::vector<std::vector<std::vector<LNLib::XYZW>>> cp_R;

		cp_L.resize(controlPoints.size());
		cp_R.resize(controlPoints.size());

		for (size_t i = 0; i < controlPoints.size(); ++i)
		{
			cp_L[i].resize(controlPoints[i].size());
			cp_R[i].resize(controlPoints[i].size());

			for (size_t j = 0; j < controlPoints[i].size(); ++j)
			{
				cp_L[i][j] = std::vector<LNLib::XYZW>(controlPoints[i][j].begin(), controlPoints[i][j].begin() + splitIdx + 1);
				cp_R[i][j] = std::vector<LNLib::XYZW>(controlPoints[i][j].begin() + splitIdx, controlPoints[i][j].end());
			}
		}

		first.DegreeU = degreeU; first.DegreeV = degreeV; first.DegreeW = degreeW;
		first.KnotVectorU = kvU; first.KnotVectorV = kvV; first.KnotVectorW = kW_L;
		first.ControlPoints = cp_L;

		second.DegreeU = degreeU; second.DegreeV = degreeV; second.DegreeW = degreeW;
		second.KnotVectorU = kvU; second.KnotVectorV = kvV; second.KnotVectorW = kW_R;
		second.ControlPoints = cp_R;
	}
	else
	{
		return false;
	}

	return true;
}

bool LNLib::NurbsVolume::Extract(const LN_NurbsVolume& volume, double parameter, VolumeDirection direction, LNLib::LN_NurbsSurface surface)
{
	int degreeU = volume.DegreeU;
	int degreeV = volume.DegreeV;
	int degreeW = volume.DegreeW;

	std::vector<double> kvU = volume.KnotVectorU;
	std::vector<double> kvV = volume.KnotVectorV;
	std::vector<double> kvW = volume.KnotVectorW;

	std::vector<std::vector<std::vector<LNLib::XYZW>>> controlPoints = volume.ControlPoints;
	std::vector<std::vector<LNLib::XYZW>> Pwc;

	LN_NurbsVolume insertVolume;
	int k, s, rr;

	int newDegreeU, newDegreeV;
	std::vector<double> newKvU, newKvV;

	if (direction == VolumeDirection::U)
	{
		if (LNLib::MathUtils::IsAlmostEqualTo(parameter, kvU[0]))
		{
			Pwc = controlPoints[0];
		}
		else if (LNLib::MathUtils::IsAlmostEqualTo(parameter, kvU[kvU.size() - 1]))
		{
			Pwc = controlPoints[controlPoints.size() - 1];
		}
		else
		{
			k = LNLib::Polynomials::GetKnotSpanIndex(degreeU, kvU, parameter);
			s = LNLib::Polynomials::GetKnotMultiplicity(kvU, parameter);
			rr = degreeU - s;

			if (rr > 0)
			{
				InsertKnot(volume, parameter, rr, direction, insertVolume);

				kvU = insertVolume.KnotVectorU;
				kvV = insertVolume.KnotVectorV;
				kvW = insertVolume.KnotVectorW;
				controlPoints = insertVolume.ControlPoints;

				k = LNLib::Polynomials::GetKnotSpanIndex(degreeU, kvU, parameter);
				s = LNLib::Polynomials::GetKnotMultiplicity(kvU, parameter);
			}
			Pwc = controlPoints[k - s];

			newDegreeU = degreeV;
			newDegreeV = degreeW;
			newKvU = kvV;
			newKvV = kvW;
		}
	}
	else if (direction == VolumeDirection::V)
	{
		if (LNLib::MathUtils::IsAlmostEqualTo(parameter, kvV[0]))
		{
			Pwc.resize(controlPoints.size());
			for (size_t i = 0; i < controlPoints.size(); ++i) {
				Pwc[i] = controlPoints[i][0];
			}
		}
		else if (LNLib::MathUtils::IsAlmostEqualTo(parameter, kvV[kvV.size() - 1]))
		{
			Pwc.resize(controlPoints.size());
			for (size_t i = 0; i < controlPoints.size(); ++i) {
				Pwc[i] = controlPoints[i][controlPoints[i].size() - 1];
			}
		}
		else
		{
			k = LNLib::Polynomials::GetKnotSpanIndex(degreeV, kvV, parameter);
			s = LNLib::Polynomials::GetKnotMultiplicity(kvV, parameter);
			rr = degreeV - s;

			if (rr > 0)
			{
				InsertKnot(volume, parameter, rr, direction, insertVolume);

				kvU = insertVolume.KnotVectorU;
				kvV = insertVolume.KnotVectorV;
				kvW = insertVolume.KnotVectorW;
				controlPoints = insertVolume.ControlPoints;

				k = LNLib::Polynomials::GetKnotSpanIndex(degreeV, kvV, parameter);
				s = LNLib::Polynomials::GetKnotMultiplicity(kvV, parameter);
			}
			Pwc.resize(controlPoints.size());
			for (size_t i = 0; i < controlPoints.size(); ++i) {
				Pwc[i] = controlPoints[i][k - s];
			}

			std::vector<std::vector<LNLib::XYZW>> tempPwc(Pwc.size());
			for (size_t i = 0; i < Pwc.size(); ++i) {
				tempPwc[i].resize(Pwc[i].size());
				for (size_t j = 0; j < Pwc[i].size(); ++j) {
					tempPwc[i][j] = Pwc[i][j];
				}
			}
			Pwc = tempPwc;

			newDegreeU = degreeW;
			newDegreeV = degreeU;
			newKvU = kvW;
			newKvV = kvU;
		}
	}
	else if (direction == VolumeDirection::W)
	{
		if (LNLib::MathUtils::IsAlmostEqualTo(parameter, kvW[0]))
		{
			Pwc.resize(controlPoints.size());
			for (size_t i = 0; i < controlPoints.size(); ++i) {
				Pwc[i].resize(controlPoints[i].size());
				for (size_t j = 0; j < controlPoints[i].size(); ++j) {
					Pwc[i][j] = controlPoints[i][j][0];
				}
			}
		}
		else if (LNLib::MathUtils::IsAlmostEqualTo(parameter, kvW[kvW.size() - 1]))
		{
			Pwc.resize(controlPoints.size());
			for (size_t i = 0; i < controlPoints.size(); ++i) {
				Pwc[i].resize(controlPoints[i].size());
				for (size_t j = 0; j < controlPoints[i].size(); ++j) {
					Pwc[i][j] = controlPoints[i][j][controlPoints[i][j].size() - 1];
				}
			}
		}
		else
		{
			k = LNLib::Polynomials::GetKnotSpanIndex(degreeW, kvW, parameter);
			s = LNLib::Polynomials::GetKnotMultiplicity(kvW, parameter);
			rr = degreeW - s;

			if (rr > 0)
			{
				InsertKnot(volume, parameter, rr, direction, insertVolume);

				kvU = insertVolume.KnotVectorU;
				kvV = insertVolume.KnotVectorV;
				kvW = insertVolume.KnotVectorW;
				controlPoints = insertVolume.ControlPoints;

				k = LNLib::Polynomials::GetKnotSpanIndex(degreeW, kvW, parameter);
				s = LNLib::Polynomials::GetKnotMultiplicity(kvW, parameter);
			}
			Pwc.resize(controlPoints.size());
			for (size_t i = 0; i < controlPoints.size(); ++i) {
				Pwc[i].resize(controlPoints[i].size());
				for (size_t j = 0; j < controlPoints[i].size(); ++j) {
					Pwc[i][j] = controlPoints[i][j][k - s];
				}
			}

			newDegreeU = degreeU;
			newDegreeV = degreeV;
			newKvU = kvU;
			newKvV = kvV;
		}
	}

	surface.DegreeU = newDegreeU;
	surface.DegreeV = newDegreeV;
	surface.KnotVectorU = newKvU;
	surface.KnotVectorV = newKvV;
	surface.ControlPoints = Pwc;

	return true;
}

void LNLib::NurbsVolume::RefineKnotVector(const LN_NurbsVolume& volume, const std::vector<double>& insertKnotElements, VolumeDirection direction, LN_NurbsVolume& result)
{
	int degreeU = volume.DegreeU;
	int degreeV = volume.DegreeV;
	int degreeW = volume.DegreeW;

	std::vector<double> kvU = volume.KnotVectorU;
	std::vector<double> kvV = volume.KnotVectorV;
	std::vector<double> kvW = volume.KnotVectorW;

	std::vector<std::vector<std::vector<LNLib::XYZW>>> controlPoints = volume.ControlPoints;

	int n = controlPoints.size() - 1;
	int m = controlPoints[0].size() - 1;
	int l = controlPoints[0][0].size() - 1;

	result.DegreeU = degreeU;
	result.DegreeV = degreeV;
	result.DegreeW = degreeW;

	std::vector<std::vector<std::vector<LNLib::XYZW>>> newControlPoints;

	int times = insertKnotElements.size() - 1;
	if (direction == VolumeDirection::U)
	{
		newControlPoints.resize(n + times + 1, std::vector<std::vector<LNLib::XYZW>>(m + 1, std::vector<LNLib::XYZW>(l + 1)));
	}
	else if (direction == VolumeDirection::V)
	{
		newControlPoints.resize(n + 1, std::vector<std::vector<LNLib::XYZW>>(m + times + 1, std::vector<LNLib::XYZW>(l + 1)));
	}
	else if (direction == VolumeDirection::W)
	{
		newControlPoints.resize(n + 1, std::vector<std::vector<LNLib::XYZW>>(m + 1, std::vector<LNLib::XYZW>(l + times + 1)));
	}

	if (direction == VolumeDirection::U)
	{
		result.KnotVectorV = kvV;
		result.KnotVectorW = kvW;

		std::vector<double> newKvU;

		for (int sep = 0; sep <= l; sep++)
		{
			for (int col = 0; col <= m; col++)
			{
				std::vector<LNLib::XYZW> uControlPoints(n + 1);
				for (int i = 0; i <= n; ++i) {
					uControlPoints[i] = controlPoints[i][col][sep];
				}

				LNLib::LN_NurbsCurve curve;
				curve.Degree = degreeU;
				curve.KnotVector = kvU;
				curve.ControlPoints = uControlPoints;

				LNLib::LN_NurbsCurve resultCurve;
				LNLib::NurbsCurve::RefineKnotVector(curve, insertKnotElements, resultCurve);

				newKvU = resultCurve.KnotVector;
				std::vector<LNLib::XYZW> iControlPoints = resultCurve.ControlPoints;

				for (int i = 0; i < iControlPoints.size(); i++)
				{
					newControlPoints[i][col][sep] = iControlPoints[i];
				}
			}
		}

		result.KnotVectorU = newKvU;
		result.ControlPoints = newControlPoints;
	}
	else if (direction == VolumeDirection::V)
	{
		result.KnotVectorU = kvU;
		result.KnotVectorW = kvW;

		std::vector<double> newKvV;

		for (int sep = 0; sep <= l; sep++)
		{
			for (int row = 0; row <= n; row++)
			{
				std::vector<LNLib::XYZW> vControlPoints(m + 1);
				for (int j = 0; j <= m; ++j) {
					vControlPoints[j] = controlPoints[row][j][sep];
				}

				LNLib::LN_NurbsCurve curve;
				curve.Degree = degreeV;
				curve.KnotVector = kvV;
				curve.ControlPoints = vControlPoints;

				LNLib::LN_NurbsCurve resultCurve;
				LNLib::NurbsCurve::RefineKnotVector(curve, insertKnotElements, resultCurve);

				newKvV = resultCurve.KnotVector;
				std::vector<LNLib::XYZW> iControlPoints = resultCurve.ControlPoints;

				for (int j = 0; j < iControlPoints.size(); j++)
				{
					newControlPoints[row][j][sep] = iControlPoints[j];
				}
			}
		}

		result.KnotVectorV = newKvV;
		result.ControlPoints = newControlPoints;
	}
	else if (direction == VolumeDirection::W)
	{
		result.KnotVectorU = kvU;
		result.KnotVectorV = kvV;

		std::vector<double> newKvW;

		for (int col = 0; col <= m; col++)
		{
			for (int row = 0; row <= n; row++)
			{
				std::vector<LNLib::XYZW> wControlPoints(l + 1);
				for (int k = 0; k <= l; ++k) {
					wControlPoints[k] = controlPoints[row][col][k];
				}

				LNLib::LN_NurbsCurve curve;
				curve.Degree = degreeW;
				curve.KnotVector = kvW;
				curve.ControlPoints = wControlPoints;

				LNLib::LN_NurbsCurve resultCurve;
				LNLib::NurbsCurve::RefineKnotVector(curve, insertKnotElements, resultCurve);

				newKvW = resultCurve.KnotVector;
				std::vector<LNLib::XYZW> iControlPoints = resultCurve.ControlPoints;

				for (int k = 0; k < iControlPoints.size(); k++)
				{
					newControlPoints[row][col][k] = iControlPoints[k];
				}
			}
		}

		result.KnotVectorW = newKvW;
		result.ControlPoints = newControlPoints;
	}
}

void LNLib::NurbsVolume::ElevateDegree(const LN_NurbsVolume& volume, int times, VolumeDirection direction, LN_NurbsVolume& result)
{
	int degreeU = volume.DegreeU;
	int degreeV = volume.DegreeV;
	int degreeW = volume.DegreeW;

	std::vector<double> kvU = volume.KnotVectorU;
	std::vector<double> kvV = volume.KnotVectorV;
	std::vector<double> kvW = volume.KnotVectorW;

	std::vector<std::vector<std::vector<LNLib::XYZW>>> controlPoints = volume.ControlPoints;

	int n = controlPoints.size() - 1;
	int m = controlPoints[0].size() - 1;
	int l = controlPoints[0][0].size() - 1;

	std::vector<std::vector<std::vector<LNLib::XYZW>>> newControlPoints;

	if (direction == VolumeDirection::U)
	{
		std::vector<double> newkvU;
		std::vector<double> newkvV = kvV;
		std::vector<double> newkvW = kvW;

		int startIdx = degreeU + 1;
		int endIdx = static_cast<int>(kvU.size()) - (degreeU + 1);
		int uniqueCount = endIdx - startIdx;
		int nh = n + (uniqueCount + 1) * times;

		newControlPoints.resize(nh + 1, std::vector<std::vector<LNLib::XYZW>>(m + 1, std::vector<LNLib::XYZW>(l + 1)));

		for (int sep = 0; sep <= l; sep++)
		{
			for (int col = 0; col <= m; col++)
			{
				std::vector<LNLib::XYZW> uControlPoints(n + 1);
				for (int i = 0; i <= n; ++i) {
					uControlPoints[i] = controlPoints[i][col][sep];
				}

				LNLib::LN_NurbsCurve curve;
				curve.Degree = degreeU;
				curve.KnotVector = kvU;
				curve.ControlPoints = uControlPoints;

				LNLib::LN_NurbsCurve resultCurve;
				LNLib::NurbsCurve::ElevateDegree(curve, times, resultCurve);

				newkvU = resultCurve.KnotVector;
				std::vector<LNLib::XYZW> iControlPoints = resultCurve.ControlPoints;

				for (int i = 0; i < iControlPoints.size(); i++)
				{
					newControlPoints[i][col][sep] = iControlPoints[i];
				}
			}
		}

		result.DegreeU = degreeU + times;
		result.DegreeV = degreeV;
		result.DegreeW = degreeW;

		result.KnotVectorU = newkvU;
		result.KnotVectorV = newkvV;
		result.KnotVectorW = newkvW;

		result.ControlPoints = newControlPoints;
	}
	else if (direction == VolumeDirection::V)
	{
		std::vector<double> newkvU = kvU;
		std::vector<double> newkvV;
		std::vector<double> newkvW = kvW;

		int startIdx = degreeV + 1;
		int endIdx = static_cast<int>(kvV.size()) - (degreeV + 1);
		int uniqueCount = endIdx - startIdx;
		int mh = m + (uniqueCount + 1) * times;

		newControlPoints.resize(n + 1, std::vector<std::vector<LNLib::XYZW>>(mh + 1, std::vector<LNLib::XYZW>(l + 1)));

		for (int sep = 0; sep <= l; sep++)
		{
			for (int row = 0; row <= n; row++)
			{
				std::vector<LNLib::XYZW> vControlPoints(m + 1);
				for (int j = 0; j <= m; ++j) {
					vControlPoints[j] = controlPoints[row][j][sep];
				}

				LNLib::LN_NurbsCurve curve;
				curve.Degree = degreeV;
				curve.KnotVector = kvV;
				curve.ControlPoints = vControlPoints;

				LNLib::LN_NurbsCurve resultCurve;
				LNLib::NurbsCurve::ElevateDegree(curve, times, resultCurve);

				newkvV = resultCurve.KnotVector;
				std::vector<LNLib::XYZW> iControlPoints = resultCurve.ControlPoints;

				for (int j = 0; j < iControlPoints.size(); j++)
				{
					newControlPoints[row][j][sep] = iControlPoints[j];
				}
			}
		}

		result.DegreeU = degreeU;
		result.DegreeV = degreeV + times;
		result.DegreeW = degreeW;

		result.KnotVectorU = newkvU;
		result.KnotVectorV = newkvV;
		result.KnotVectorW = newkvW;

		result.ControlPoints = newControlPoints;
	}
	else if (direction == VolumeDirection::W)
	{
		std::vector<double> newkvU = kvU;
		std::vector<double> newkvV = kvV;
		std::vector<double> newkvW;

		int startIdx = degreeW + 1;
		int endIdx = static_cast<int>(kvW.size()) - (degreeW + 1);
		int uniqueCount = endIdx - startIdx;
		int lh = l + (uniqueCount + 1) * times;

		newControlPoints.resize(n + 1, std::vector<std::vector<LNLib::XYZW>>(m + 1, std::vector<LNLib::XYZW>(lh + 1)));

		for (int col = 0; col <= m; col++)
		{
			for (int row = 0; row <= n; row++)
			{
				std::vector<LNLib::XYZW> wControlPoints(l + 1);
				for (int k = 0; k <= l; ++k) {
					wControlPoints[k] = controlPoints[row][col][k];
				}

				LNLib::LN_NurbsCurve curve;
				curve.Degree = degreeW;
				curve.KnotVector = kvW;
				curve.ControlPoints = wControlPoints;

				LNLib::LN_NurbsCurve resultCurve;
				LNLib::NurbsCurve::ElevateDegree(curve, times, resultCurve);

				newkvW = resultCurve.KnotVector;
				std::vector<LNLib::XYZW> iControlPoints = resultCurve.ControlPoints;

				for (int k = 0; k < iControlPoints.size(); k++)
				{
					newControlPoints[row][col][k] = iControlPoints[k];
				}
			}
		}

		result.DegreeU = degreeU;
		result.DegreeV = degreeV;
		result.DegreeW = degreeW + times;

		result.KnotVectorU = newkvU;
		result.KnotVectorV = newkvV;
		result.KnotVectorW = newkvW;

		result.ControlPoints = newControlPoints;
	}
}

bool LNLib::NurbsVolume::IsClosed(const LN_NurbsVolume& volume, VolumeDirection direction)
{
	int numU = volume.ControlPoints.size();
	if (numU == 0) return false;
	int numV = volume.ControlPoints[0].size();
	if (numV == 0) return false;
	int numW = volume.ControlPoints[0][0].size();

	switch (direction)
	{
	case VolumeDirection::U:
	{
		for (int v = 0; v < numV; ++v)
		{
			for (int w = 0; w < numW; ++w)
			{
				std::vector<XYZW> uCurveCPs;
				uCurveCPs.reserve(numU);
				for (int u = 0; u < numU; ++u)
				{
					uCurveCPs.push_back(volume.ControlPoints[u][v][w]);
				}

				LN_NurbsCurve curve;
				curve.Degree = volume.DegreeU;
				curve.KnotVector = volume.KnotVectorU;
				curve.ControlPoints = uCurveCPs;

				if (!NurbsCurve::IsClosed(curve))
				{
					return false;
				}
			}
		}
		return true;
	}

	case VolumeDirection::V:
	{
		for (int u = 0; u < numU; ++u)
		{
			for (int w = 0; w < numW; ++w)
			{
				std::vector<XYZW> vCurveCPs;
				vCurveCPs.reserve(numV);
				for (int v = 0; v < numV; ++v)
				{
					vCurveCPs.push_back(volume.ControlPoints[u][v][w]);
				}

				LN_NurbsCurve curve;
				curve.Degree = volume.DegreeV;
				curve.KnotVector = volume.KnotVectorV;
				curve.ControlPoints = vCurveCPs;

				if (!NurbsCurve::IsClosed(curve))
				{
					return false;
				}
			}
		}
		return true;
	}

	case VolumeDirection::W:
	{
		for (int u = 0; u < numU; ++u)
		{
			for (int v = 0; v < numV; ++v)
			{
				std::vector<XYZW> wCurveCPs;
				wCurveCPs.reserve(numW);
				for (int w = 0; w < numW; ++w)
				{
					wCurveCPs.push_back(volume.ControlPoints[u][v][w]);
				}

				LN_NurbsCurve curve;
				curve.Degree = volume.DegreeW;
				curve.KnotVector = volume.KnotVectorW;
				curve.ControlPoints = wCurveCPs;

				if (!NurbsCurve::IsClosed(curve))
				{
					return false;
				}
			}
		}
		return true;
	}

	default:
		return false;
	}
}

LNLib::UVW LNLib::NurbsVolume::GetParamOnVolume(const LN_NurbsVolume& volume, const LNLib::XYZ& givenPoint)
{
	std::vector<double> kvU = volume.KnotVectorU;
	std::vector<double> kvV = volume.KnotVectorV;
	std::vector<double> kvW = volume.KnotVectorW;

	double minU = kvU[0], maxU = kvU[kvU.size() - 1];
	double minV = kvV[0], maxV = kvV[kvV.size() - 1];
	double minW = kvW[0], maxW = kvW[kvW.size() - 1];

	bool isClosedU = IsClosed(volume, VolumeDirection::U);
	bool isClosedV = IsClosed(volume, VolumeDirection::V);
	bool isClosedW = IsClosed(volume, VolumeDirection::W);

	int num = 40;
	double stepU = (maxU - minU) / (num - 1);
	double stepV = (maxV - minV) / (num - 1);
	double stepW = (maxW - minW) / (num - 1);

	double minDistSq = std::numeric_limits<double>::max();
	UVW bestParam(minU, minV, minW);

	for (int i = 0; i < num; ++i)
	{
		for (int j = 0; j < num; ++j)
		{
			for (int k = 0; k < num; ++k)
			{
				double currentU = minU + i * stepU;
				double currentV = minV + j * stepV;
				double currentW = minW + k * stepW;

				UVW current(currentU, currentV, currentW);
				LNLib::XYZ point = GetPointOnVolume(volume, current);

				double distSq =(point-givenPoint).SqrLength();
				if (distSq < minDistSq)
				{
					minDistSq = distSq;
					bestParam = current;
				}
			}
		}
	}


	UVW param = bestParam;
	int maxIterations = 20;

	for (int iter = 0; iter < maxIterations; ++iter)
	{
		std::vector<std::vector<std::vector<LNLib::XYZ>>> ders = ComputeVolumeRationalDerivatives(volume, 2, param);

		XYZ V = ders[0][0][0];
		XYZ Vu = ders[1][0][0];
		XYZ Vv = ders[0][1][0];
		XYZ Vw = ders[0][0][1];
		XYZ Vuu = ders[2][0][0];
		XYZ Vvv = ders[0][2][0];
		XYZ Vww = ders[0][0][2];
		XYZ Vuv = ders[1][1][0];
		XYZ Vuw = ders[1][0][1];
		XYZ Vvw = ders[0][1][1];

		XYZ diff = V - givenPoint;
		double currentDist = diff.Length();

		if (MathUtils::IsAlmostEqualTo(currentDist, 0.0, LNLib::Constants::DistanceEpsilon)) {
			return param;
		}

		double J00 = Vu.DotProduct(Vu) + diff.DotProduct(Vuu);
		double J01 = Vu.DotProduct(Vv) + diff.DotProduct(Vuv);
		double J02 = Vu.DotProduct(Vw) + diff.DotProduct(Vuw);

		double J10 = J01;
		double J11 = Vv.DotProduct(Vv) + diff.DotProduct(Vvv);
		double J12 = Vv.DotProduct(Vw) + diff.DotProduct(Vvw);

		double J20 = J02;
		double J21 = J12;
		double J22 = Vw.DotProduct(Vw) + diff.DotProduct(Vww);

		double K0 = Vu.DotProduct(diff);
		double K1 = Vv.DotProduct(diff);
		double K2 = Vw.DotProduct(diff);

		double det = J00 * (J11 * J22 - J12 * J21)
			- J01 * (J10 * J22 - J12 * J20)
			+ J02 * (J10 * J21 - J11 * J20);

		if (MathUtils::IsAlmostEqualTo(std::abs(det), 0.0)) {
			break;
		}

		double deltaU = ((-K0) * (J11 * J22 - J12 * J21)
			- J01 * ((-K1) * J22 - J12 * (-K2))
			+ J02 * ((-K1) * J21 - J11 * (-K2))) / det;

		double deltaV = (J00 * ((-K1) * J22 - J12 * (-K2))
			- (-K0) * (J10 * J22 - J12 * J20)
			+ J02 * (J10 * (-K2) - (-K1) * J20)) / det;

		double deltaW = (J00 * (J11 * (-K2) - (-K1) * J21)
			- J01 * (J10 * (-K2) - (-K1) * J20)
			+ (-K0) * (J10 * J21 - J11 * J20)) / det;


		double newU = param.GetU() + deltaU;
		double newV = param.GetV() + deltaV;
		double newW = param.GetW() + deltaW;

		XYZ deltaWorld = Vu * deltaU + Vv * deltaV + Vw * deltaW;
		if (MathUtils::IsAlmostEqualTo(deltaWorld.Length(), 0.0, LNLib::Constants::DistanceEpsilon)) {
			param = UVW(newU, newV, newW);
			break;
		}

		if (isClosedU) {
			double range = maxU - minU;
			while (newU < minU) newU += range;
			while (newU > maxU) newU -= range;
		}
		else {
			newU = std::max(minU, std::min(maxU, newU));
		}

		if (isClosedV) {
			double range = maxV - minV;
			while (newV < minV) newV += range;
			while (newV > maxV) newV -= range;
		}
		else {
			newV = std::max(minV, std::min(maxV, newV));
		}

		if (isClosedW) {
			double range = maxW - minW;
			while (newW < minW) newW += range;
			while (newW > maxW) newW -= range;
		}
		else {
			newW = std::max(minW, std::min(maxW, newW));
		}

		param = UVW(newU, newV, newW);
	}

	return param;
}

void LNLib::NurbsVolume::Reverse(const LN_NurbsVolume& volume, VolumeDirection direction, LN_NurbsVolume& result)
{
	int degreeU = volume.DegreeU;
	int degreeV = volume.DegreeV;
	int degreeW = volume.DegreeW;

	std::vector<double> kvU = volume.KnotVectorU;
	std::vector<double> kvV = volume.KnotVectorV;
	std::vector<double> kvW = volume.KnotVectorW;

	std::vector<std::vector<std::vector<LNLib::XYZW>>> controlPoints = volume.ControlPoints;

	int n = controlPoints.size() - 1;
	int m = controlPoints[0].size() - 1;
	int l = controlPoints[0][0].size() - 1;

	std::vector<std::vector<std::vector<LNLib::XYZW>>> newControlPoints;

	result.DegreeU = degreeU;
	result.DegreeV = degreeV;
	result.DegreeW = degreeW;

	if (direction == VolumeDirection::U)
	{
		std::vector<double> kvS = kvV;
		std::vector<double> kvT = kvW;

		newControlPoints.resize(n + 1, std::vector<std::vector<LNLib::XYZW>>(m + 1, std::vector<LNLib::XYZW>(l + 1)));

		std::vector<double> kvR;

		for (int sep = 0; sep <= l; sep++)
		{
			for (int col = 0; col <= m; col++)
			{
				std::vector<LNLib::XYZW> uControlPoints(n + 1);
				for (int i = 0; i <= n; ++i) {
					uControlPoints[i] = controlPoints[i][col][sep];
				}

				LNLib::LN_NurbsCurve curve;
				curve.Degree = degreeU;
				curve.KnotVector = kvU;
				curve.ControlPoints = uControlPoints;

				LNLib::LN_NurbsCurve resultCurve;
				LNLib::NurbsCurve::Reverse(curve, resultCurve);

				kvR = resultCurve.KnotVector;
				std::vector<LNLib::XYZW> iControlPoints = resultCurve.ControlPoints;

				for (int i = 0; i < iControlPoints.size(); i++)
				{
					newControlPoints[i][col][sep] = iControlPoints[i];
				}
			}
		}

		result.KnotVectorU = kvR;
		result.KnotVectorV = kvS;
		result.KnotVectorW = kvT;
		result.ControlPoints = newControlPoints;
	}
	else if (direction == VolumeDirection::V)
	{
		std::vector<double> kvR = kvU;
		std::vector<double> kvT = kvW;

		std::vector<double> kvS;

		for (int sep = 0; sep <= l; sep++)
		{
			for (int row = 0; row <= n; row++)
			{
				std::vector<LNLib::XYZW> vControlPoints(m + 1);
				for (int j = 0; j <= m; ++j) {
					vControlPoints[j] = controlPoints[row][j][sep];
				}

				LNLib::LN_NurbsCurve curve;
				curve.Degree = degreeV;
				curve.KnotVector = kvV;
				curve.ControlPoints = vControlPoints;

				LNLib::LN_NurbsCurve resultCurve;
				LNLib::NurbsCurve::Reverse(curve, resultCurve);

				kvS = resultCurve.KnotVector;
				std::vector<LNLib::XYZW> iControlPoints = resultCurve.ControlPoints;

				for (int j = 0; j < iControlPoints.size(); j++)
				{
					newControlPoints[row][j][sep] = iControlPoints[j];
				}
			}
		}

		result.KnotVectorU = kvR;
		result.KnotVectorV = kvS;
		result.KnotVectorW = kvT;
		result.ControlPoints = newControlPoints;
	}
	else if (direction == VolumeDirection::W)
	{
		std::vector<double> kvR = kvU;
		std::vector<double> kvS = kvV;

		std::vector<double> kvT;

		for (int col = 0; col <= m; col++)
		{
			for (int row = 0; row <= n; row++)
			{
				std::vector<LNLib::XYZW> wControlPoints(l + 1);
				for (int k = 0; k <= l; ++k) {
					wControlPoints[k] = controlPoints[row][col][k];
				}

				LNLib::LN_NurbsCurve curve;
				curve.Degree = degreeW;
				curve.KnotVector = kvW;
				curve.ControlPoints = wControlPoints;

				LNLib::LN_NurbsCurve resultCurve;
				LNLib::NurbsCurve::Reverse(curve, resultCurve);

				kvT = resultCurve.KnotVector;
				std::vector<LNLib::XYZW> iControlPoints = resultCurve.ControlPoints;

				for (int k = 0; k < iControlPoints.size(); k++)
				{
					newControlPoints[row][col][k] = iControlPoints[k];
				}
			}
		}

		result.KnotVectorU = kvR;
		result.KnotVectorV = kvS;
		result.KnotVectorW = kvT;
		result.ControlPoints = newControlPoints;
	}
}

void LNLib::NurbsVolume::Swap(const LN_NurbsVolume& volume, VolumeDirection direction, VolumeDirection swap, LN_NurbsVolume& result)
{
	int degreeU = volume.DegreeU;
	int degreeV = volume.DegreeV;
	int degreeW = volume.DegreeW;

	std::vector<double> kvU = volume.KnotVectorU;
	std::vector<double> kvV = volume.KnotVectorV;
	std::vector<double> kvW = volume.KnotVectorW;

	std::vector<std::vector<std::vector<LNLib::XYZW>>> controlPoints = volume.ControlPoints;
	std::vector<std::vector<std::vector<LNLib::XYZW>>> newContrlPoints;

	int u = controlPoints.size();
	int v = controlPoints[0].size();
	int w = controlPoints[0][0].size();

	if (direction == VolumeDirection::U && swap == VolumeDirection::V)
	{
		result.DegreeU = degreeV;
		result.DegreeV = degreeU;
		result.DegreeW = degreeW;

		result.KnotVectorU = kvV;
		result.KnotVectorV = kvU;
		result.KnotVectorW = kvW;

		std::vector<std::vector<std::vector<LNLib::XYZW>>> newContrlPoints(v);

		for (int j = 0; j < v; ++j) {
			newContrlPoints[j].resize(u);

			for (int i = 0; i < u; ++i) {
				newContrlPoints[j][i].resize(w);

				for (int k = 0; k < w; ++k) {
					newContrlPoints[j][i][k] = controlPoints[i][j][k];
				}
			}
		}

		result.ControlPoints = newContrlPoints;
	}
	else if (direction == VolumeDirection::U && swap == VolumeDirection::W)
	{
		result.DegreeU = degreeW;
		result.DegreeV = degreeV;
		result.DegreeW = degreeU;

		result.KnotVectorU = kvW;
		result.KnotVectorV = kvV;
		result.KnotVectorW = kvU;

		std::vector<std::vector<std::vector<LNLib::XYZW>>> newContrlPoints(w);

		for (int k = 0; k < w; ++k) {
			newContrlPoints[k].resize(v);
			for (int j = 0; j < v; ++j) {
				newContrlPoints[k][j].resize(u);
				for (int i = 0; i < u; ++i) {
					newContrlPoints[k][j][i] = controlPoints[i][j][k];
				}
			}
		}

		result.ControlPoints = newContrlPoints;
	}
	else if (direction == VolumeDirection::V && swap == VolumeDirection::W)
	{
		result.DegreeU = degreeU;
		result.DegreeV = degreeW;
		result.DegreeW = degreeV;

		result.KnotVectorU = kvU;
		result.KnotVectorV = kvW;
		result.KnotVectorW = kvV;

		std::vector<std::vector<std::vector<LNLib::XYZW>>> newContrlPoints(u);

		for (int i = 0; i < u; ++i) {
			newContrlPoints[i].resize(w);

			for (int k = 0; k < w; ++k) {
				newContrlPoints[i][k].resize(v);

				for (int j = 0; j < v; ++j) {
					newContrlPoints[i][k][j] = controlPoints[i][j][k];
				}
			}
		}

		result.ControlPoints = newContrlPoints;
	}
}

bool LNLib::NurbsVolume::Check(const LN_NurbsVolume& volume)
{
	std::vector<double> kvU = volume.KnotVectorU;
	std::vector<double> kvV = volume.KnotVectorV;
	std::vector<double> kvW = volume.KnotVectorW;

	double minU = kvU[0], maxU = kvU[kvU.size() - 1];
	double minV = kvV[0], maxV = kvV[kvV.size() - 1];
	double minW = kvW[0], maxW = kvW[kvW.size() - 1];

	int num = 50;

	double stepU = (maxU - minU) / (num - 1);
	double stepV = (maxV - minV) / (num - 1);
	double stepW = (maxW - minW) / (num - 1);

	int invalidCount = 0;

	for (int i = 0; i < num; ++i)
	{
		for (int j = 0; j < num; ++j)
		{
			for (int k = 0; k < num; ++k)
			{
				double u = minU + i * stepU;
				double v = minV + j * stepV;
				double w = minW + k * stepW;
				UVW param(u, v, w);

				std::vector<std::vector<std::vector<XYZ>>> ders = ComputeVolumeRationalDerivatives(volume, 1, param);

				XYZ Vu = ders[1][0][0];
				XYZ Vv = ders[0][1][0];
				XYZ Vw = ders[0][0][1];

				double jacobian = Vu.X() * Vv.Y() * Vw.Z() +
					Vu.Y() * Vv.Z() * Vw.X() +
					Vu.Z() * Vv.X() * Vw.Y() -
					Vu.Z() * Vv.Y() * Vw.X() -
					Vu.Y() * Vv.X() * Vw.Z() -
					Vu.X() * Vv.Z() * Vw.Y();

				if (jacobian <= 0.0)
				{
					return false;
				}
			}
		}
	}
	return true;
}
