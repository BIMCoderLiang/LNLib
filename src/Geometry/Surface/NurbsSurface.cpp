#include "NurbsSurface.h"
#include "Polynomials.h"
#include "UV.h"
#include "XYZ.h"
#include "XYZW.h"
#include "Matrix.h"
#include "MathUtils.h"
#include "NurbsCurve.h"
#include "BsplineSurface.h"

void LNLib::NurbsSurface::GetPointOnSurface(const std::vector<std::vector<XYZW>>& controlPoints, const std::vector<double>& knotVectorU, const std::vector<double>& knotVectorV, unsigned int degreeU, unsigned int degreeV, UV uv, XYZ& point)
{
	unsigned int n = static_cast<int>(knotVectorU.size() - degreeU - 2);
	unsigned int uSpanIndex = Polynomials::GetKnotSpanIndex(n, degreeU, uv.GetU(), knotVectorU);
	std::vector<double> basisFunctionsU;
	basisFunctionsU.resize(degreeU + 1);
	Polynomials::BasisFunctions(uSpanIndex, degreeU, uv.GetU(), knotVectorU, basisFunctionsU);

	unsigned int m = static_cast<int>(knotVectorV.size() - degreeV - 2);
	unsigned int vSpanIndex = Polynomials::GetKnotSpanIndex(m, degreeV, uv.GetV(), knotVectorV);
	std::vector<double> basisFunctionsV;
	basisFunctionsV.resize(degreeV + 1);
	Polynomials::BasisFunctions(vSpanIndex, degreeV, uv.GetV(), knotVectorV, basisFunctionsV);


	std::vector<XYZW> temp;
	for (unsigned int l = 0; l <= degreeV; l++)
	{
		temp[l] = XYZW(0, 0, 0, 0);
		for (unsigned int k = 0; k <= degreeU; k++)
		{
			temp[l] += basisFunctionsU[k] * controlPoints[uSpanIndex - degreeU + k][vSpanIndex - degreeV + l];
		}
	}

	XYZW Sw = XYZW(0, 0, 0, 0);
	for (unsigned int l = 0; l <= degreeV; l++)
	{
		Sw += basisFunctionsV[l] * temp[l];
	}
	point = Sw.ToXYZ(true);
}


void LNLib::NurbsSurface::ComputeRationalSurfaceDerivatives(const std::vector<std::vector<XYZW>>& controlPoints, const std::vector<double>& knotVectorU, const std::vector<double>& knotVectorV, unsigned int degreeU, unsigned int degreeV, UV uv, unsigned int derivative, std::vector<std::vector<XYZ>>& derivatives)
{
	
	std::vector<std::vector<XYZW>> ders;
	ders.resize(degreeU + 1);
	for (unsigned int i = 0; i <= degreeU; i++)
	{
		ders[i].resize(degreeV + 1);
	}
	BsplineSurface::ComputeDerivatives(controlPoints, knotVectorU, knotVectorV, degreeU, degreeV, uv, derivative, ders);

	std::vector<std::vector<XYZ>> Aders;
	for (unsigned int i = 0; i <= degreeU; i++)
	{
		for (unsigned int j = 0; j <= degreeV; j++)
		{
			Aders[i][j] = ders[i][j].ToXYZ(false);
		}
	}

	std::vector<std::vector<double>> wders;
	for (unsigned int i = 0; i <= degreeU; i++)
	{
		for (unsigned int j = 0; j <= degreeV; j++)
		{
			wders[i][j] = ders[i][j].GetW();
		}
	}

	for (unsigned int k = 0; k <= derivative; k++)
	{
		for (unsigned int l = 0; l <= derivative - k; l++)
		{
			XYZ v = Aders[k][l];
			for (unsigned int j = 1; j <= l; j++)
			{
				v = v - MathUtils::Binomial(l, j) * wders[0][j] * derivatives[k][l - j];
			}

			for (unsigned int i = 1; i <= k; i++)
			{
				v = v - MathUtils::Binomial(k, i) * wders[i][0] * derivatives[k - i][l];

				XYZ v2 = XYZ();
				for (unsigned int j = 1; j <= l; j++)
				{
					v2 = v2 + MathUtils::Binomial(l, j) * wders[i][j] * derivatives[k - i][l - j];
				}

				v = v - MathUtils::Binomial(k, i) * v2;
			}

			derivatives[k][l] = v / wders[0][0];
		}
	}
}

void LNLib::NurbsSurface::InsertKnot(const std::vector<std::vector<XYZW>>& controlPoints, const std::vector<double>& knotVector, unsigned int degree, double insertKnot, unsigned int times, bool isUDirection, std::vector<double>& insertedKnotVector, std::vector<std::vector<XYZW>>& updatedControlPoints)
{
	int n = static_cast<int>(knotVector.size()) - degree - 2;
	int knotSpanIndex = Polynomials::GetKnotSpanIndex(n, degree, insertKnot, knotVector);
	unsigned int multiplicity = Polynomials::GetKnotMultiplicity(degree, insertKnot, knotVector);

	if (multiplicity == degree)
	{
		insertedKnotVector = knotVector;
		updatedControlPoints = controlPoints;
	}

	if ((times + multiplicity) > degree)
	{
		times = degree - multiplicity;
	}

	for (int i = 0; i <= knotSpanIndex; i++)
	{
		insertedKnotVector[i] = knotVector[i];
	}

	for (unsigned int i = 1; i <= times; i++)
	{
		insertedKnotVector[knotSpanIndex + i] = insertKnot;
	}

	for (int i = knotSpanIndex + 1; i < knotVector.size(); i++)
	{
		insertedKnotVector[i + times] = knotVector[i];
	}

	std::vector<std::vector<double>> alpha;
	alpha.resize(degree - multiplicity);
	for (unsigned int i = 0; i < degree - multiplicity; i++)
	{
		alpha[i].resize(times + 1);
	}

	for (unsigned int j = 1; j <= times; j++)
	{
		int L = knotSpanIndex - degree + j;
		for (unsigned int i = 0; i <= degree - j - multiplicity; i++)
		{
			alpha[i][j] = (insertKnot - knotVector[L + i]) / (knotVector[i + knotSpanIndex + 1] - knotVector[L + i]);
		}
	}

	

	std::vector<XYZW> temp;
	temp.resize(degree + 1);

	if (isUDirection)
	{
		updatedControlPoints.resize(controlPoints.size() + times);
		for (int i = 0; i < controlPoints.size() + times; i++)
		{
			updatedControlPoints[i].resize(controlPoints[i].size());
		}

		for (int col = 0; col < controlPoints[0].size(); col++)
		{
			for (unsigned int i = 0; i <= knotSpanIndex - degree; i++)
			{
				updatedControlPoints[i][col] = controlPoints[i][col];
			}

			for (int i = knotSpanIndex - multiplicity; i < controlPoints.size(); i++)
			{
				updatedControlPoints[i + times][col] = controlPoints[i][col];
			}

			for (unsigned int i = 0; i < degree - multiplicity + 1; i++)
			{
				temp[i] = controlPoints[knotSpanIndex - degree + i][col];
			}


			int L = 0;
			for (unsigned int j = 1; j <= times; j++)
			{
				int L = knotSpanIndex - degree + j;
				for (unsigned int i = 0; i <= degree - j - multiplicity; i++)
				{
					double a = alpha[i][j];
					temp[i] = a * temp[i + 1] + (1.0 - a) * temp[i];
				}

				updatedControlPoints[L][col] = temp[0];
				updatedControlPoints[knotSpanIndex + times - j - multiplicity][col] = temp[degree - j - multiplicity];
			}

			for (unsigned int i = L + 1; i < knotSpanIndex - multiplicity; i++)
			{
				updatedControlPoints[i][col] = temp[i - L];
			}
		}
	}
	else
	{
		updatedControlPoints.resize(controlPoints.size());
		for (int i = 0; i < controlPoints.size(); i++)
		{
			updatedControlPoints[i].resize(controlPoints[i].size() + times);
		}

		for (int row = 0; row < controlPoints.size(); row++)
		{
			for (unsigned int i = 0; i <= knotSpanIndex - degree; i++)
			{
				updatedControlPoints[row][i] = controlPoints[row][i];
			}

			for (int i = knotSpanIndex - multiplicity; i < controlPoints.size(); i++)
			{
				updatedControlPoints[row][i+times] = controlPoints[row][i];
			}

			for (unsigned int i = 0; i < degree - multiplicity + 1; i++)
			{
				temp[i] = controlPoints[row][knotSpanIndex - degree + i];
			}


			int L = 0;
			for (unsigned int j = 1; j <= times; j++)
			{
				int L = knotSpanIndex - degree + j;
				for (unsigned int i = 0; i <= degree - j - multiplicity; i++)
				{
					double a = alpha[i][j];
					temp[i] = a * temp[i + 1] + (1.0 - a) * temp[i];
				}

				updatedControlPoints[row][L] = temp[0];
				updatedControlPoints[row][knotSpanIndex + times - j - multiplicity] = temp[degree - j - multiplicity];
			}

			for (unsigned int i = L + 1; i < knotSpanIndex - multiplicity; i++)
			{
				updatedControlPoints[row][i] = temp[i - L];
			}
		}
	}
}

void LNLib::NurbsSurface::RefineKnotVector(const std::vector<std::vector<XYZW>>& controlPoints, const std::vector<double>& knotVectorU, const std::vector<double>& knotVectorV, unsigned int degreeU, unsigned int degreeV, std::vector<double>& insertKnotElements, bool isUDirection, std::vector<double>& insertedKnotVectorU, std::vector<double>& insertedKnotVectorV, std::vector<std::vector<XYZW>>& updatedControlPoints)
{
	
	if (isUDirection)
	{
		std::vector<std::vector<XYZW>> transposedControlPoints = Matrix::Transpose(controlPoints);
		for (int i = 0; i < static_cast<int>(transposedControlPoints.size()); i++)
		{
			std::vector<XYZW> temp;
			NurbsCurve::RefineKnotVector(degreeU, knotVectorU, transposedControlPoints[i], insertKnotElements, insertedKnotVectorU, temp);
			updatedControlPoints[i] = temp;
		}
		for (int i = 0; i < static_cast<int>(knotVectorV.size()); i++)
		{
			insertedKnotVectorV[i] = knotVectorV[i];
		}
		updatedControlPoints = Matrix::Transpose(updatedControlPoints);
	}
	else
	{
		for (int i = 0; i < static_cast<int>(controlPoints.size()); i++)
		{
			NurbsCurve::RefineKnotVector(degreeV, knotVectorV, controlPoints[i], insertKnotElements, insertedKnotVectorV, updatedControlPoints[i]);
		}
		for (int i = 0; i < static_cast<int>(knotVectorU.size()); i++)
		{
			insertedKnotVectorU[i] = knotVectorU[i];
		}
	}
}

void LNLib::NurbsSurface::ToBezierPatches(const std::vector<std::vector<XYZW>>& controlPoints, const std::vector<double>& knotVectorU, const std::vector<double>& knotVectorV, unsigned int degreeU, unsigned int degreeV, bool isUDirection, int& bezierCount, std::vector<std::vector<std::vector<XYZW>>>& decomposedControlPoints)
{
	if (isUDirection)
	{
		int n = static_cast<int>(knotVectorU.size() - degreeU - 2);
		int m = n + degreeU + 1;

		int a = static_cast<int>(degreeU);
		int b = static_cast<int>(degreeU) + 1;

		bezierCount = 0;
		for (unsigned i = 0; i <= degreeU; i++)
		{
			for (int row = 0; row <= m; row++)
			{
				decomposedControlPoints[bezierCount][i][row] = controlPoints[i][row];
			}
		}

		while (b < m)
		{
			int i = b;
			while (b < m && MathUtils::IsAlmostEqualTo(knotVectorU[b + 1], knotVectorU[b]))
			{
				b++;
			}
			int multi = b - i + 1;
			if (multi < static_cast<int>(degreeU))
			{
				double numerator = knotVectorU[b] - knotVectorU[a];

				std::vector<double> alhpaVector;
				alhpaVector.resize(degreeU - multi);
				for (int j = degreeU; j > multi; j--)
				{
					alhpaVector[j - multi - 1] = numerator / (knotVectorU[a + j] - knotVectorU[a]);
				}

				int r = degreeU - multi;
				for (int j = 1; j <= r; j++)
				{
					int save = r - j;
					int s = multi + j;
					for (int k = degreeU; k >= s; k--)
					{
						double alpha = alhpaVector[k - s];
						for (int row = 0; row <= m; row++)
						{
							decomposedControlPoints[bezierCount][k][row] = alpha * decomposedControlPoints[bezierCount][k][row] + 
																			(1.0 - alpha) * decomposedControlPoints[bezierCount][k - 1][row];
						}
						
					}

					if (b < m)
					{
						for (int row = 0; row <= m; row++)
						{
							decomposedControlPoints[bezierCount + 1][save][row] = decomposedControlPoints[bezierCount][degreeU][row];
						}
						
					}
				}

				bezierCount += 1;
				if (b < m)
				{
					for (unsigned int i = degreeU - multi; i <= degreeU; i++)
					{
						for (int row = 0; row <= m; row++)
						{
							decomposedControlPoints[bezierCount][i][row] = controlPoints[b - degreeU + i][row];
						}
					}

					a = b;
					b += 1;
				}
			}
		}
	}
	else
	{
		int n = static_cast<int>(knotVectorV.size() - degreeV - 2);
		int m = n + degreeV + 1;

		int a = static_cast<int>(degreeV);
		int b = static_cast<int>(degreeV) + 1;

		bezierCount = 0;

		std::vector<std::vector<XYZW>> transposedControlPoints = Matrix::Transpose(controlPoints);
		for (unsigned i = 0; i <= degreeV; i++)
		{
			for (int row = 0; row <= m; row++)
			{
				decomposedControlPoints[bezierCount][i][row] = transposedControlPoints[i][row];
			}
		}

		while (b < m)
		{
			int i = b;
			while (b < m && MathUtils::IsAlmostEqualTo(knotVectorV[b + 1], knotVectorV[b]))
			{
				b++;
			}
			int multi = b - i + 1;
			if (multi < static_cast<int>(degreeV))
			{
				double numerator = knotVectorV[b] - knotVectorV[a];

				std::vector<double> alhpaVector;
				alhpaVector.resize(degreeV - multi);
				for (int j = degreeV; j > multi; j--)
				{
					alhpaVector[j - multi - 1] = numerator / (knotVectorV[a + j] - knotVectorV[a]);
				}

				int r = degreeV - multi;
				for (int j = 1; j <= r; j++)
				{
					int save = r - j;
					int s = multi + j;
					for (int k = degreeV; k >= s; k--)
					{
						double alpha = alhpaVector[k - s];
						for (int row = 0; row <= m; row++)
						{
							decomposedControlPoints[bezierCount][k][row] = alpha * decomposedControlPoints[bezierCount][k][row] +
								(1.0 - alpha) * decomposedControlPoints[bezierCount][k - 1][row];
						}

					}

					if (b < m)
					{
						for (int row = 0; row <= m; row++)
						{
							decomposedControlPoints[bezierCount + 1][save][row] = decomposedControlPoints[bezierCount][degreeV][row];
						}

					}
				}

				bezierCount += 1;
				if (b < m)
				{
					for (unsigned int i = degreeV - multi; i <= degreeV; i++)
					{
						for (int row = 0; row <= m; row++)
						{
							decomposedControlPoints[bezierCount][i][row] = transposedControlPoints[b - degreeV + i][row];
						}
					}

					a = b;
					b += 1;
				}
			}
		}
	}
}


