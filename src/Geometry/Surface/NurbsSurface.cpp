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
	XYZW sw = XYZW();
	BsplineSurface::GetPointOnSurface(controlPoints, knotVectorU, knotVectorV, degreeU, degreeV, uv, sw);
	point = sw.ToXYZ(true);
}


void LNLib::NurbsSurface::ComputeRationalSurfaceDerivatives(const std::vector<std::vector<XYZW>>& controlPoints, const std::vector<double>& knotVectorU, const std::vector<double>& knotVectorV, unsigned int degreeU, unsigned int degreeV, UV uv, unsigned int derivative, std::vector<std::vector<XYZ>>& derivatives)
{
	
	std::vector<std::vector<XYZW>> ders;
	BsplineSurface::ComputeDerivatives(controlPoints, knotVectorU, knotVectorV, degreeU, degreeV, uv, derivative, ders);

	std::vector<std::vector<XYZ>> Aders;
	for (int i = 0; i <= static_cast<int>(degreeU); i++)
	{
		for (int j = 0; j <= static_cast<int>(degreeV); j++)
		{
			Aders[i][j] = ders[i][j].ToXYZ(false);
		}
	}

	std::vector<std::vector<double>> wders;
	for (int i = 0; i <= static_cast<int>(degreeU); i++)
	{
		for (int j = 0; j <= static_cast<int>(degreeV); j++)
		{
			wders[i][j] = ders[i][j].GetW();
		}
	}

	for (int k = 0; k <= static_cast<int>(derivative); k++)
	{
		for (int l = 0; l <= static_cast<int>(derivative - k); l++)
		{
			XYZ v = Aders[k][l];
			for ( int j = 1; j <= l; j++)
			{
				v = v - MathUtils::Binomial(l, j) * wders[0][j] * derivatives[k][l - j];
			}

			for (int i = 1; i <= k; i++)
			{
				v = v - MathUtils::Binomial(k, i) * wders[i][0] * derivatives[k - i][l];

				XYZ v2 = XYZ();
				for (int j = 1; j <= l; j++)
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
	int multiplicity = Polynomials::GetKnotMultiplicity(insertKnot, knotVector);

	if (multiplicity == degree)
	{
		insertedKnotVector = knotVector;
		updatedControlPoints = controlPoints;
		return;
	}

	if ((times + multiplicity) > degree)
	{
		times = degree - multiplicity;
	}

	insertedKnotVector.resize(knotVector.size() + times);
	for (int i = 0; i <= knotSpanIndex; i++)
	{
		insertedKnotVector[i] = knotVector[i];
	}

	for (int i = 1; i <= static_cast<int>(times); i++)
	{
		insertedKnotVector[knotSpanIndex + i] = insertKnot;
	}

	for (int i = knotSpanIndex + 1; i < static_cast<int>(knotVector.size()); i++)
	{
		insertedKnotVector[i + times] = knotVector[i];
	}

	std::vector<std::vector<double>> alpha;
	alpha.resize(degree - multiplicity);
	for (int i = 0; i < static_cast<int>(degree) - multiplicity; i++)
	{
		alpha[i].resize(times + 1);
	}

	for (int j = 1; j <= static_cast<int>(times); j++)
	{
		int L = knotSpanIndex - degree + j;
		for (int i = 0; i <= static_cast<int>(degree) - j - multiplicity; i++)
		{
			alpha[i][j] = (insertKnot - knotVector[L + i]) / (knotVector[i + knotSpanIndex + 1] - knotVector[L + i]);
		}
	}

	std::vector<XYZW> temp;
	temp.resize(degree + 1);

	int controlPointsRows = static_cast<int>(controlPoints.size());
	int controlPointsColumns = static_cast<int>(controlPoints[0].size());

	if (isUDirection)
	{
		updatedControlPoints.resize(controlPointsRows + times);
		for (int i = 0; i < controlPointsRows + static_cast<int>(times); i++)
		{
			updatedControlPoints[i].resize(controlPointsColumns);
		}

		for (int col = 0; col < controlPointsColumns; col++)
		{
			for (int i = 0; i <= knotSpanIndex - static_cast<int>(degree); i++)
			{
				updatedControlPoints[i][col] = controlPoints[i][col];
			}

			for (int i = knotSpanIndex - multiplicity; i < controlPointsRows; i++)
			{
				updatedControlPoints[i + times][col] = controlPoints[i][col];
			}

			for (int i = 0; i < static_cast<int>(degree) - multiplicity + 1; i++)
			{
				temp[i] = controlPoints[knotSpanIndex - degree + i][col];
			}


			int L = 0;
			for (int j = 1; j <= static_cast<int>(times); j++)
			{
				L = knotSpanIndex - degree + j;
				for (int i = 0; i <= static_cast<int>(degree) - j - multiplicity; i++)
				{
					double a = alpha[i][j];
					temp[i] = a * temp[i + 1] + (1.0 - a) * temp[i];
				}

				updatedControlPoints[L][col] = temp[0];
				updatedControlPoints[knotSpanIndex + times - j - multiplicity][col] = temp[degree - j - multiplicity];
			}

			for (int i = L + 1; i < knotSpanIndex - multiplicity; i++)
			{
				updatedControlPoints[i][col] = temp[i - L];
			}
		}
	}
	else
	{
		updatedControlPoints.resize(controlPointsRows);
		for (int i = 0; i < controlPointsRows; i++)
		{
			updatedControlPoints[i].resize(controlPointsColumns + times);
		}

		for (int row = 0; row < controlPointsRows; row++)
		{
			for (int i = 0; i <= knotSpanIndex - static_cast<int>(degree); i++)
			{
				updatedControlPoints[row][i] = controlPoints[row][i];
			}

			for (int i = knotSpanIndex - multiplicity; i < controlPointsColumns; i++)
			{
				updatedControlPoints[row][i+times] = controlPoints[row][i];
			}

			for (int i = 0; i < static_cast<int>(degree) - multiplicity + 1; i++)
			{
				temp[i] = controlPoints[row][knotSpanIndex - degree + i];
			}


			int L = 0;
			for (int j = 1; j <= static_cast<int>(times); j++)
			{
				L = knotSpanIndex - degree + j;
				for (int i = 0; i <= static_cast<int>(degree) - j - multiplicity; i++)
				{
					double a = alpha[i][j];
					temp[i] = a * temp[i + 1] + (1.0 - a) * temp[i];
				}

				updatedControlPoints[row][L] = temp[0];
				updatedControlPoints[row][knotSpanIndex + times - j - multiplicity] = temp[degree - j - multiplicity];
			}

			for (int i = L + 1; i < knotSpanIndex - multiplicity; i++)
			{
				updatedControlPoints[row][i] = temp[i - L];
			}
		}
	}
}

void LNLib::NurbsSurface::RefineKnotVector(const std::vector<std::vector<XYZW>>& controlPoints, const std::vector<double>& knotVectorU, const std::vector<double>& knotVectorV, unsigned int degreeU, unsigned int degreeV, std::vector<double>& insertKnotElements, bool isUDirection, std::vector<double>& insertedKnotVectorU, std::vector<double>& insertedKnotVectorV, std::vector<std::vector<XYZW>>& updatedControlPoints)
{
	std::vector<double> knots;
	unsigned int degree;
	std::vector<std::vector<XYZW>> tempControlPoints;
	std::vector<std::vector<XYZW>> newControlPoints;

	if (isUDirection)
	{
		Matrix::Transpose(controlPoints, tempControlPoints);
		int size = static_cast<int>(knotVectorU.size());
		knots.resize(size);
		for (int i = 0; i < size; i++)
		{
			knots[i] = knotVectorU[i];
		}
		degree = degreeU;
	}
	else
	{
		tempControlPoints = controlPoints;
		int size = static_cast<int>(knotVectorV.size());
		knots.resize(size);
		for (int i = 0; i < size; i++)
		{
			knots[i] = knotVectorV[i];
		}
		degree = degreeV;
	}

	std::vector<double> insertedKnotVector;
	
	for (int i = 0; i < static_cast<int>(tempControlPoints.size()); i++)
	{
		insertedKnotVector.clear();
		std::vector<XYZW> updatedCurveControlPoints;
		NurbsCurve::RefineKnotVector(degree, knots, tempControlPoints[i], insertKnotElements, insertedKnotVector, updatedCurveControlPoints);
		newControlPoints[i] = updatedCurveControlPoints;
	}

	if (isUDirection)
	{
		insertedKnotVectorU = insertedKnotVector;
		insertedKnotVectorV = knotVectorV;
		Matrix::Transpose(newControlPoints, updatedControlPoints);
	}
	else
	{
		insertedKnotVectorU = knotVectorU;
		insertedKnotVectorV = insertedKnotVector;
		updatedControlPoints = newControlPoints;
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

		std::vector<std::vector<XYZW>> transposedControlPoints; 
		Matrix::Transpose(controlPoints, transposedControlPoints);
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

void LNLib::NurbsSurface::ElevateDegree(const std::vector<std::vector<XYZW>>& controlPoints, const std::vector<double>& knotVectorU, const std::vector<double>& knotVectorV, unsigned int degreeU, unsigned int degreeV, unsigned int times, bool isUDirection, std::vector<double>& insertedKnotVectorU, std::vector<double>& insertedKnotVectorV, std::vector<std::vector<XYZW>>& updatedControlPoints)
{
	std::vector<double> knots;
	unsigned int degree;
	std::vector<std::vector<XYZW>> tempControlPoints;
	std::vector<std::vector<XYZW>> newControlPoints;

	if (isUDirection)
	{
		Matrix::Transpose(controlPoints, tempControlPoints);
		int size = static_cast<int>(knotVectorU.size());
		knots.resize(size);
		for (int i = 0; i < size; i++)
		{
			knots[i] = knotVectorU[i];
		}
		degree = degreeU;
	}
	else
	{
		tempControlPoints = controlPoints;
		int size = static_cast<int>(knotVectorV.size());
		knots.resize(size);
		for (int i = 0; i < size; i++)
		{
			knots[i] = knotVectorV[i];
		}
		degree = degreeV;
	}

	std::vector<double> insertedKnotVector;

	for (int i = 0; i < static_cast<int>(tempControlPoints.size()); i++)
	{
		insertedKnotVector.clear();
		std::vector<XYZW> updatedCurveControlPoints;
		NurbsCurve::ElevateDegree(degree, knots, tempControlPoints[i], times, insertedKnotVector, updatedCurveControlPoints);
		newControlPoints[i] = updatedCurveControlPoints;
	}

	if (isUDirection)
	{
		insertedKnotVectorU = insertedKnotVector;
		insertedKnotVectorV = knotVectorV;
		Matrix::Transpose(newControlPoints, updatedControlPoints);
	}
	else
	{
		insertedKnotVectorU = knotVectorU;
		insertedKnotVectorV = insertedKnotVector;
		updatedControlPoints = newControlPoints;
	}
}


