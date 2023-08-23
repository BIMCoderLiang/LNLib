#include "gtest/gtest.h"
#include "Polynomials.h"
#include "MathUtils.h"

using namespace LNLib;

TEST(Test_Polynomials, All)
{
	int degree = 2;
	int error_degree = -4;
	int error_degree1 = 4;
	std::vector<double> coefficients = { 1,2,3 };
	EXPECT_THROW(Polynomials::Horner(error_degree, coefficients, 0.5), std::invalid_argument);
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(Polynomials::Horner(degree, coefficients, 0.5),2.75));
	EXPECT_THROW(Polynomials::Horner(error_degree1, coefficients, 0.5),std::invalid_argument);

	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(Polynomials::Bernstein(-1, degree, 0.5), 0));
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(Polynomials::Bernstein(degree + 1, degree, 0.5), 0));
	EXPECT_THROW(Polynomials::Bernstein(1, degree, 1.5), std::out_of_range);
	EXPECT_THROW(Polynomials::Bernstein(1, degree, -2), std::out_of_range);
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(Polynomials::Bernstein(0, degree, 0.5),1));
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(Polynomials::Bernstein(degree, degree, 0.5), 1));

	EXPECT_THROW(Polynomials::AllBernstein(-1,0.5), std::invalid_argument);
	EXPECT_THROW(Polynomials::AllBernstein(degree, 1.5), std::out_of_range);
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(Polynomials::AllBernstein(degree, 0.5)[0], 0.25));
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(Polynomials::AllBernstein(degree, 0.5)[degree], 0.25));

	std::vector<double> error_knotVector = { 0,0,1,2,1,0 };
	EXPECT_THROW(Polynomials::GetKnotMultiplicity(error_knotVector, 1.0),std::invalid_argument);
	std::vector<double> knotVector = { 0,0,0,1,1,1 };
	EXPECT_TRUE(Polynomials::GetKnotMultiplicity(knotVector, 1.0) == 3);
	EXPECT_TRUE(Polynomials::GetKnotSpanIndex(2, knotVector, 1) == 2);

	knotVector = { 0,0,0,1,2,3,4,4,5,5,5 };
	int spanIndex = Polynomials::GetKnotSpanIndex(degree, knotVector, 5.0 / 2);
	std::vector<double> basis = Polynomials::BasisFunctions(spanIndex, degree, knotVector, 5.0 / 2);
	std::vector<double> check = {1.0/8, 6.0/8, 1.0/8};
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(basis[0], check[0]) && MathUtils::IsAlmostEqualTo(basis[1], check[1]) && MathUtils::IsAlmostEqualTo(basis[2], check[2]));

	auto ders = Polynomials::BasisFunctionsDerivatives(spanIndex, degree, 2, knotVector, 5.0 / 2);
	std::vector<std::vector<double>> checkders = { {0.125,0.75,0.125},{-0.5,0.0,0.5},{1.0,-2.0,1.0} };
	EXPECT_TRUE(MathUtils::IsAlmostEqualTo(ders[0][0], checkders[0][0]) && MathUtils::IsAlmostEqualTo(ders[0][1], checkders[0][1]) && MathUtils::IsAlmostEqualTo(ders[0][2], checkders[0][2]) &&
			    MathUtils::IsAlmostEqualTo(ders[1][0], checkders[1][0]) && MathUtils::IsAlmostEqualTo(ders[1][1], checkders[1][1]) && MathUtils::IsAlmostEqualTo(ders[1][2], checkders[1][2]) &&
		        MathUtils::IsAlmostEqualTo(ders[2][0], checkders[2][0]) && MathUtils::IsAlmostEqualTo(ders[2][1], checkders[2][1]) && MathUtils::IsAlmostEqualTo(ders[2][2], checkders[2][2]));

	checkders = { {0.125,-0.5,1.0},{0.75,0.0,-2.0},{0.125,0.5,1.0} };
	for (int i = degree; i >= 0; i--)
	{
		double b = Polynomials::OneBasisFunction(spanIndex - i, degree, knotVector, 5.0 / 2);
		EXPECT_TRUE(MathUtils::IsAlmostEqualTo(check[2 - i],b));

		auto oneders = Polynomials::OneBasisFunctionDerivative(spanIndex - i, degree, 2, knotVector, 5.0 / 2);
		EXPECT_TRUE(MathUtils::IsAlmostEqualTo(oneders[0], checkders[2 - i][0]) && 
					MathUtils::IsAlmostEqualTo(oneders[1], checkders[2 - i][1]) && 
					MathUtils::IsAlmostEqualTo(oneders[2], checkders[2 - i][2]));
	}
}