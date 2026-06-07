#include "gtest/gtest.h"
#include "UV.h"
#include "XYZ.h"
#include "XYZW.h"
#include "NurbsSurface.h"
#include "LNObject.h"
using namespace LNLib;

TEST(Test_NurbsSurface, All)
{
	int degreeU = 2; 
	int degreeV = 2;
	std::vector<double> kvU = { 0,0,0,1,2,3,4,4,5,5,5 };
	std::vector<double> kvV = { 0,0,0,1,2,3,3,3 };
	UV uv = UV(5.0 / 2, 1);

	XYZW P20 = XYZW(-1, 2, 4, 1);
	XYZW P21 = XYZW(0, 2, 4, 1);
	XYZW P22 = XYZW(0, 6, 4, 2);
	XYZW P23 = XYZW(0, 2, 0, 1);
	XYZW P24 = XYZW(1, 2, 0, 1);

	XYZW P10 = 0.9 * P20;
	XYZW P11 = 0.9 * P21;
	XYZW P12 = 0.9 * P22;
	XYZW P13 = 0.9 * P23;
	XYZW P14 = 0.9 * P24;

	XYZW P00 = 0.9 * P10;
	XYZW P01 = 0.9 * P11;
	XYZW P02 = 0.9 * P12;
	XYZW P03 = 0.9 * P13;
	XYZW P04 = 0.9 * P14;

	XYZW P30 = XYZW(3, 6, 8, 2);
	XYZW P31 = XYZW(4, 6, 8, 2);
	XYZW P32 = XYZW(12, 24, 12, 6);
	XYZW P33 = XYZW(4, 6, 0, 2);
	XYZW P34 = XYZW(5, 6, 0, 2);

	XYZW P40 = XYZW(3, 2, 4, 1);
	XYZW P41 = XYZW(4, 2, 4, 1);
	XYZW P42 = XYZW(8, 6, 4, 2);
	XYZW P43 = XYZW(4, 2, 0, 1);
	XYZW P44 = XYZW(5, 2, 0, 1);

	XYZW P50 = 1.5 * P40;
	XYZW P51 = 1.5 * P41;
	XYZW P52 = 1.5 * P42;
	XYZW P53 = 1.5 * P43;
	XYZW P54 = 1.5 * P44;

	XYZW P60 = 1.5 * P50;
	XYZW P61 = 1.5 * P51;
	XYZW P62 = 1.5 * P52;
	XYZW P63 = 1.5 * P53;
	XYZW P64 = 1.5 * P54;

	XYZW P70 = 1.5 * P60;
	XYZW P71 = 1.5 * P61;
	XYZW P72 = 1.5 * P62;
	XYZW P73 = 1.5 * P63;
	XYZW P74 = 1.5 * P64;

	std::vector<std::vector<XYZW>> cps = {
	
		{P00, P01, P02, P03, P04},
		{P10, P11, P12, P13, P14},
		{P20, P21, P22, P23, P24},
		{P30, P31, P32, P33, P34},
		{P40, P41, P42, P43, P44},
		{P50, P51, P52, P53, P54},
		{P60, P61, P62, P63, P64},
		{P70, P71, P72, P73, P74},
		
	};
	LN_NurbsSurface surface;
	surface.DegreeU = degreeU;
	surface.DegreeV = degreeV;
	surface.KnotVectorU = kvU;
	surface.KnotVectorV = kvV;
	surface.ControlPoints = cps;
	XYZ result = NurbsSurface::GetPointOnSurface(surface, uv);
	EXPECT_TRUE(result.IsAlmostEqualTo(XYZ(2,98.0/27,68.0/27)));

	std::vector<std::vector<XYZ>> ders =  NurbsSurface::ComputeRationalSurfaceDerivatives(surface,1,uv);
	EXPECT_TRUE(ders[0][0].IsAlmostEqualTo(XYZ(2, 98.0 / 27, 68.0 / 27)));
}

TEST(Test_NurbsSurface, Box)
{
	int degreeU = 3;
	int degreeV = 3;
	std::vector<double> kvU = { 0,0,0,0,0.4,0.6,1,1,1,1 };
	std::vector<double> kvV = { 0,0,0,0,0.4,0.6,1,1,1,1 };
	std::vector<std::vector<XYZW>> controlPoints(6, std::vector<XYZW>(6));

	controlPoints[0][0] = XYZW(0, 0, 0, 1);
	controlPoints[0][1] = XYZW(6.666666, 0, 4, 1);
	controlPoints[0][2] = XYZW(16.666666, 0, 22, 1);
	controlPoints[0][3] = XYZW(33.333333, 0, 22, 1);
	controlPoints[0][4] = XYZW(43.333333, 0, 4, 1);
	controlPoints[0][5] = XYZW(50, 0, 0, 1);

	controlPoints[1][0] = XYZW(0, 6.666667, 0, 1);
	controlPoints[1][1] = XYZW(6.6666667, 6.666667, 9.950068, 1);
	controlPoints[1][2] = XYZW(16.6666666, 6.666667, 9.65541838, 1);
	controlPoints[1][3] = XYZW(33.3333333, 6.666667, 47.21371742, 1);
	controlPoints[1][4] = XYZW(43.3333333, 6.666667, -11.56982167, 1);
	controlPoints[1][5] = XYZW(50, 6.6666667, 0, 1);

	controlPoints[2][0] = XYZW(0, 16.666666, 0, 1);
	controlPoints[2][1] = XYZW(6.6666667, 16.666666, -7.9001371, 1);
	controlPoints[2][2] = XYZW(16.6666666, 16.666666, 46.6891632, 1);
	controlPoints[2][3] = XYZW(33.3333333, 16.666667, -28.4274348, 1);
	controlPoints[2][4] = XYZW(43.3333333, 16.666667, 35.1396433, 1);
	controlPoints[2][5] = XYZW(50, 16.6666667, 0, 1);

	controlPoints[3][0] = XYZW(0, 33.3333333, 0, 1);
	controlPoints[3][1] = XYZW(6.6666667, 33.3333333, 29.2877911, 1);
	controlPoints[3][2] = XYZW(16.6666666, 33.3333333, -30.4644718, 1);
	controlPoints[3][3] = XYZW(33.3333333, 33.3333333, 129.1582990, 1);
	controlPoints[3][4] = XYZW(43.3333333, 33.3333333, -62.1717142, 1);
	controlPoints[3][5] = XYZW(50, 33.333333, 0, 1);

	controlPoints[4][0] = XYZW(0, 43.333333, 0, 1);
	controlPoints[4][1] = XYZW(6.6666667, 43.333333, -10.384636, 1);
	controlPoints[4][2] = XYZW(16.6666666, 43.333333, 59.21371742, 1);
	controlPoints[4][3] = XYZW(33.3333333, 43.333333, -37.7272976, 1);
	controlPoints[4][4] = XYZW(43.3333333, 43.333333, 45.1599451, 1);
	controlPoints[4][5] = XYZW(50, 43.333333, 0, 1);

	controlPoints[5][0] = XYZW(0, 50, 0, 1);
	controlPoints[5][1] = XYZW(6.6666667, 50, 0, 1);
	controlPoints[5][2] = XYZW(16.6666666, 50, 0, 1);
	controlPoints[5][3] = XYZW(33.3333333, 50, 0, 1);
	controlPoints[5][4] = XYZW(43.3333333, 50, 0, 1);
	controlPoints[5][5] = XYZW(50, 50, 0, 1);

	LN_NurbsSurface surface;
	surface.DegreeU = degreeU;
	surface.DegreeV = degreeV;
	surface.KnotVectorU = kvU;
	surface.KnotVectorV = kvV;
	surface.ControlPoints = controlPoints;

	LN_BoundingBox3d box3d =  NurbsSurface::GetBoundingBox(surface);
	EXPECT_NEAR(box3d.MinPoint.X(),0.0, Constants::DoubleEpsilon);
	EXPECT_NEAR(box3d.MaxPoint.X(), 50.0, Constants::DoubleEpsilon);
	EXPECT_NEAR(box3d.MinPoint.Y(), 0.0, Constants::DoubleEpsilon);
	EXPECT_NEAR(box3d.MaxPoint.Y(), 50.0, Constants::DoubleEpsilon);
	EXPECT_NEAR(box3d.MinPoint.Z(), -62.1717, Constants::DoubleEpsilon);
	EXPECT_NEAR(box3d.MaxPoint.Z(), 129.1583, Constants::DoubleEpsilon);

	LN_OrientedBoundingBox3d obox3d = NurbsSurface::GetOrientedBoundingBox(surface);
}