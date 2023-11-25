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
	LN_Surface surface;
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