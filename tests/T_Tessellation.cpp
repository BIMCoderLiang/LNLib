#include "gtest/gtest.h"
#include "NurbsCurve.h"
#include "NurbsSurface.h"
#include "XYZ.h"
#include "XYZW.h"
#include "MathUtils.h"
#include "Voronoi.h"
#include "LNObject.h"

using namespace LNLib;

TEST(Test_Tessellation, Curve)
{
	XYZ center = XYZ(0, 0, 0);
	XYZ xAxis = XYZ(1, 0, 0);
	XYZ yAxis = XYZ(0, 1, 0);
	double radius = 100;
	LN_NurbsCurve curve;
	bool createArc = NurbsCurve::CreateArc(center, xAxis, yAxis, 0, 2 * Constants::Pi, radius, radius, curve);
	std::vector<double> knots = curve.KnotVector;
	std::vector<XYZ> points = NurbsCurve::Tessellate(curve);
	EXPECT_TRUE(points.size() > 0);
	EXPECT_TRUE(points[0].IsAlmostEqualTo(NurbsCurve::GetPointOnCurve(curve, knots[0])));
	EXPECT_TRUE(points[points.size()-1].IsAlmostEqualTo(NurbsCurve::GetPointOnCurve(curve, knots[knots.size()-1])));
}

TEST(Test_Tessellation, Voronoi)
{
	float xValues[4] = { -22, -17, 4,22 };
	float yValues[4] = { -9, 31,13,-5 };

	long count = 4;

	VoronoiDiagramGenerator vdg;
	vdg.setGenerateVoronoi(true);
	vdg.setGenerateDelaunay(true);
	vdg.generateVoronoi(xValues, yValues, count, -100, 100, -100, 100, 0);
	vdg.resetDelaunayEdgesIterator();

	float x1, y1, x2, y2;
	int edges = 0;
	while (vdg.getNextDelaunay(x1, y1, x2, y2))
	{
		edges++;
	}
	EXPECT_TRUE(edges == 5);

	std::vector<std::vector<int>> triangles =  vdg.getTriangles();
	EXPECT_TRUE(triangles.size() == 2);
}

TEST(Test_Tessellation, Surface)
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

	LNLib::LN_Mesh mesh1 = NurbsSurface::Triangulate(surface, 1000, 1000, true);
	EXPECT_TRUE(mesh1.Faces.size() != 0);

	LNLib::LN_Mesh mesh2 = NurbsSurface::Triangulate(surface, 1000, 1000, false);
	EXPECT_TRUE(mesh2.Faces.size() != 0);
}