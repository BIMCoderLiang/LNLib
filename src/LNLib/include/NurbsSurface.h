/*
 * Author:
 * 2023/06/08 - Yuqing Liang (BIMCoder Liang)
 * bim.frankliang@foxmail.com
 * 微信公众号：BIMCoder梁老师
 *
 * Use of this source code is governed by a GPL-3.0 license that can be found in
 * the LICENSE file.
 */

#pragma once
#include "LNLibDefinitions.h"
#include <vector>

namespace LNLib
{
	class UV;
	class XYZ;
	class XYZW;
	class LNLIB_EXPORT NurbsSurface
	{
	public:

		/// <summary>
		/// The NURBS Book 2nd Edition Page134
		/// Algorithm A4.3
		/// Compute point on rational B-spline surface.
		/// </summary>
		static XYZ GetPointOnSurface(int degreeU, int degreeV, const std::vector<double>& knotVectorU, const std::vector<double>& knotVectorV, UV uv, const std::vector<std::vector<XYZW>>& controlPoints);

		/// <summary>
		/// The NURBS Book 2nd Edition Page137
		/// Algorithm A4.4
		/// Compute S(paramU,paramV) derivatives.
		/// </summary>
		static std::vector<std::vector<XYZ>> ComputeRationalSurfaceDerivatives(int degreeU, int degreeV, int derivative, const std::vector<double>& knotVectorU, const std::vector<double>& knotVectorV, UV uv, const std::vector<std::vector<XYZW>>& controlPoints);

		/// <summary>
		/// The NURBS Book 2nd Edition Page137
		/// Algorithm A5.3
		/// Surface knot insertion along U or V direction.
		/// </summary>
		static void InsertKnot(int degree, const std::vector<double>& knotVector, const std::vector<std::vector<XYZW>>& controlPoints, double insertKnot, int times, bool isUDirection, std::vector<double>& insertedKnotVector, std::vector<std::vector<XYZW>>& updatedControlPoints);

		/// <summary>
		/// The NURBS Book 2nd Edition Page167
		/// Algorithm A5.5
		/// Refine surface knot vector.
		/// </summary>
		static void RefineKnotVector(int degreeU, int degreeV, const std::vector<double>& knotVectorU, const std::vector<double>& knotVectorV, const std::vector<std::vector<XYZW>>& controlPoints, std::vector<double>& insertKnotElements, bool isUDirection, std::vector<double>& insertedKnotVectorU, std::vector<double>& insertedKnotVectorV, std::vector<std::vector<XYZW>>& updatedControlPoints);

		/// <summary>
		/// The NURBS Book 2nd Edition Page177
		/// Algorithm A5.7
		/// Decompose surface into Bezier patches.
		/// decomposedControlPoints[i][j][k] means ith patch jth row kth column control point.
		/// </summary>
		static std::vector<std::vector<std::vector<XYZW>>> DecomposeToBeziers(int degreeU, int degreeV, const std::vector<double>& knotVectorU, const std::vector<double>& knotVectorV, const std::vector<std::vector<XYZW>>& controlPoints);

		/// <summary>
		/// The NURBS Book 2nd Edition Page186
		/// Surface knot removal.
		/// </summary>
		static void RemoveKnot(int degreeU, int degreeV, const std::vector<double>& knotVectorU, const std::vector<double>& knotVectorV, const std::vector<std::vector<XYZW>>& controlPoints, double removeKnot, int times, bool isUDirection, std::vector<double>& restKnotVectorU, std::vector<double>& restKnotVectorV, std::vector<std::vector<XYZW>>& updatedControlPoints);

		/// <summary>
		/// The NURBS Book 2nd Edition Page209
		/// Algorithm A5.10
		/// Degree elevate a surface t times.
		/// </summary>
		static void ElevateDegree(int degreeU, int degreeV, const std::vector<double>& knotVectorU, const std::vector<double>& knotVectorV, const std::vector<std::vector<XYZW>>& controlPoints, int times, bool isUDirection, std::vector<double>& updatedKnotVectorU, std::vector<double>& updatedKnotVectorV, std::vector<std::vector<XYZW>>& updatedControlPoints);

		/// <summary>
		/// The NURBS Book 2nd Edition Page227
		/// Degree reduce U or V Direction Bezier-shape nurbs curve from degree to degree - 1.
		/// </summary>
		static bool ReduceDegree(int degreeU, int degreeV, const std::vector<double>& knotVectorU, const std::vector<double>& knotVectorV, const std::vector<std::vector<XYZW>>& controlPoints, bool isUDirection, std::vector<double>& updatedKnotVectorU, std::vector<double>& updatedKnotVectorV, std::vector<std::vector<XYZW>>& updatedControlPoints);

		/// <summary>
		/// The NURBS Book 2nd Edition Page232
		/// Equally spaced parameter values on each candidate span.
		/// </summary>
		static void EquallyTessellate(int degreeU, int degreeV, const std::vector<double>& knotVectorU, const std::vector<double>& knotVectorV, const std::vector<std::vector<XYZW>>& controlPoints, std::vector<XYZ>& tessellatedPoints, std::vector<UV>& correspondingKnots);

		/// <summary>
		/// The NURBS Book 2nd Edition Page232
		/// Point inversion:finding the corresponding parameter make S(u,v) = P.
		/// </summary>
		static UV GetParamOnSurface(int degreeU, int degreeV, const std::vector<double>& knotVectorU, const std::vector<double>& knotVectorV, const std::vector<std::vector<XYZW>>& controlPoints, const XYZ& givenPoint);

		/// <summary>
		/// The NURBS Book 2nd Edition Page235
		/// Surface Tangent Vector Inversion: finding the corresponding UV tangent [du dv] make T = Su*du+Sv*dv.
		/// </summary>
		static bool GetUVTangent(int degreeU, int degreeV, const std::vector<double>& knotVectorU, const std::vector<double>& knotVectorV, const std::vector<std::vector<XYZW>>& controlPoints, const UV param, const XYZ& tangent, UV& uvTangent);
		
		/// <summary>
		/// The NURBS Book 2nd Edition Page334
		/// construct a NURBS surface by four counter-clockwise points.
		/// point1 ~ point4 are counter-clock placement.
		/// </summary>
		static void CreateBilinearSurface(const XYZ& point1, const XYZ& point2, const XYZ& point3, const XYZ& point4, int& degreeU, int& degreeV, std::vector<double>& knotVectorU, std::vector<double>& knotVectorV, std::vector<std::vector<XYZW>>& controlPoints);

		/// <summary>
		/// The NURBS Book 2nd Edition Page336
		/// Create a right cirular cylinder.
		/// </summary>
		static bool CreateCylindricalSurface(const XYZ& origin, const XYZ& xAxis, const XYZ& yAxis, double startRad, double endRad, double radius, double height, int& degreeU, int& degreeV, std::vector<double>& knotVectorU, std::vector<double>& knotVectorV, std::vector<std::vector<XYZW>>& controlPoints);

		/// <summary>
		/// The NURBS Book 2nd Edition Page337
		/// Create a ruled surface.
		/// </summary>
		static void CreateRuledSurface(int degree0, const std::vector<double>& knotVector0, const std::vector<XYZW> controlPoints0, int degree1, const std::vector<double>& knotVector1, const std::vector<XYZW>& controlPoints1, int& degreeU, int& degreeV, std::vector<double>& knotVectorU, std::vector<double>& knotVectorV, std::vector<std::vector<XYZW>>& controlPoints);

		/// <summary>
		/// The NURBS Book 2nd Edition Page346
		/// Algorithm A8.1
		/// Create a revolved surface.
		/// </summary>
		static bool CreateRevolvedSurface(const XYZ& origin, const XYZ& axis, double rad, const std::vector<XYZW>& generatrixControlPoints, int& degreeU, std::vector<double>& knotVectorU, std::vector<std::vector<XYZW>>& controlPoints);

		/// <summary>
		/// The NURBS Book 2nd Edition Page348
		/// Nonuniform scaling of surface.
		/// </summary>
		static std::vector<std::vector<XYZW>> NonuniformScaling(const std::vector<std::vector<XYZW>>& controlPoints, double xFactor, double yFactor, double zFactor, const XYZ& referencePoint);

		/// <summary>
		/// The NURBS Book 2nd Edition Page359
		/// Create Nurbs corner fillet surface.
		/// </summary>
		static void MakeCornerFilletSurface();

		/// <summary>
		/// The NURBS Book 2nd Edition Page380
		/// Algorithm A9.4
		/// Global surface interpolation.
		/// </summary>
		static void GlobalInterpolation(const std::vector<std::vector<XYZ>>& throughPoints, int degreeU, int degreeV, std::vector<double>& knotVectorU, std::vector<double>& knotVectorV, std::vector<std::vector<XYZW>>& controlPoints);

		/// <summary>
		/// The NURBS Book 2nd Edition Page404
		/// Algorithm A9.5
		/// Local surface interpolation through (n+1)*(m+1) points.
		/// </summary>
		static bool BicubicLocalInterpolation(const std::vector<std::vector<XYZ>>& throughPoints, std::vector<double>& knotVectorU, std::vector<double>& knotVectorV, std::vector<std::vector<XYZW>>& controlPoints);

		/// <summary>
		/// The NURBS Book 2nd Edition Page422
		/// Algorithm A9.7
		/// Global surface approximation with fixed number of control points.
		/// </summary>
		static void GlobalApproximation(const std::vector<std::vector<XYZ>>& throughPoints, int degreeU, int degreeV, int controlPointsRows, int controlPointsColumns, std::vector<double>& knotVectorU, std::vector<double>& knotVectorV, std::vector<std::vector<XYZW>>& controlPoints);

		/// <summary>
		/// The NURBS Book 2nd Edition Page456
		/// </summary>
		static bool CreateSwungSurface();

		/// <summary>
		/// The NURBS Book 2nd Edition Page457
		/// </summary>
		static bool CreateLoftSurface();

		/// <summary>
		/// The NURBS Book 2nd Edition Page472
		/// </summary>
		static bool CreateSweepSurface();

		/// <summary>
		/// The NURBS Book 2nd Edition Page485
		/// </summary>
		static bool CreateGordonSurface();

		/// <summary>
		/// The NURBS Book 2nd Edition Page485
		/// </summary>
		static bool CreateCoonsSurface();
	};
}
