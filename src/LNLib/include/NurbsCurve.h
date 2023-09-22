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
	class XYZ;
	class XYZW;
	class Matrix4d;
	class LNLIB_EXPORT NurbsCurve
	{
	public:

		/// <summary>
		/// The NURBS Book 2nd Edition Page124
		/// Algorithm A4.1
		/// Compute point on rational B-spline curve.
		/// </summary>
		static XYZ GetPointOnCurve(int degree, const std::vector<double>& knotVector, double paramT, const std::vector<XYZW>& controlPoints);

		/// <summary>
		/// The NURBS Book 2nd Edition Page127
		/// Algorithm A4.2
		/// Compute C(paramT) derivatives from Cw(paramT) deraivatives.
		/// </summary>
		static std::vector<XYZ> ComputeRationalCurveDerivatives(int degree, int derivative, const std::vector<double>& knotVector, double paramT, const std::vector<XYZW>& controlPoints);

		/// <summary>
		/// The NURBS Book 2nd Edition Page151
		/// Algorithm A5.1
		/// Curve knot insertion.
		/// Note that multiplicity + times <= degree.
		/// </summary>
		static void InsertKnot(int degree, const std::vector<double>& knotVector, const std::vector<XYZW>& controlPoints, double insertKnot, int times, std::vector<double>& insertedKnotVector, std::vector<XYZW>& updatedControlPoints);

		/// <summary>
		/// The NURBS Book 2nd Edition Page155
		/// Algorithm A5.2
		/// Computes point on rational B-spline curve.
		/// </summary>
		static XYZ GetPointOnCurveByCornerCut(int degree, const std::vector<double>& knotVector, double paramT, std::vector<XYZW>& controlPoints);

		/// <summary>
		/// The NURBS Book 2nd Edition Page164
		/// Algorithm A5.4
		/// Refine curve knot vector.
		/// </summary>
		static void RefineKnotVector(int degree, const std::vector<double>& knotVector, const std::vector<XYZW>& controlPoints, std::vector<double>& insertKnotElements, std::vector<double>& insertedKnotVector, std::vector<XYZW>& updatedControlPoints);

		/// <summary>
		/// The NURBS Book 2nd Edition Page173
		/// Algorithm A5.6
		/// Decompose curve into Bezier segements.
		/// decomposedControlPoints[i][j] means ith segement jth control point.
		/// </summary>
		static std::vector<std::vector<XYZW>> DecomposeToBeziers(int degree, const std::vector<double>& knotVector, const std::vector<XYZW>& controlPoints);

		/// <summary>
		/// The NURBS Book 2nd Edition Page185
		/// Algorithm A5.8
		/// Curve knot removal.
		/// </summary>
		static void RemoveKnot(int degree, const std::vector<double>& knotVector, const std::vector<XYZW>& controlPoints, double removeKnot, int times, std::vector<double>& restKnotVector, std::vector<XYZW>& updatedControlPoints);

		/// <summary>
		/// The NURBS Book 2nd Edition Page206
		/// Algorithm A5.9
		/// Degree elevate a curve t times.
		/// </summary>
		static void ElevateDegree(int degree, const std::vector<double>& knotVector, const std::vector<XYZW>& controlPoints, int times, std::vector<double>& updatedKnotVector, std::vector<XYZW>& updatedControlPoints);

		/// <summary>
		/// The NURBS Book 2nd Edition Page223
		/// Algorithm A5.11
		/// Degree reduce a curve from degree to degree - 1.
		/// 
		/// This algorithm has more restricts: it is usually return false since maxError greater than tolerance.
		/// </summary>
		static bool ReduceDegree(int degree, const std::vector<double>& knotVector, const std::vector<XYZW>& controlPoints, std::vector<double>& updatedKnotVector, std::vector<XYZW>& updatedControlPoints);

		/// <summary>
		/// The NURBS Book 2nd Edition Page230
		/// Equally spaced parameter values on each candidate span.
		/// </summary>
		static void EquallyTessellate(int degree, const std::vector<double>& knotVector, const std::vector<XYZW>& controlPoints, std::vector<XYZ>& tessellatedPoints, std::vector<double>& correspondingKnots);

		/// <summary>
		/// The NURBS Book 2nd Edition Page230
		/// Point inversion:finding the corresponding parameter make C(u) = P.
		/// </summary>
		static double GetParamOnCurve(int degree, const std::vector<double>& knotVector, const std::vector<XYZW>& controlPoints, const XYZ& givenPoint);

		/// <summary>
		/// The NURBS Book 2nd Edition Page236
		/// Curve make Transform.
		/// </summary>
		static void CreateTransform(const std::vector<XYZW>& controlPoints, const Matrix4d& matrix, std::vector<XYZW>& transformedControlPoints);

		/// <summary>
		/// The NURBS Book 2nd Edition Page255
		/// Reparameterization using a linear rational function : (alpha * u + beta)/(gamma * u + delta)
		/// </summary>
		static void Reparameterization(int degree, const std::vector<double>& knotVector, const std::vector<XYZW>& controlPoints, double alpha, double beta, double gamma, double delta, std::vector<double>& updatedKnotVector, std::vector<XYZW>& updatedControlPoints);

		/// <summary>
		/// The NURBS Book 2nd Edition Page263
		/// Knot reverse for curve reverse operation.
		/// </summary>
		static void ReverseKnotVector(const std::vector<double>& knotVector, std::vector<double> reversedKnotVector);

		/// <summary>
		/// The NURBS Book 2nd Edition Page263
		/// ControlPoints reverse for curve reverse operation.
		/// </summary>
		static void ReverseControlPoints(const std::vector<XYZW>& controlPoints, std::vector<XYZW> reversedControlPoints);

		/// <summary>
		/// The NURBS Book 2nd Edition Page263
		/// Curve reverse,but not use reparameterization.
		/// </summary>
		static void Reverse(const std::vector<double>& knotVector, const std::vector<XYZW>& controlPoints, std::vector<double>& reversedKnotVector, std::vector<XYZW>& reversedControlPoints);

		/// <summary>
		/// The NURBS Book 2nd Edition Page308
		/// Algorithm A7.1
		/// Create arbitrary NURBS arc.
		/// </summary>
		static bool CreateArc(const XYZ& center, const XYZ& xAxis, const XYZ& yAxis, double xRadius, double yRadius, double startRad, double endRad, int& degree, std::vector<double>& knotVector, std::vector<XYZW>& controlPoints);

		/// <summary>
		/// The NURBS Book 2nd Edition Page314
		/// Algorithm A7.2
		/// Create one Bezier conic arc.
		/// </summary>
		static bool CreateOneConicArc(const XYZ& start, const XYZ& startTangent, const XYZ& end, const XYZ& endTangent, const XYZ& pointOnConic, XYZ& projectPoint, double& projectPointWeight);

		/// <summary>
		/// The NURBS Book 2nd Edition Page317
		/// Split arc.
		/// </summary>
		static void SplitArc(const XYZ& start, const XYZ& projectPoint, double projectPointWeight, const XYZ& end, XYZ& insertPointAtStartSide, XYZ& splitPoint, XYZ& insertPointAtEndSide, double insertWeight);

		/// <summary>
		/// The NURBS Book 2nd Edition Page317
		/// Algorithm A7.3
		/// Construct open conic arc.
		/// </summary>
		static bool CreateOpenConic(const XYZ& start, const XYZ& startTangent, const XYZ& end, const XYZ& endTangent, const XYZ& pointOnConic, int& degree, std::vector<double>& knotVector, std::vector<XYZW>& controlPoints);

		/// <summary>
		/// The NURBS Book 2nd Edition Page369
		/// Algorithm A9.1
		/// Global interpolation through n+1 points.
		/// </summary>
		static void GlobalInterpolation(unsigned int degree, const std::vector<XYZ>& throughPoints, std::vector<double>& knotVector, std::vector<XYZW>& controlPoints);

		/// <summary>
		/// The NURBS Book 2nd Edition Page369 - 374
		/// Global interpolation by through points and tangents. (including Algorithm A9.2)
		/// </summary>
		static void GlobalInterpolationByTangents(unsigned int degree, const std::vector<XYZ>& throughPoints, const std::vector<XYZ>& tangents, const std::vector<int>& tangentsIndices,  std::vector<double>& knotVector, std::vector<XYZW>& controlPoints);

		/// <summary>
		/// The NURBS Book 2nd Edition Page395
		/// Local cubic curve interpolation by through points.
		/// </summary>
		static void LocalCubicCurveInterpolation(const std::vector<XYZ>& throughPoints, std::vector<double>& knotVector, std::vector<XYZW>& controlPoints);

		/// <summary>
		/// The NURBS Book 2nd Edition Page410
		/// Least square curve approximation.
		/// </summary>
		static bool LeastSquaresApproximation(unsigned int degree, const std::vector<XYZ>& throughPoints, int controlPointsCount, std::vector<double>& knotVector, std::vector<XYZW>& controlPoints);

		/// <summary>
		/// The NURBS Book 2nd Edition Page413
		/// Algorithm A9.6
		/// Weighted and contrained least squares approximation.
		/// </summary>
		static bool WeightedAndContrainedLeastSquaresApproximation(unsigned int degree, const std::vector<XYZ>& throughPoints, const std::vector<double>& weights, const std::vector<XYZ>& tangents, const std::vector<int>& tangentIndices, const std::vector<double>& weightedTangents, int controlPointsCount, std::vector<double>& knotVector, std::vector<XYZW>& controlPoints);

		/// <summary>
		/// The NURBS Book 2nd Edition Page428
		/// Algorithm A9.8
		/// Get knot removal error bound (nonrational).
		/// </summary>
		static double ComputerRemoveKnotErrorBound(unsigned int degree, const std::vector<double>& knotVector, const std::vector<XYZW>& controlPoints, int removalIndex);

		/// <summary>
		/// The NURBS Book 2nd Edition Page429
		/// Algorithm A9.9
		/// Remove knots from curve by given bound.
		/// </summary>
		static void RemoveKnotsByGivenBound(unsigned int degree, const std::vector<double>& knotVector, const std::vector<XYZW>& controlPoints, const std::vector<double> params, std::vector<double>& error, double maxError, std::vector<double>& updatedKnotVector, std::vector<XYZW>& updatedControlPoints);

		/// <summary>
		/// The NURBS Book 2nd Edition Page431
		/// Algorithm A9.10
		/// Global curve approximation to within bound maxError.
		/// </summary>
		static void GlobalCurveApproximationByErrorBound(unsigned int degree, const std::vector<XYZ>& throughPoints, double maxError, std::vector<double>& knotVector, std::vector<XYZW>& controlPoints);

		/// <summary>
		/// The NURBS Book 2nd Edition Page440
		/// Algorithm A9.11
		/// Fit to tolerance with conic segment.
		/// </summary>
		static bool LocalRationalQuadraticCurveApproximation(int startPointIndex, int endPointIndex, const std::vector<XYZ>& throughPoints, const XYZ& startTangent, const XYZ& endTangent, double maxError, std::vector<XYZW>& middleControlPoints);

		/// <summary>
		/// The NURBS Book 2nd Edition Page448
		/// Algorithm A9.12
		/// Fit to tolerance with cubic segment.
		/// </summary>
		static bool LocalNonRationalCubicCurveApproximation(int startPointIndex, int endPointIndex, const std::vector<XYZ>& throughPoints, const XYZ& startTangent, const XYZ& endTangent, double maxError, std::vector<XYZW>& middleControlPoints);

		/// <summary>
		/// The NURBS Book 2nd Edition Page577
		/// Algorithm A12.1
		/// Unclamp a clamped curve.
		/// </summary>
		static void ToUnclampCurve(unsigned int degree, std::vector<double>& knotVector, std::vector<XYZW>& controlPoints);
	};
}


