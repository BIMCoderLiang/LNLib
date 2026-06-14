/*
 * Author:
 * 2024/01/13 - Yuqing Liang (BIMCoder Liang)
 * bim.frankliang@foxmail.com
 * 
 *
 * Use of this source code is governed by a LGPL-2.1 license that can be found in
 * the LICENSE file.
 */

#pragma once

namespace LNLib
{
	enum class CurveCurveIntersectionType : int
	{
		Intersecting = 0,
		Parallel = 1,
		Coincident = 2,
		Skew = 3,
	};

	enum class LinePlaneIntersectionType : int
	{
		Intersecting = 0,
		Parallel = 1,
		On = 2,
	};

	enum class CurveNormal :int
	{
		Normal = 0,
		Binormal = 1,
	};

	enum class SurfaceDirection : int
	{
		All = 0,
		UDirection = 1,
		VDirection = 2,
	};

	enum class SurfaceCurvature : int
	{
		Maximum = 0,
		Minimum = 1,
		Gauss = 2,
		Mean = 3,
		Abs = 4,
		Rms = 5
	};

	enum class IntegratorType :int
	{
		Simpson = 0,
		GaussLegendre = 1,
		Chebyshev = 2,
	};

	enum class OffsetType :int
	{
		// Tiller & Hanson Algorithm for C0 profile.
		// 
		// 1. Tiller & Hanson Algorithm should iterative use. 
		// When diff is larger than tolerance, should subdivide curve until less than tolerance.
		// Finally use Merge curve.
		// 
		// 2. Tiller & Hanson Algorithm is not good in negative offset & high degree curve.
		TillerAndHanson = 0,

		// Piegl & Tiller Algorithm for high degree profile,
		// which had better controlled by error tolerance and yet self-intersection had not dealt with.
		PieglAndTiller = 1,
	};

	enum class ExtensionType : int
	{
		// The linear way. 
		// Create the line segment tangent to the source curve at the extension point. 
		// After that, shift by the given length from this point along it (the segment). 
		// This option guarantees the G1 smoothness at the extension point.
		Tangent = 0,

		// In terms of mathematics, it is more correct to call it "the osculating circle method"
		// Create a circle with a radius equal to the curvature radius at the extension point. 
		// The plane, in which this circle is located, is defined by the vectors t − n of the Frenet trihedron
		// (the orthonormal system of vectors: tangent (t), normal (n) and binormal (b)).
		// The plane passes through the point of extension. This method of extension provides at least C2 smoothness.
		Arc = 1,

		// Maintain the equation of the original curve while modifying only the range of the parameter (domain).
		Natural = 2,
	};

	enum class VolumeDirection : int
	{
		U = 0,
		V = 1,
		W = 2,
	};
}



