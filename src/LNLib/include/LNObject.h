/*
 * Author:
 * 2023/11/15 - Yuqing Liang (BIMCoder Liang)
 * bim.frankliang@foxmail.com
 * 
 *
 * Use of this source code is governed by a LGPL-2.1 license that can be found in
 * the LICENSE file.
 */

#pragma once
#include "LNLibDefinitions.h"
#include "UV.h"
#include "XYZ.h"
#include "XYZW.h"
#include <vector>

namespace LNLib
{
	template <typename T>
	struct LN_BezierCurve
	{
		int Degree = 1;
		std::vector<T> ControlPoints;
	};

	template <typename T>
	struct LN_BezierSurface
	{
		int DegreeU;
		int DegreeV;
		std::vector<std::vector<T>> ControlPoints;
	};

	template <typename T>
	struct LN_BsplineCurve
	{
		int Degree = 1;
		std::vector<double> KnotVector;
		std::vector<T> ControlPoints;
	};

	typedef LN_BsplineCurve<XYZW> LNLIB_EXPORT LN_NurbsCurve;

	template <typename T>
	struct LN_BsplineSurface
	{
		int DegreeU;
		int DegreeV;
		std::vector<double> KnotVectorU;
		std::vector<double> KnotVectorV;
		std::vector<std::vector<T>> ControlPoints;
	};

	typedef LN_BsplineSurface<XYZW> LNLIB_EXPORT LN_NurbsSurface;

	template <typename T>
	struct LN_BsplineVolume
	{
		int DegreeU;
		int DegreeV;
		int DegreeW;

		std::vector<double> KnotVectorU;
		std::vector<double> KnotVectorV;
		std::vector<double> KnotVectorW;

		std::vector<std::vector<std::vector<T>>> ControlPoints;
	};

	typedef LN_BsplineVolume<XYZW> LNLIB_EXPORT LN_NurbsVolume;

	struct LNLIB_EXPORT LN_Mesh
	{
		std::vector<XYZ> Vertices;
		std::vector<std::vector<int>> Faces;
		std::vector<UV> UVs;
		std::vector<int> UVIndices;
		std::vector<XYZ> Normals;
		std::vector<int> NormalIndices;
	};

	struct LNLIB_EXPORT LN_ArcInfo
	{
		double Radius = 0.0;
		XYZ Center = XYZ(0,0,0);
		bool IsCircle = false;
	};

	//AABB Box
	struct LN_BoundingBox3d {
		XYZ MinPoint;
		XYZ MaxPoint;

		bool Intersects(const LN_BoundingBox3d& other) const {
			return (MinPoint.X() <= other.MaxPoint.X() && MaxPoint.X() >= other.MinPoint.X()) &&
					(MinPoint.Y() <= other.MaxPoint.Y() && MaxPoint.Y() >= other.MinPoint.Y()) &&
					(MinPoint.Z() <= other.MaxPoint.Z() && MaxPoint.Z() >= other.MinPoint.Z());
		}
	};

	//OBB Box
	struct LN_OrientedBoundingBox3d {

		XYZ Center;        
		XYZ Axes[3];
		double HalfExtents[3];
	};
}



