#pragma once

#include "LNLibDefinitions.h"
#include <vector>


namespace LNLib {

	class XYZ;
	class UV;
	class LNLIB_EXPORT BezierSurface
	{

	public:

		/// <summary>
		/// The NURBS Book 2nd Edition Page39
		/// Algorithm A1.7
		/// Compute a point on a Bezier surface by the deCasteljau.
		/// </summary>
		static void GetPointOnSurfaceByDeCasteljau(const std::vector<std::vector<XYZ>>& controlPoints, unsigned int n, unsigned int m, UV uv, XYZ& point);
	};
}



