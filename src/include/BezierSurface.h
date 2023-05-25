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
		/// 
		/// controlpoints with (n+1) * (m+1)
		///  
		///  [0][0]  [0][1] ... ...  [0][m]     ------- v direction
		///  [1][0]  [1][1] ... ...  [1][m]    |
		///    .                               |
		///    .                               u direction
		///    .							   
		///  [n][0]  [n][1] ... ...  [n][m]      
		/// 
		/// </summary>
		static void GetPointOnSurfaceByDeCasteljau(const std::vector<std::vector<XYZ>>& controlPoints, unsigned int n, unsigned int m, UV uv, XYZ& point);
	};
}



