/*
 * Author:
 * 2026/05/24 - Yuqing Liang (BIMCoder Liang)
 * bim.frankliang@foxmail.com
 * 
 *
 * Use of this source code is governed by a LGPL-2.1 license that can be found in
 * the LICENSE file.
 */

#pragma once
#include "LNLibDefinitions.h"
#include "LNObject.h"
#include "LNEnums.h"
#include <vector>

namespace LNLib
{
	class UVW;
	class XYZ;
	class LNLIB_EXPORT NurbsVolume
	{
	public:

        static bool Check(const LN_NurbsVolume& volume);

        /// <summary>
        /// One member of uvw must be -1, it means not be fixed, namely variation direction.
        /// </summary>
        static LN_NurbsCurve GetIsoCurve(const LN_NurbsVolume& volume, UVW uvw);

        static LNLib::XYZ GetPointOnVolume(const LN_NurbsVolume& volume, UVW uvw);

        static std::vector<std::vector<std::vector<LNLib::XYZ>>> ComputeVolumeRationalDerivatives(const LN_NurbsVolume& volume, int derivative, UVW uvw);

        static void InsertKnot(const LN_NurbsVolume& volume, double insertKnot, int times, VolumeDirection direction, LN_NurbsVolume& result);

        static bool SplitAt(const LN_NurbsVolume& volume, double parameter, VolumeDirection direction, LN_NurbsVolume& first, LN_NurbsVolume& second);

        static bool Extract(const LN_NurbsVolume& volume, double parameter, VolumeDirection direction, LNLib::LN_NurbsSurface surface);

        static void RefineKnotVector(const LN_NurbsVolume& volume, const std::vector<double>& insertKnotElements, VolumeDirection direction, LN_NurbsVolume& result);

        static void ElevateDegree(const LN_NurbsVolume& volume, int times, VolumeDirection direction, LN_NurbsVolume& result);

        static bool IsClosed(const LN_NurbsVolume& volume, VolumeDirection direction);

        static UVW GetParamOnVolume(const LN_NurbsVolume& volume, const LNLib::XYZ& givenPoint);

        static void Reverse(const LN_NurbsVolume& volume, VolumeDirection direction, LN_NurbsVolume& result);

        static void Swap(const LN_NurbsVolume& volume, VolumeDirection direction, VolumeDirection swap, LN_NurbsVolume& result);
	};
}
