/*
 * Author:
 * 2025/11/30 - Yuqing Liang (BIMCoder Liang)
 * bim.frankliang@foxmail.com
 * 
 *
 * Use of this source code is governed by a LGPL-2.1 license that can be found in
 * the LICENSE file.
 */

using System.Runtime.InteropServices;
using static LNLibSharp.LNLibDefinitions;

namespace LNLibSharp
{
    public static partial class LNLibIntersection
    {
        [DllImport(
            LNLIB_CAPI_DLL, 
            EntryPoint = "LNLIB_INTERSECT_compute_rays",
            CallingConvention = CallingConvention.Cdecl)]
        public static extern CurveCurveIntersectionType ComputeRays(
            XYZ point0, XYZ vector0,
            XYZ point1, XYZ vector1,
            out double param0,
            out double param1,
            out XYZ intersectPoint);

        [DllImport(
           LNLIB_CAPI_DLL,
           EntryPoint = "LNLIB_INTERSECT_compute_line_and_plane",
           CallingConvention = CallingConvention.Cdecl)]
        public static extern LinePlaneIntersectionType ComputeLineAndPlane(
            XYZ normal,
            XYZ pointOnPlane,
            XYZ pointOnLine,
            XYZ lineDirection,
            out XYZ intersectPoint);
    }
}