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
    public static partial class LNLibAPI
    {
        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
        public static extern CurveCurveIntersectionType intersection_compute_rays(
            XYZ point0, XYZ vector0,
            XYZ point1, XYZ vector1,
            out double out_param0,
            out double out_param1,
            out XYZ out_intersect_point);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
        public static extern LinePlaneIntersectionType intersection_compute_line_and_plane(
            XYZ plane_normal,
            XYZ point_on_plane,
            XYZ point_on_line,
            XYZ line_direction,
            out XYZ out_intersect_point);
    }
}