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
        public static extern XYZ projection_point_to_ray(
            XYZ origin,
            XYZ direction,
            XYZ point);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
        public static extern int projection_point_to_line(
            XYZ start,
            XYZ end,
            XYZ point,
            out XYZ out_project_point);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
        public static extern XYZ projection_stereographic(
            XYZ point_on_sphere,
            double radius);
    }
}