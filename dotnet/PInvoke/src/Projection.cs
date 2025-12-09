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
    public static partial class LNLibProjection
    {
        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_PROJ_point_to_ray")]
        public static extern XYZ PointToRay(
        XYZ origin,
        XYZ direction,
        XYZ point);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_PROJ_point_to_line")]
        public static extern int PointToLine(
            XYZ start,
            XYZ end,
            XYZ point,
            out XYZ project_point);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_PROJ_stereographic")]
        public static extern XYZ Stereographic(
            XYZ point_on_sphere,
            double radius);
    }
}