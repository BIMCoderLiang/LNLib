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
    [StructLayout(LayoutKind.Sequential)]
    public struct Matrix4d
    {
        [MarshalAs(UnmanagedType.ByValArray, SizeConst = 16)]
        public double[] m;

        public Matrix4d(int dummy) { m = new double[16]; }
    }

    public static partial class LNLibAPI
    {
        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
        public static extern Matrix4d matrix4d_identity();

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
        public static extern Matrix4d matrix4d_create_translation(XYZ vector);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
        public static extern Matrix4d matrix4d_create_rotation(XYZ axis, double rad);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
        public static extern Matrix4d matrix4d_create_scale(XYZ scale);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
        public static extern Matrix4d matrix4d_create_reflection(XYZ normal);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
        public static extern XYZ matrix4d_get_basis_x(Matrix4d m);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
        public static extern XYZ matrix4d_get_basis_y(Matrix4d m);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
        public static extern XYZ matrix4d_get_basis_z(Matrix4d m);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
        public static extern XYZ matrix4d_get_basis_w(Matrix4d m);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
        public static extern XYZ matrix4d_of_point(Matrix4d m, XYZ point);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
        public static extern XYZ matrix4d_of_vector(Matrix4d m, XYZ vector);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
        public static extern Matrix4d matrix4d_multiply(Matrix4d a, Matrix4d b);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
        public static extern int matrix4d_get_inverse(Matrix4d m, out Matrix4d out_inverse);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
        public static extern double matrix4d_get_determinant(Matrix4d m);
    }
}