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
        public double m00, m01, m02, m03,
                      m10, m11, m12, m13,
                      m20, m21, m22, m23,
                      m30, m31, m32, m33;
    }
    public static partial class LNLibMatrix4d
    {
        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_MATRIX_create")]
        public static extern Matrix4d Create();

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_MATRIX_create_from_xyz")]
        public static extern Matrix4d CreateFromXYZ(
            XYZ basis_x,
            XYZ basis_y,
            XYZ basis_z,
            XYZ origin);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_MATRIX_create_reflection")]
        public static extern Matrix4d CreateReflection(
            XYZ normal);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_MATRIX_create_reflection_by_distance")]
        public static extern Matrix4d CreateReflectionByDistance(
            XYZ normal,
            double distance_from_origin);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_MATRIX_create_rotation")]
        public static extern Matrix4d CreateRotation(
            XYZ axis,
            double rad);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_MATRIX_create_rotation_at_point")]
        public static extern Matrix4d CreateRotationAtPoint(
            XYZ origin,
            XYZ axis,
            double rad);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_MATRIX_create_translation")]
        public static extern Matrix4d CreateTranslation(
            XYZ vector);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_MATRIX_create_scale")]
        public static extern Matrix4d CreateScale(
            XYZ scale);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_MATRIX_get_basis_x")]
        public static extern XYZ GetBasisX(
            Matrix4d m);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_MATRIX_set_basis_x")]
        public static extern Matrix4d SetBasisX(
            Matrix4d m,
            XYZ basis_x);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_MATRIX_get_basis_y")]
        public static extern XYZ GetBasisY(
            Matrix4d m);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_MATRIX_set_basis_y")]
        public static extern Matrix4d SetBasisY(
            Matrix4d m,
            XYZ basis_y);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_MATRIX_get_basis_z")]
        public static extern XYZ GetBasisZ(
            Matrix4d m);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_MATRIX_set_basis_z")]
        public static extern Matrix4d SetBasisZ(
            Matrix4d m,
            XYZ basis_z);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_MATRIX_get_basis_w")]
        public static extern XYZ GetBasisW(
            Matrix4d m);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_MATRIX_set_basis_w")]
        public static extern Matrix4d SetBasisW(
            Matrix4d m,
            XYZ basis_w);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_MATRIX_multiply")]
        public static extern Matrix4d Multiply(
            Matrix4d a,
            Matrix4d b);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_MATRIX_of_point")]
        public static extern XYZ OfPoint(
            Matrix4d m,
            XYZ point);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_MATRIX_of_vector")]
        public static extern XYZ OfVector(
            Matrix4d m,
            XYZ vector);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_MATRIX_get_inverse")]
        public static extern int GetInverse(
            Matrix4d m,
            out Matrix4d out_inverse);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_MATRIX_get_transpose")]
        public static extern Matrix4d GetTranspose(
            Matrix4d m);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_MATRIX_get_scale")]
        public static extern XYZ GetScale(
            Matrix4d m);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_MATRIX_get_determinant")]
        public static extern double GetDeterminant(
            Matrix4d m);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_MATRIX_is_identity")]
        public static extern int IsIdentity(
            Matrix4d m);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_MATRIX_has_reflection")]
        public static extern int HasReflection(
            Matrix4d m);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_MATRIX_is_translation")]
        public static extern int IsTranslation(
            Matrix4d m);
    }
}