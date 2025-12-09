/*
 * Author:
 * 2025/11/30 - Yuqing Liang (BIMCoder Liang)
 * bim.frankliang@foxmail.com
 * 
 *
 * Use of this source code is governed by a LGPL-2.1 license that can be found in
 * the LICENSE file.
 */

using System;
using System.Runtime.InteropServices;

namespace LNLibSharp
{
    [StructLayout(LayoutKind.Sequential)]
    public struct LN_ArcInfo
    {
        public double Radius;
        public XYZ Center;
    }

    [StructLayout(LayoutKind.Sequential)]
    public struct LN_Mesh
    {
        public IntPtr Vertices;
        public int VerticesCount;

        public IntPtr Facets;
        public int FacetsCount;

        public IntPtr UVs;
        public int UVsCount;

        public IntPtr UVIndices;
        public int UVIndicesCount;

        public IntPtr Normals;
        public int NormalsCount;

        public IntPtr NormalIndices;
        public int NormalIndicesCount;
    }
}