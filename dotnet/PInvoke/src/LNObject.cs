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

namespace LNLibSharp
{
[StructLayout(LayoutKind.Sequential)]
public struct LN_ArcInfo
{
    public double radius;
    public XYZ center;
}

[StructLayout(LayoutKind.Sequential)]
public struct LN_Mesh
{
    public IntPtr vertices;
    public int vertices_count;

    public IntPtr faces;
    public int faces_data_count;

    public IntPtr uvs;
    public int uvs_count;

    public IntPtr uv_indices;
    public int uv_indices_data_count;

    public IntPtr normals;
    public int normals_count;

    public IntPtr normal_indices;
    public int normal_indices_data_count;
}
}