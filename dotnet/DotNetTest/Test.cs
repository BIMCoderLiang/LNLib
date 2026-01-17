/*
 * Author:
 * 2026/01/17 - Yuqing Liang (BIMCoder Liang)
 * bim.frankliang@foxmail.com
 * 
 *
 * Use of this source code is governed by a LGPL-2.1 license that can be found in
 * the LICENSE file.
 */

using LNLibSharp;
using System;
using System.Runtime.InteropServices;

static IntPtr AllocAndPin<T>(T[] array, out GCHandle handle) where T : struct
{
    handle = GCHandle.Alloc(array, GCHandleType.Pinned);
    return handle.AddrOfPinnedObject();
}

var ctrlPts = new XYZW[]
{
            new XYZW { wx = 0, wy = 0, wz = 0, w = 1 },
            new XYZW { wx = 1, wy = 0, wz = 0, w = 1 }
};

var knots = new double[] { 0.0, 0.0, 1.0, 1.0 };

GCHandle ctrlHandle = default;
GCHandle knotHandle = default;

try
{
    IntPtr ctrlPtr = AllocAndPin(ctrlPts, out ctrlHandle);
    IntPtr knotPtr = AllocAndPin(knots, out knotHandle);

    var curve = new LN_NurbsCurve
    {
        degree = 1,
        knot_vector = knotPtr,
        knot_count = knots.Length,
        control_points = ctrlPtr,
        control_point_count = ctrlPts.Length
    };

    int MAX_POINTS = 10;
    XYZ[] tessPoints = new XYZ[MAX_POINTS];

    int pointCount = LNLibNurbsCurve.Tessellate(curve, tessPoints);

    if (pointCount <= 0 || pointCount > MAX_POINTS)
    {
        Console.WriteLine($"Tessellate failed or returned invalid count: {pointCount}");
        return;
    }

    Console.WriteLine($"Tessellated {pointCount} points:");
    for (int i = 0; i < pointCount; i++)
    {
        Console.WriteLine($"  [{i}] = {tessPoints[i]}");
    }
}
finally
{
    if (ctrlHandle.IsAllocated) ctrlHandle.Free();
    if (knotHandle.IsAllocated) knotHandle.Free();
}