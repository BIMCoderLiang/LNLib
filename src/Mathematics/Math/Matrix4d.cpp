/*
 * Author:
 * 2023/06/08 - Yuqing Liang (BIMCoder Liang)
 * bim.frankliang@foxmail.com
 * 微信公众号：BIMCoder梁老师
 *
 * Use of this source code is governed by a GPL-3.0 license license that can be found in
 * the LICENSE file.
 */

#include "Matrix4d.h"
#include "XYZ.h"
#include "XYZW.h"

using namespace LNLib;

XYZ LNLib::Matrix4d::OfPoint(const XYZ& point)
{
    double x = point.GetX() * m_matrix4d[0][0] + point.GetY() * m_matrix4d[1][0] + point.GetZ() * m_matrix4d[2][0] + m_matrix4d[3][0];
    double y = point.GetX() * m_matrix4d[0][1] + point.GetY() * m_matrix4d[1][1] + point.GetZ() * m_matrix4d[2][1] + m_matrix4d[3][1];
    double z = point.GetX() * m_matrix4d[0][2] + point.GetY() * m_matrix4d[1][2] + point.GetZ() * m_matrix4d[2][2] + m_matrix4d[3][2];
    
    double w = point.GetX() * m_matrix4d[0][3] + point.GetY() * m_matrix4d[1][3] + point.GetZ() * m_matrix4d[2][3] + m_matrix4d[3][3];
    
    return XYZW(x,y,z,w).ToXYZ(true);
}
