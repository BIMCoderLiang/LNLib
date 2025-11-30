/*
 * Author:
 * 2025/11/27 - Yuqing Liang (BIMCoder Liang)
 * bim.frankliang@foxmail.com
 *
 *
 * Use of this source code is governed by a LGPL-2.1 license that can be found in
 * the LICENSE file.
 */

#include "Matrix4d_CAPI.h"
#include "Matrix4d.h"
#include "XYZ.h"
#include <cstring>

extern "C" {

	static LNLib::XYZ ToXYZ(XYZ_C c) {
		return LNLib::XYZ(c.x, c.y, c.z);
	}

	static XYZ_C FromXYZ(const LNLib::XYZ& v) {
		return { v.GetX(), v.GetY(), v.GetZ() };
	}

	static LNLib::Matrix4d ToMatrix4d(Matrix4d_C c) {
		LNLib::Matrix4d m;
		for (int i = 0; i < 4; ++i)
			for (int j = 0; j < 4; ++j)
				m.SetElement(i, j, c.m[i * 4 + j]);
		return m;
	}

	static Matrix4d_C FromMatrix4d(const LNLib::Matrix4d& m) {
		Matrix4d_C c;
		for (int i = 0; i < 4; ++i)
			for (int j = 0; j < 4; ++j)
				c.m[i * 4 + j] = m.GetElement(i, j);
		return c;
	}

	LNLIB_EXPORT Matrix4d_C matrix4d_identity()
	{
		return FromMatrix4d(LNLib::Matrix4d());
	}

	LNLIB_EXPORT Matrix4d_C matrix4d_create_translation(XYZ_C vector)
	{
		return FromMatrix4d(LNLib::Matrix4d::CreateTranslation(ToXYZ(vector)));
	}

	LNLIB_EXPORT Matrix4d_C matrix4d_create_rotation(XYZ_C axis, double rad)
	{
		return FromMatrix4d(LNLib::Matrix4d::CreateRotation(ToXYZ(axis), rad));
	}

	LNLIB_EXPORT Matrix4d_C matrix4d_create_scale(XYZ_C scale)
	{
		return FromMatrix4d(LNLib::Matrix4d::CreateScale(ToXYZ(scale)));
	}

	LNLIB_EXPORT Matrix4d_C matrix4d_create_reflection(XYZ_C normal)
	{
		return FromMatrix4d(LNLib::Matrix4d::CreateReflection(ToXYZ(normal)));
	}

	LNLIB_EXPORT XYZ_C matrix4d_get_basis_x(Matrix4d_C m)
	{
		return FromXYZ(ToMatrix4d(m).GetBasisX());
	}

	LNLIB_EXPORT XYZ_C matrix4d_get_basis_y(Matrix4d_C m)
	{
		return FromXYZ(ToMatrix4d(m).GetBasisY());
	}

	LNLIB_EXPORT XYZ_C matrix4d_get_basis_z(Matrix4d_C m)
	{
		return FromXYZ(ToMatrix4d(m).GetBasisZ());
	}

	LNLIB_EXPORT XYZ_C matrix4d_get_basis_w(Matrix4d_C m)
	{
		return FromXYZ(ToMatrix4d(m).GetBasisW());
	}

	LNLIB_EXPORT XYZ_C matrix4d_of_point(Matrix4d_C m, XYZ_C point)
	{
		return FromXYZ(ToMatrix4d(m).OfPoint(ToXYZ(point)));
	}

	LNLIB_EXPORT XYZ_C matrix4d_of_vector(Matrix4d_C m, XYZ_C vector)
	{
		return FromXYZ(ToMatrix4d(m).OfVector(ToXYZ(vector)));
	}

	LNLIB_EXPORT Matrix4d_C matrix4d_multiply(Matrix4d_C a, Matrix4d_C b)
	{
		return FromMatrix4d(ToMatrix4d(a).Multiply(ToMatrix4d(b)));
	}

	LNLIB_EXPORT int matrix4d_get_inverse(Matrix4d_C m, Matrix4d_C* out_inverse)
	{
		LNLib::Matrix4d inv;
		bool success = ToMatrix4d(m).GetInverse(inv);
		if (success && out_inverse) {
			*out_inverse = FromMatrix4d(inv);
		}
		return success ? 1 : 0;
	}

	LNLIB_EXPORT double matrix4d_get_determinant(Matrix4d_C m)
	{
		return ToMatrix4d(m).GetDeterminant();
	}

}