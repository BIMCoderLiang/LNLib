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
#include "XYZ_CAPI.h"
#include <cstring>

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

Matrix4d_C LNLIB_MATRIX_create()
{
	return FromMatrix4d(LNLib::Matrix4d());
}
Matrix4d_C LNLIB_MATRIX_create_from_xyz(XYZ_C basis_x, XYZ_C basis_y, XYZ_C basis_z, XYZ_C origin)
{
	return FromMatrix4d(LNLib::Matrix4d::Matrix4d(ToXYZ(basis_x), ToXYZ(basis_y), ToXYZ(basis_z), ToXYZ(origin)));
}

Matrix4d_C LNLIB_MATRIX_create_translation(XYZ_C vector)
{
	return FromMatrix4d(LNLib::Matrix4d::CreateTranslation(ToXYZ(vector)));
}

Matrix4d_C LNLIB_MATRIX_create_rotation(XYZ_C axis, double rad)
{
	return FromMatrix4d(LNLib::Matrix4d::CreateRotation(ToXYZ(axis), rad));
}

Matrix4d_C LNLIB_MATRIX_create_rotation_at_point(const XYZ_C& origin, const XYZ_C& axis, double rad)
{
	return FromMatrix4d(LNLib::Matrix4d::CreateRotationAtPoint(ToXYZ(origin), ToXYZ(axis), rad));
}

Matrix4d_C LNLIB_MATRIX_create_scale(XYZ_C scale)
{
	return FromMatrix4d(LNLib::Matrix4d::CreateScale(ToXYZ(scale)));
}

Matrix4d_C LNLIB_MATRIX_create_reflection(XYZ_C normal)
{
	return FromMatrix4d(LNLib::Matrix4d::CreateReflection(ToXYZ(normal)));
}

Matrix4d_C LNLIB_MATRIX_create_reflection_by_distance(XYZ_C normal, double distance_from_origin)
{
	return FromMatrix4d(LNLib::Matrix4d::CreateReflection(ToXYZ(normal), distance_from_origin));
}

XYZ_C LNLIB_MATRIX_get_basis_x(Matrix4d_C m)
{
	return FromXYZ(ToMatrix4d(m).GetBasisX());
}

Matrix4d_C LNLIB_MATRIX_set_basis_x(Matrix4d_C m, XYZ_C basisX)
{
	LNLib::XYZ x = ToXYZ(basisX);
	LNLib::Matrix4d mat = ToMatrix4d(m);
	mat.SetBasisX(x);
	return FromMatrix4d(mat);
}

XYZ_C LNLIB_MATRIX_get_basis_y(Matrix4d_C m)
{
	return FromXYZ(ToMatrix4d(m).GetBasisY());
}

Matrix4d_C LNLIB_MATRIX_set_basis_y(Matrix4d_C m, XYZ_C basisY)
{
	LNLib::XYZ y = ToXYZ(basisY);
	LNLib::Matrix4d mat = ToMatrix4d(m);
	mat.SetBasisY(y);
	return FromMatrix4d(mat);
}

XYZ_C LNLIB_MATRIX_get_basis_z(Matrix4d_C m)
{
	return FromXYZ(ToMatrix4d(m).GetBasisZ());
}

Matrix4d_C LNLIB_MATRIX_set_basis_z(Matrix4d_C m, XYZ_C basisZ)
{
	LNLib::XYZ z = ToXYZ(basisZ);
	LNLib::Matrix4d mat = ToMatrix4d(m);
	mat.SetBasisZ(z);
	return FromMatrix4d(mat);
}

XYZ_C LNLIB_MATRIX_get_basis_w(Matrix4d_C m)
{
	return FromXYZ(ToMatrix4d(m).GetBasisW());
}

Matrix4d_C LNLIB_MATRIX_set_basis_w(Matrix4d_C m, XYZ_C basisW)
{
	LNLib::XYZ w = ToXYZ(basisW);
	LNLib::Matrix4d mat = ToMatrix4d(m);
	mat.SetBasisW(w);
	return FromMatrix4d(mat);
}

XYZ_C LNLIB_MATRIX_of_point(Matrix4d_C m, XYZ_C point)
{
	return FromXYZ(ToMatrix4d(m).OfPoint(ToXYZ(point)));
}

XYZ_C LNLIB_MATRIX_of_vector(Matrix4d_C m, XYZ_C vector)
{
	return FromXYZ(ToMatrix4d(m).OfVector(ToXYZ(vector)));
}

Matrix4d_C LNLIB_MATRIX_multiply(Matrix4d_C a, Matrix4d_C b)
{
	return FromMatrix4d(ToMatrix4d(a).Multiply(ToMatrix4d(b)));
}

int LNLIB_MATRIX_get_inverse(Matrix4d_C m, Matrix4d_C* out_inverse)
{
	LNLib::Matrix4d inv;
	bool success = ToMatrix4d(m).GetInverse(inv);
	if (success && out_inverse) {
		*out_inverse = FromMatrix4d(inv);
	}
	return success ? 1 : 0;
}

Matrix4d_C LNLIB_MATRIX_get_transpose(Matrix4d_C m)
{
	return FromMatrix4d(ToMatrix4d(m).GetTranspose());
}

XYZ_C LNLIB_MATRIX_get_scale(Matrix4d_C m)
{
	return FromXYZ(ToMatrix4d(m).GetScale());
}

double LNLIB_MATRIX_get_determinant(Matrix4d_C m)
{
	return ToMatrix4d(m).GetDeterminant();
}

int LNLIB_MATRIX_is_identity(Matrix4d_C m)
{
	bool success = ToMatrix4d(m).IsIdentity();
	return success ? 1 : 0;
}

int LNLIB_MATRIX_has_reflection(Matrix4d_C m)
{
	bool success = ToMatrix4d(m).HasReflection();
	return success ? 1 : 0;
}

int LNLIB_MATRIX_is_translation(Matrix4d_C m)
{
	bool success = ToMatrix4d(m).IsTranslation();
	return success ? 1 : 0;
}