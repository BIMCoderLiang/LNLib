/*
 * Author:
 * 2025/11/27 - Yuqing Liang (BIMCoder Liang)
 * bim.frankliang@foxmail.com
 * 
 *
 * Use of this source code is governed by a LGPL-2.1 license that can be found in
 * the LICENSE file.
 */

#pragma once
#include "LNLibDefinitions.h"
#include "XYZ_CAPI.h"
#include "XYZW_CAPI.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct { double m[16]; } Matrix4d_C;

LNLIB_EXPORT Matrix4d_C LNLIB_MATRIX_create();
LNLIB_EXPORT Matrix4d_C LNLIB_MATRIX_create_from_xyz(XYZ_C basis_x, XYZ_C basis_y, XYZ_C basis_z, XYZ_C origin);

LNLIB_EXPORT Matrix4d_C LNLIB_MATRIX_create_reflection(XYZ_C normal);
LNLIB_EXPORT Matrix4d_C LNLIB_MATRIX_create_reflection_by_distance(XYZ_C normal, double distance_from_origin);
LNLIB_EXPORT Matrix4d_C LNLIB_MATRIX_create_rotation(XYZ_C axis, double rad);
LNLIB_EXPORT Matrix4d_C LNLIB_MATRIX_create_rotation_at_point(const XYZ_C& origin, const XYZ_C& axis, double rad);
LNLIB_EXPORT Matrix4d_C LNLIB_MATRIX_create_translation(XYZ_C vector);
LNLIB_EXPORT Matrix4d_C LNLIB_MATRIX_create_scale(XYZ_C scale);

LNLIB_EXPORT XYZ_C		LNLIB_MATRIX_get_basis_x(Matrix4d_C m);
LNLIB_EXPORT Matrix4d_C LNLIB_MATRIX_set_basis_x(Matrix4d_C m, XYZ_C basis_x);
LNLIB_EXPORT XYZ_C		LNLIB_MATRIX_get_basis_y(Matrix4d_C m);
LNLIB_EXPORT Matrix4d_C LNLIB_MATRIX_set_basis_y(Matrix4d_C m, XYZ_C basis_y);
LNLIB_EXPORT XYZ_C		LNLIB_MATRIX_get_basis_z(Matrix4d_C m);
LNLIB_EXPORT Matrix4d_C LNLIB_MATRIX_set_basis_z(Matrix4d_C m, XYZ_C basis_z);
LNLIB_EXPORT XYZ_C		LNLIB_MATRIX_get_basis_w(Matrix4d_C m);
LNLIB_EXPORT Matrix4d_C LNLIB_MATRIX_set_basis_w(Matrix4d_C m, XYZ_C basis_w);

LNLIB_EXPORT Matrix4d_C LNLIB_MATRIX_multiply(Matrix4d_C a, Matrix4d_C b);
LNLIB_EXPORT XYZ_C		LNLIB_MATRIX_of_point(Matrix4d_C m, XYZ_C point);
LNLIB_EXPORT XYZ_C		LNLIB_MATRIX_of_vector(Matrix4d_C m, XYZ_C vector);

LNLIB_EXPORT int		LNLIB_MATRIX_get_inverse(Matrix4d_C m, Matrix4d_C* out_inverse);
LNLIB_EXPORT Matrix4d_C LNLIB_MATRIX_get_transpose(Matrix4d_C m);
LNLIB_EXPORT XYZ_C		LNLIB_MATRIX_get_scale(Matrix4d_C m);
LNLIB_EXPORT double		LNLIB_MATRIX_get_determinant(Matrix4d_C m);
LNLIB_EXPORT int		LNLIB_MATRIX_is_identity(Matrix4d_C m);
LNLIB_EXPORT int		LNLIB_MATRIX_has_reflection(Matrix4d_C m);
LNLIB_EXPORT int		LNLIB_MATRIX_is_translation(Matrix4d_C m);

#ifdef __cplusplus
}
#endif