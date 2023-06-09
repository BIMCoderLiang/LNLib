/*
 * Author:
 * 2023/06/22 - Yuqing Liang (BIMCoder Liang)
 * bim.frankliang@foxmail.com
 * 微信公众号：BIMCoder梁老师
 *
 * Use of this source code is governed by a GPL-3.0 license license that can be found in
 * the LICENSE file.
 */

#pragma once
#include "LNLibDefinitions.h"

namespace LNLib
{
	class XYZ;
	class LNLIB_EXPORT Matrix4d
	{
	public:

		Matrix4d();

		Matrix4d(XYZ basisX, XYZ basisY, XYZ basisZ, XYZ origin);

		Matrix4d(double a00, double a01, double a02, double a03,
					double a10, double a11, double a12, double a13,
						double a20, double a21, double a22, double a23,
							double a30, double a31, double a32, double a33);

	public:
		static Matrix4d CreateReflection(const XYZ& origin, const XYZ& normal);
		static Matrix4d CreateRotation(const XYZ& axis, double rad);
		static Matrix4d CreateTranslation(const XYZ& vector);
		static Matrix4d CreateScale(double scale, bool isScaleOrigin);
		static Matrix4d CreateScale(const XYZ& scale);

	public:
		static Matrix4d CreateCamera(const XYZ& eyePoint, const XYZ& lookDirection, const XYZ& upDirection, const XYZ& rightDirection);
		static Matrix4d CreateLookAt(const XYZ& eyePoint, const XYZ& targetPoint, const XYZ& upDirection);
		static Matrix4d Orthogonal(double width, double height, double zNear, double zFar);
		static Matrix4d Perspective(double width, double height, double zNear, double zFar);
		static Matrix4d PerspectiveFov(double fov, double aspect, double zNear, double zFar);
		static Matrix4d PerspectiveMultiFovs(double fovX, double fovY, double zNear, double zFar);

	public:
		void SetBasisX(const XYZ basisX);
		XYZ GetBasisX() const;
		void SetBasisY(const XYZ basisY);
		XYZ GetBasisY() const;
		void SetBasisZ(const XYZ basisZ);
		XYZ GetBasisZ() const;

	public:
		Matrix4d Multiply(const Matrix4d& right);
		XYZ OfPoint(const XYZ& point);
		XYZ OfVector(const XYZ& vector);

	public:
		bool GetInverse(Matrix4d& inverse);
		Matrix4d GetTranspose();
		XYZ GetScale();
		double GetDeterminant();
		bool IsIdentity();
		bool HasReflection();
		bool IsTranslation();

	public:
		Matrix4d& operator =(const Matrix4d & another);

	private:
		double m_matrix4d[4][4];

	};
}

