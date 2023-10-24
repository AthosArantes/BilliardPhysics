/*
**	Physics code for billiard balls
**
**	Copyright (C) 2001  Florian Berger
**	Email: harpin_floh@yahoo.de, florian.berger@jk.uni-linz.ac.at
**
**	Updated Version foobillard++ started at 12/2010
**	Copyright (C) 2010 - 2013 Holger Schaekel (foobillardplus@go4more.de)
**
**	Refactor and small improvements - oct/2023
**	Copyright (C) 2023 Athos Arantes Pereira
**
**	This program is free software; you can redistribute it and/or modify
**	it under the terms of the GNU General Public License Version 2 as
**	published by the Free Software Foundation;
**
**	This program is distributed in the hope that it will be useful,
**	but WITHOUT ANY WARRANTY; without even the implied warranty of
**	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
**	GNU General Public License for more details.
**
**	You should have received a copy of the GNU General Public License
**	along with this program; if not, write to the Free Software
**	Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
*/

#pragma once

#include <cmath>
#include <limits>

namespace BilliardPhysics
{
	using scalar_t = float;

	inline constexpr scalar_t scalar_t_infinity() noexcept
	{
		return std::numeric_limits<scalar_t>::infinity();
	}

	inline constexpr scalar_t scalar_t_epsilon() noexcept
	{
		return std::numeric_limits<scalar_t>::epsilon();
	}

	inline scalar_t sqrt(scalar_t v) { return std::sqrt(v); }
	inline scalar_t acos(scalar_t v) { return std::acos(v); }
	inline scalar_t abs(scalar_t v) { return std::fabs(v); }
	inline scalar_t sin(scalar_t v) { return std::sin(v); }
	inline scalar_t cos(scalar_t v) { return std::cos(v); }

	// ==========================================================================================
	struct Vector
	{
		Vector operator-(const Vector& rhs) const noexcept
		{
			return {x - rhs.x, y - rhs.y, z - rhs.z};
		}
		Vector& operator-=(const Vector& rhs) noexcept
		{
			x -= rhs.x;
			y -= rhs.y;
			z -= rhs.z;
			return *this;
		}

		Vector operator+(const Vector& rhs) const noexcept
		{
			return {x + rhs.x, y + rhs.y, z + rhs.z};
		}
		Vector& operator+=(const Vector& rhs) noexcept
		{
			x += rhs.x;
			y += rhs.y;
			z += rhs.z;
			return *this;
		}

		Vector operator*(scalar_t rhs) const noexcept
		{
			return {x * rhs, y * rhs, z * rhs};
		}
		Vector& operator*=(scalar_t rhs) noexcept
		{
			x *= rhs;
			y *= rhs;
			z *= rhs;
			return *this;
		}

		Vector operator/(scalar_t rhs) const noexcept
		{
			return {x / rhs, y / rhs, z / rhs};
		}
		Vector& operator/=(scalar_t rhs) noexcept
		{
			x /= rhs;
			y /= rhs;
			z /= rhs;
			return *this;
		}

		Vector Cross(const Vector& rhs) const noexcept
		{
			return {y * rhs.z - rhs.y * z, z * rhs.x - rhs.z * x, x * rhs.y - rhs.x * y};
		}

		scalar_t Mul(const Vector& rhs) const noexcept
		{
			return x * rhs.x + y * rhs.y + z * rhs.z;
		}

		// LengthSqr
		scalar_t Abssq() const noexcept
		{
			return x * x + y * y + z * z;
		}

		// Length
		scalar_t Abs() const noexcept
		{
			scalar_t lensq = Abssq();
			return (lensq >= scalar_t(0)) ? sqrt(Abssq()) : scalar_t(0);
		}

		// Normalize
		Vector Unit() const noexcept
		{
			scalar_t l = Abs();
			return (l != scalar_t(0)) ? Vector {x / l, y / l, z / l} : Vector {scalar_t(0), scalar_t(0), scalar_t(0)};
		}

		// Returns positive angle between 0 and M_PI
		scalar_t Angle(const Vector& rhs) const noexcept
		{
			return acos(Unit().Mul(rhs));
		}

		Vector Proj(const Vector& rhs) const noexcept
		{
			scalar_t v2ls = rhs.Mul(rhs);
			return (v2ls != scalar_t(0)) ? rhs * (Mul(rhs) / v2ls) : Vector {x, y, z};
		}

		Vector NComp(const Vector& rhs) const noexcept
		{
			return *this - Proj(rhs);
		}

		// Normal distance to line(v1,v2)
		scalar_t NDist(const Vector& v1, const Vector& v2) const noexcept
		{
			return (*this - v1).NComp(v2 - v1).Abs();
		}

		static Vector ZERO;

		scalar_t x;
		scalar_t y;
		scalar_t z;
	};

	// ==========================================================================================
	// Used to calculate rotation matrix
	struct Matrix
	{
		void Rotate(const Vector& w, scalar_t phi)
		{
			Vector bx, by, bz; // base
			Vector dp, dp2;
			Vector nax {0, 0, 0};

			if (phi != 0 && w.Abs() > 0) {
				bz = w.Unit();
				if (abs(bz.x) <= abs(bz.y) && abs(bz.x) <= abs(bz.z)) nax = Vector {1, 0, 0};
				if (abs(bz.y) <= abs(bz.z) && abs(bz.y) <= abs(bz.x)) nax = Vector {0, 1, 0};
				if (abs(bz.z) <= abs(bz.x) && abs(bz.z) <= abs(bz.y)) nax = Vector {0, 0, 1};
				bx = (nax - (bz * nax.Mul(bz))).Unit();
				by = bz.Cross(bx);

				scalar_t sinphi = sin(phi);
				scalar_t cosphi = cos(phi);

				for (int i = 0; i < 3; i++) {
					// transform into axis-system
					dp.x = v[i].Mul(bx);
					dp.y = v[i].Mul(by);
					dp.z = v[i].Mul(bz);

					dp2.x = dp.x * cosphi - dp.y * sinphi;
					dp2.y = dp.y * cosphi + dp.x * sinphi;
					dp2.z = dp.z;

					// retransform back
					v[i].x = dp2.x * bx.x + dp2.y * by.x + dp2.z * bz.x;
					v[i].y = dp2.x * bx.y + dp2.y * by.y + dp2.z * bz.y;
					v[i].z = dp2.x * bx.z + dp2.y * by.z + dp2.z * bz.z;
				}
			}
		}

		static Matrix IDENTITY;

		Vector v[3];
	};
}
