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

	inline constexpr scalar_t scalar_t_max_value() noexcept
	{
		return std::numeric_limits<scalar_t>::max();
	}

	inline constexpr scalar_t scalar_t_min_value() noexcept
	{
		return std::numeric_limits<scalar_t>::min();
	}

	inline scalar_t sqrt(scalar_t v) { return std::sqrt(v); }
	inline scalar_t acos(scalar_t v) { return std::acos(v); }
	inline scalar_t abs(scalar_t v) { return std::fabs(v); }
	inline scalar_t sin(scalar_t v) { return std::sin(v); }
	inline scalar_t cos(scalar_t v) { return std::cos(v); }
	inline scalar_t exp(scalar_t v) { return std::exp(v); }

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

		scalar_t Dot(const Vector& rhs) const noexcept
		{
			return x * rhs.x + y * rhs.y + z * rhs.z;
		}

		scalar_t LengthSqr() const noexcept
		{
			return x * x + y * y + z * z;
		}

		scalar_t Length() const noexcept
		{
			return sqrt(LengthSqr());
		}

		// Normalize
		Vector Unit() const noexcept
		{
			scalar_t l = Length();
			return (l != scalar_t(0)) ? Vector {x / l, y / l, z / l} : Vector {scalar_t(0), scalar_t(0), scalar_t(0)};
		}

		// Returns positive angle between 0 and M_PI
		scalar_t Angle(const Vector& rhs) const noexcept
		{
			return acos(Unit().Dot(rhs));
		}

		Vector Proj(const Vector& rhs) const noexcept
		{
			scalar_t v2ls = rhs.Dot(rhs);
			return (v2ls != scalar_t(0)) ? rhs * (Dot(rhs) / v2ls) : Vector {x, y, z};
		}

		Vector NComp(const Vector& rhs) const noexcept
		{
			return *this - Proj(rhs);
		}

		// Normal distance to line(v1,v2)
		scalar_t NDist(const Vector& v1, const Vector& v2) const noexcept
		{
			return (*this - v1).NComp(v2 - v1).Length();
		}

		static Vector ZERO;

		scalar_t x;
		scalar_t y;
		scalar_t z;
	};

	// ==========================================================================================
	struct BoundingBox
	{
		// Test if another bounding box intersects.
		bool Intersects(const BoundingBox& rhs) const
		{
			if (rhs.max.x < min.x || rhs.min.x > max.x || rhs.max.y < min.y || rhs.min.y > max.y || rhs.max.z < min.z || rhs.min.z > max.z) {
				return false;
			}
			return true;
		}

		// Merge another bounding box.
		void Merge(const BoundingBox& box)
		{
			if (box.min.x < min.x) { min.x = box.min.x; }
			if (box.min.y < min.y) { min.y = box.min.y; }
			if (box.min.z < min.z) { min.z = box.min.z; }
			if (box.max.x > max.x) { max.x = box.max.x; }
			if (box.max.y > max.y) { max.y = box.max.y; }
			if (box.max.z > max.z) { max.z = box.max.z; }
		}

		Vector min;
		Vector max;
	};
}
