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

#include "Math.h"
#include <vector>

namespace BilliardPhysics
{
	class Collider
	{
	public:
		struct Shape
		{
			struct Point
			{
				Vector position;
			};
			struct Line
			{
				Vector position[2];
			};
			struct Triangle
			{
				Vector position[3];
				Vector normal;
			};
			struct Cylinder
			{
				Vector position;
				scalar_t height;
				scalar_t radius;
			};

			enum class Type
			{
				Point,
				Line,
				Triangle,
				Cylinder
			};
			Type type;

			union
			{
				Point point;
				Line line;
				Triangle triangle;
				Cylinder cylinder;
			};

			// Updated by Collider::Update()
			BoundingBox bbox;
		};

		struct Property
		{
			// Friction const
			scalar_t mu;
			// Const loss per hit (0th order in speed)
			scalar_t loss0;
			// Max loss
			scalar_t loss_max;
			// Width of higher order loss curve
			scalar_t loss_wspeed;
		};

	public:
		virtual ~Collider()
		{
		}

		// Update the bounding box of all collider shapes.
		void Update();

	public:
		// Collider shapes.
		std::vector<Shape> shapes;
		// Physics properties of the shapes.
		const Property* property;

		// Bounding box of all collider shapes.
		BoundingBox bbox;
	};

	// ==========================================================================================
	class Pocket
	{
	public:
		virtual ~Pocket()
		{
		}

	public:
		Vector position;
		scalar_t radius;
		scalar_t depth;
	};

	// ==========================================================================================
	class Ball
	{
		friend class Engine;

	public:
		Ball();
		virtual ~Ball()
		{
		}

		void Define(scalar_t radius, scalar_t mass);

		Vector PerimeterSpeed() const noexcept
		{
			return angularVelocity.Cross(Vector {scalar_t(0), scalar_t(0), -radius});
		}
		Vector PerimeterSpeed(const Vector& normal) const noexcept
		{
			return angularVelocity.Cross(normal * -radius);
		}

		bool IsPocketed() const noexcept { return inPocket != nullptr; }

	protected:
		virtual void ApplyRotation(scalar_t dt);
		virtual void OnCollided(const Ball* ball, scalar_t dt);
		virtual void OnCollided(const Collider::Shape* shape, const Collider* collider, scalar_t dt);
		virtual void OnPocketed(const Pocket* pocket);

	protected:
		// Radius [m]
		scalar_t radius;
		// Mass [kg]
		scalar_t mass;
		// Mass [kg*m^2]
		scalar_t massMom;

		scalar_t diameter;
		scalar_t invMass;
		scalar_t invMassMom;

		Vector position;
		Vector velocity;
		// Rotation speed and axe [rad./s] in table coords
		Vector angularVelocity;

		// Bounding box for collider filtering.
		BoundingBox bbox;

		// The pocket if it's fully inside one, nullptr otherwise.
		Pocket* inPocket;
		// Set to false to disable all physics interactions.
		bool enabled;
	};

	// ==========================================================================================
	class Engine
	{
		struct Collision
		{
			size_t hash;

			Ball* ball1;
			Ball* ball2;
			const Collider::Shape* shape;
			const Collider* collider;
		};

	public:
		Engine();
		virtual ~Engine();

	protected:
		// Return the pocket which the ball is within it's area.
		Pocket* GetPocketBallInside(const Ball* ball) const;

		// Returns true if collision shape and ball strobe away from each other, false else (at time dt)
		bool BallCollided(const Ball* ball, const Collider::Shape* shape, scalar_t dt) const;
		// Returns true if balls strobe away from each other, false else (at time dt)
		bool BallCollided(const Ball* b1, const Ball* b2, scalar_t dt) const;

		bool IsBallInColliderRange(const Ball* ball, const Collider::Shape* shape) const;
		scalar_t CalcCollisionTime(const Ball* ball, const Collider::Shape* shape) const;
		scalar_t CalcCollisionTime(const Ball* b1, const Ball* b2) const;

		void MoveBalls(scalar_t dt);
		void CollectColliders(const Ball* ball);

		// Ball interaction with a shape collider
		void BallInteraction(Ball* ball, const Collider::Shape* shape, const Collider* collider);
		// Ball interaction with another ball
		void BallInteraction(Ball* b1, Ball* b2);
		// Ball interaction with table
		void BallInteraction(Ball* ball);

		// Apply collision threshold to intersecting balls.
		void ApplyContactThreshold(Ball* ball);

		// this one does not remove fallen balls
		void StepSimulationEuler(scalar_t dt, int depth);
		int StepSimulation(scalar_t dt);

	private:
		// Collected colliders to test against a ball.
		std::vector<const Collider*> colliders;
		// Collisions that happened exactly at the same time.
		std::vector<Collision> collisions;

	protected:
		// m/s^2
		scalar_t Gravity;
		// Table roll-friction
		scalar_t MuRoll;
		// Table slide-friction
		scalar_t MuSlide;
		// Friction const between ball and ball
		scalar_t MuBall;

		// cm/s
		scalar_t SlideThreshSpeed;
		// Ball spin deceleration rate
		scalar_t SpotR;

		// The distance used to keep balls from touching each other.
		// Very small value.
		scalar_t ContactThreshold;

		std::vector<Ball*> balls;
		std::vector<Pocket*> pockets;

		// Colliders that are always checked.
		// Usually the cushions.
		std::vector<const Collider*> fieldColliders;
		// Colliders that are checked when the ball is in air.
		// Usually table borders, pocket borders, etc.
		std::vector<const Collider*> envColliders;

		// Physics properties for the table plane.
		Collider::Property fieldProperty;
	};
}
