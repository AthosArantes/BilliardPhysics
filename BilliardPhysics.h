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
	struct Ball
	{
		// Diameter [m]
		scalar_t diameter;
		// Mass [kg]
		scalar_t mass;
		// Mass [kg*m^2]
		scalar_t massMom;

		Vector position;
		Vector velocity;

		// Rotation speed and axe [rad./s] in table coords
		Vector angularVelocity;
		// Rotation matrix
		Matrix rotation;

		// Set to 0 to disable all physics interaction.
		int enabled;

		// Zero based pocket index if fully inside a pocket.
		int pocketIndex;
	};

	struct PocketHole
	{
		Vector position;
		scalar_t radius;
		scalar_t depth;
	};

	struct ColliderShape
	{
		enum class Type
		{
			Point,
			Line,
			Triangle,
		};

		Type type;
		Vector r1; // pos
		Vector r2; // pos
		Vector r3; // pos
		Vector normal;
	};

	struct ColliderShapeProperty
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

	struct ColliderShapeGroup
	{
		std::vector<ColliderShape> shapes;
		ColliderShapeProperty property;
	};

	// ==========================================================================================
	class Engine
	{
	public:
		Engine();
		virtual ~Engine();

	public:
		// Returns the pocket index if the ball can be pocketed, or -1 if the ball cannot be pocketed.
		int CanPocketBall(const Ball& ball) const;

		bool IsBallInColliderRange(const Ball& ball, const ColliderShape& shape) const;
		scalar_t BallDistance(const Ball& ball, const ColliderShape& shape) const;

		// Returns 1 if collision shape and ball strobe away from each other, 0 else (at time dt)
		bool BallCollided(const Ball& ball, const ColliderShape& shape, scalar_t dt) const;
		// Returns 1 if balls strobe away from each other, 0 else (at time dt)
		bool BallCollided(const Ball& b1, const Ball& b2, scalar_t dt) const;

		Vector PerimeterSpeed(const Ball& ball) const;
		Vector PerimeterSpeedNormal(const Ball& ball, const Vector& normal) const;

		scalar_t CalcCollisionTime(const Ball& b1, const Ball& b2) const;
		scalar_t CalcCollisionTime(const Ball& ball, const ColliderShape& shape) const;

		void MoveBalls(scalar_t dt);
		void ApplyRotationMatrix(Ball& ball, scalar_t dt);

		void BallInteraction(Ball& ball, const ColliderShape& shape, const ColliderShapeProperty& shapeProp);
		void BallInteraction(Ball& b1, Ball& b2);
		void BallInteraction(Ball& ball);

		// this one does not remove fallen balls
		void StepSimulationEuler(scalar_t dt, int depth);
		int StepSimulation(scalar_t dt);

	public:
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
		// 3mm radius der auflageflaeche - not used for rollmom (only rotational-friction around spot)
		scalar_t SpotR;
		scalar_t OmegaMin;
		scalar_t AirResistance;

		// Tolerance used to snap ball on plane zero.
		scalar_t ThreshPosition;

		// Collision properties for the table
		ColliderShapeProperty slateProp;

		std::vector<Ball> balls;
		std::vector<PocketHole> pockets;
		std::vector<ColliderShapeGroup> cushions;
		std::vector<ColliderShapeGroup> scene;
	};
}
