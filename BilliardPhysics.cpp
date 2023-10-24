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

#include "BilliardPhysics.h"

namespace BilliardPhysics
{
	static const scalar_t SQRTM1 {100000};

	static scalar_t plane_dist(Vector r, Vector rp, Vector n)
	{
		return (r - rp).Mul(n);
	}

	// ==========================================================================================
	Ball::Ball()
	{
		diameter = scalar_t(57) / scalar_t(1000);
		mass = scalar_t(17) / scalar_t(100);
		massMom = mass * diameter * diameter / scalar_t(5);
		
		position = Vector::ZERO;
		velocity = Vector::ZERO;
		angularVelocity = Vector::ZERO;
		rotation = Matrix::IDENTITY;

		pocketIndex = -1;
		enabled = 1;
	}

	// ==========================================================================================
	Engine::Engine()
	{
		Gravity = scalar_t(981) / scalar_t(100);
		MuRoll = scalar_t(2) / scalar_t(100);
		MuSlide = scalar_t(2) / scalar_t(10);
		MuBall = scalar_t(1) / scalar_t(10);

		SlideThreshSpeed = scalar_t(5) / scalar_t(1000);
		SpotR = scalar_t(12) / scalar_t(1000);
		OmegaMin = SlideThreshSpeed / SpotR;
		AirResistance = scalar_t(1) / scalar_t(1000);

		MinColDist = scalar_t_epsilon();

		slateProp.mu = scalar_t(2) / scalar_t(10);
		slateProp.loss0 = scalar_t(6) / scalar_t(10);
		slateProp.loss_max = scalar_t(80) / scalar_t(100);
		slateProp.loss_wspeed = scalar_t(2);
	}

	Engine::~Engine() = default;

	int Engine::CanPocketBall(const Ball& ball) const
	{
		for (size_t i = 0; i < pockets.size(); ++i) {
			const PocketHole& pocket = pockets[i];
			Vector r = ball.position - pocket.position;
			scalar_t a = pocket.radius * scalar_t(2);
			if (abs(r.x) < a || abs(r.y) < a) {
				r.z = scalar_t(0);
				if (r.Abs() < pocket.radius) {
					return i;
				}
			}
		}
		return -1;
	}

	bool Engine::IsBallInColliderRange(const Ball& ball, const ColliderShape& shape) const
	{
		switch (shape.type) {
			case ColliderShape::Type::Triangle:
			{
				Vector dr1 = shape.r2 - shape.r1;
				Vector dr2 = shape.r3 - shape.r2;
				Vector dr3 = shape.r1 - shape.r3;
				Vector n = dr1.Cross(dr2).Unit();
				return (
					plane_dist(ball.position, shape.r1, n.Cross(dr1).Unit()) >= scalar_t(0) &&
					plane_dist(ball.position, shape.r2, n.Cross(dr2).Unit()) >= scalar_t(0) &&
					plane_dist(ball.position, shape.r3, n.Cross(dr3).Unit()) >= scalar_t(0)
				);
			}
			case ColliderShape::Type::Line:
			{
				Vector r = ball.position - shape.r1;
				Vector dr = shape.r2 - shape.r1;
				scalar_t dra = dr.Abs();
				scalar_t d = r.Mul(dr) / dra;
				return (d >= scalar_t(0) && d < dra);
			}
			case ColliderShape::Type::Point:
			{
				return true;
			}
		}
		return true;
	}

	scalar_t Engine::BallDistance(const Ball& ball, const ColliderShape& shape) const
	{
		if (IsBallInColliderRange(ball, shape)) {
			switch (shape.type) {
				case ColliderShape::Type::Triangle:
				{
					return (ball.position - shape.r1).Mul(shape.normal);
				}
				case ColliderShape::Type::Line:
				{
					Vector r = ball.position - shape.r1;
					Vector dr = shape.r2 - shape.r1;
					return (r - r.Proj(dr)).Abs();
				}
				case ColliderShape::Type::Point:
				{
					return (ball.position - shape.r1).Abs();
				}
			}
		}
		return scalar_t_infinity(); //old return = -1.0E20
	}

	bool Engine::BallCollided(const Ball& ball, const ColliderShape& shape, scalar_t dt) const
	{
		switch (shape.type) {
			case ColliderShape::Type::Triangle:
			{
				return (ball.velocity.Mul(shape.normal) > scalar_t(0));
			}
			case ColliderShape::Type::Line:
			{
				Vector ballpos(ball.position + (ball.velocity * dt));
				return (ball.velocity.Mul((ballpos - shape.r1).NComp(shape.r2 - shape.r1)) > scalar_t(0));
			}
			case ColliderShape::Type::Point:
			{
				Vector ballpos(ball.position + (ball.velocity * dt));
				return ((ballpos - shape.r1).Mul(ball.velocity) > scalar_t(0));
			}
		}
		return true;
	}

	bool Engine::BallCollided(const Ball& b1, const Ball& b2, scalar_t dt) const
	{
		Vector b1pos(b1.position + (b1.velocity * dt));
		Vector b2pos(b2.position + (b2.velocity * dt));
		return (b2pos - b1pos).Mul(b2.velocity - b1.velocity) > scalar_t(0);
	}

	Vector Engine::PerimeterSpeed(const Ball& ball) const
	{
		return ball.angularVelocity.Cross(Vector {scalar_t(0), scalar_t(0), -ball.diameter / scalar_t(2)});
	}

	Vector Engine::PerimeterSpeedNormal(const Ball& ball, const Vector& normal) const
	{
		return ball.angularVelocity.Cross(normal * (-ball.diameter / scalar_t(2)));
	}

	scalar_t Engine::CalcCollisionTime(const Ball& b1, const Ball& b2) const
	{
		Vector dv = b1.velocity - b2.velocity;
		Vector dr = b1.position - b2.position;
		scalar_t vs = dv.x * dv.x + dv.y * dv.y + dv.z * dv.z;
		scalar_t rs = dr.x * dr.x + dr.y * dr.y + dr.z * dr.z;
		scalar_t ds = (b1.diameter + b2.diameter) / 2;
		ds *= ds;

		// Avoid division by zero
		if (vs == scalar_t(0)) {
			return scalar_t_infinity();
		}

		scalar_t p = (dv.x * dr.x + dv.y * dr.y + dv.z * dr.z) / vs;
		scalar_t q = (rs - ds) / vs;
		q = (p * p > q) ? sqrt(p * p - q) : SQRTM1;
		scalar_t t1 = -p + q;
		scalar_t t2 = -p - q;

		return (t1 < t2) ? t1 : t2;
	}

	scalar_t Engine::CalcCollisionTime(const Ball& ball, const ColliderShape& shape) const
	{
		constexpr scalar_t inf = scalar_t_infinity();
		scalar_t rval {};

		switch (shape.type) {
			case ColliderShape::Type::Triangle:
			{
				Vector dr(ball.position - shape.r1);
				scalar_t h = dr.Mul(shape.normal) - (ball.diameter / scalar_t(2));
				scalar_t vn = ball.velocity.Mul(shape.normal);
				if (vn == scalar_t(0)) {
					rval = inf;
					break;
				}
				rval = -h / vn;
				break;
			}
			case ColliderShape::Type::Line:
			{
				// del all comps par to cylinder
				Vector dr(shape.r2 - shape.r1);
				Vector r(ball.position - shape.r1);
				r -= r.Proj(dr);
				Vector v = ball.velocity;
				v -= v.Proj(dr);

				scalar_t vls = v.Abssq();
				if (vls == scalar_t(0)) {
					rval = inf;
					break;
				}

				scalar_t ph = v.Mul(r) / vls;
				scalar_t q = (r.Abssq() - (ball.diameter * ball.diameter / scalar_t(4))) / vls;

				scalar_t t1, t2;
				if (ph * ph > q) {
					t1 = -ph + sqrt(ph * ph - q);
					t2 = -ph - sqrt(ph * ph - q);
				} else {
					t1 = SQRTM1;
					t2 = SQRTM1;
				}
				// solve |r+vt|=d/2
				rval = (t1 < t2) ? t1 : t2;
				break;
			}
			case ColliderShape::Type::Point:
			{
				Vector r(ball.position - shape.r1);
				scalar_t ph = ball.velocity.Mul(r) / ball.velocity.Abssq();
				Vector v = ball.velocity; // ??
				scalar_t q = (r.Abssq() - (ball.diameter * ball.diameter / scalar_t(4))) / v.Abssq();
				scalar_t t1, t2;
				if (ph * ph > q) {
					t1 = -ph + sqrt(ph * ph - q);
					t2 = -ph - sqrt(ph * ph - q);
				} else {
					t1 = SQRTM1;
					t2 = SQRTM1;
				}
				rval = (t1 < t2) ? t1 : t2;
				break;
			}
		}

		if (!IsBallInColliderRange(ball, shape)) {
			rval = inf;
		}

		return rval;
	}

	void Engine::MoveBalls(scalar_t dt)
	{
		Vector dx;

		for (Ball& ball : balls) {
			if (!ball.enabled) {
				continue;
			}

			// translate ball
			dx = ball.velocity * dt;
			ball.position += dx;

			// perform rotation
			ApplyRotationMatrix(ball, dt);
		}
	}

	void Engine::ApplyRotationMatrix(Ball& ball, scalar_t dt)
	{
		scalar_t dphi = ball.angularVelocity.Abs() * dt;
		if (dphi <= scalar_t(0)) {
			return;
		}
		ball.rotation.Rotate(ball.angularVelocity, dphi);
	}

	void Engine::BallInteraction(Ball& ball, const ColliderShape& shape, const ColliderShapeProperty& shapeProp)
	{
		Vector hit_normal;

		switch (shape.type) {
			case ColliderShape::Type::Triangle:
			{
				hit_normal = shape.normal;
				break;
			}
			case ColliderShape::Type::Line:
			{
				Vector dr(shape.r2 - shape.r1);
				hit_normal = ball.position - shape.r1;
				hit_normal = (hit_normal - hit_normal.Proj(dr)).Unit();
				break;
			}
			case ColliderShape::Type::Point:
			{
				hit_normal = (ball.position - shape.r1).Unit();
				break;
			}
		}

		Vector vn = ball.velocity.Proj(hit_normal);
		Vector vp = ball.velocity - vn;

		// normal component
		scalar_t loss = shapeProp.loss0 + (shapeProp.loss_max - shapeProp.loss0) * (scalar_t(1) - exp(-vn.Abs() / shapeProp.loss_wspeed));
		Vector dv(vn * -(scalar_t(1) + sqrt(scalar_t(1) - loss)));
		ball.velocity += dv;

		// parallel component
		Vector ps = PerimeterSpeedNormal(ball, hit_normal);
		dv = (ps + vp).Unit() * (-dv.Abs() * shapeProp.mu);
		Vector dw = (dv * (ball.mass / scalar_t(2) / ball.massMom)).Cross(hit_normal * ball.diameter);
		Vector dw2 = dw + ball.angularVelocity.Proj(dw);

		if (dw2.Mul(ball.angularVelocity) < scalar_t(0)) {
			dw = dw - dw2;
			dv = dv.Unit() * (dw.Abs() * scalar_t(2) * ball.massMom / ball.mass / ball.diameter);
		}

		ball.angularVelocity += dw;
		ball.velocity += dv;

		// optionally disable z velocity if no jump shots desired.
		// maybe some angular momentum loss has to be implemented here
	}

	void Engine::BallInteraction(Ball& b1, Ball& b2)
	{
		const scalar_t half = scalar_t(1) / scalar_t(2);

		Vector dvec(b1.position - b2.position);
		Vector duvec = dvec.Unit();

		// balls in coord system of b1
		// stoss
		Ball b1s = b1;
		Ball b2s = b2;
		b2s.velocity.x -= b1s.velocity.x;
		b2s.velocity.y -= b1s.velocity.y;
		b2s.velocity.z -= b1s.velocity.z;
		b1s.velocity = Vector::ZERO;

		Vector dvn = duvec * duvec.Mul(b2s.velocity);
		Vector dvp = b2s.velocity - dvn;

		b2s.velocity = b2s.velocity - dvn;
		b2s.velocity += dvn * ((b2s.mass - b1s.mass) / (b1s.mass + b2s.mass));
		b1s.velocity = dvn * (scalar_t(2) * b2s.mass / (b1s.mass + b2s.mass)); // (momentum transfer)/m1
		b2.velocity = b1.velocity + b2s.velocity;
		b1.velocity = b1.velocity + b1s.velocity;

		// angular momentum transfer
		Vector dpn = b1s.velocity * b1s.mass; // momentum transfer from ball2 to ball1
		Vector perimeter_speed_b1 = b1.angularVelocity.Cross(duvec * (-b1.diameter / scalar_t(2)));
		Vector perimeter_speed_b2 = b2.angularVelocity.Cross(duvec * (b2.diameter / scalar_t(2)));
		Vector fric_dir = ((perimeter_speed_b2 - perimeter_speed_b1) + dvp).Unit();
		Vector dpp = fric_dir * (-dpn.Abs() * MuBall); // dp parallel of ball2
		Vector dw2 = dpp.Cross(duvec) * (b2.diameter / b2.massMom);
		Vector dw1 = dpp.Cross(duvec) * (b1.diameter / b1.massMom);
		Vector dw2max = (b2.angularVelocity - b1.angularVelocity).Proj(dw2) * half;
		Vector dw1max = (b1.angularVelocity - b2.angularVelocity).Proj(dw2) * half;
		if (dw1.Abs() > dw1max.Abs() || dw2.Abs() > dw2max.Abs()) {
			// Correct momentum transfer to max
			dpp *= dw2max.Abs() / dw2.Abs();
			// Correct amg mom transfer to max
			dw2 = dw2max;
			dw1 = dw1max;
		}
		b1.angularVelocity = b1.angularVelocity - dw1;
		b2.angularVelocity = b2.angularVelocity - dw2;

		// parallel momentum transfer due to friction between balls
		Vector dv1 = dpp * -b1.mass;
		Vector dv2 = dpp * b2.mass;
		dv1.z = scalar_t(0);
		dv2.z = scalar_t(0);
		b1.velocity += dv1;
		b2.velocity += dv2;
	}

	void Engine::BallInteraction(Ball& ball)
	{
		Vector hit_normal {scalar_t(0), scalar_t(0), scalar_t(1)};

		Vector vn = ball.velocity.Proj(hit_normal);
		Vector vp = ball.velocity - vn;

		// normal component
		scalar_t loss = slateProp.loss0 + (slateProp.loss_max - slateProp.loss0) * (scalar_t(1) - exp(-vn.Abs() / slateProp.loss_wspeed));
		Vector dv(vn * -(scalar_t(1) + sqrt(scalar_t(1) - loss)));
		ball.velocity += dv;

		// parallel component
		Vector ps = PerimeterSpeedNormal(ball, hit_normal);
		dv = (ps + vp).Unit() * (-dv.Abs() * slateProp.mu);
		Vector dw = (dv * (ball.mass / scalar_t(2) / ball.massMom)).Cross(hit_normal * ball.diameter);
		Vector dw2 = dw + ball.angularVelocity.Proj(dw);

		if (dw2.Mul(ball.angularVelocity) < scalar_t(0)) {
			dw = dw - dw2;
			dv = dv.Unit() * (dw.Abs() * scalar_t(2) * ball.massMom / ball.mass / ball.diameter);
		}

		ball.angularVelocity += dw;

		// REVISION: Check this condition as it creates an unrealistic table bounce.
		if (ps.Abs() != scalar_t(0)) {
			ball.velocity += dv;
		}
	}

	void Engine::StepSimulationEuler(scalar_t dt, int depth)
	{
		scalar_t dt1 {};
		scalar_t dtmin {};

		Ball* colBall = nullptr;
		Ball* colBall2 = nullptr;
		const ColliderShape* colShape = nullptr;
		const ColliderShapeProperty* colShapeProp = nullptr;

		MoveBalls(dt);

		// lambda for collision detection
		auto stepCollisionDetection = [&](Ball& ball, const ColliderShapeGroup& shapeGroup)
		{
			// IMPROVE: Perhaps an aabb check here ?

			for (const ColliderShape& shape : shapeGroup.shapes) {
				dt1 = CalcCollisionTime(ball, shape);
				if (dt1 < dtmin && dt1 > -dt) {
					if (!BallCollided(ball, shape, -dt)) {
						// dont strobe apart
						dtmin = dt1;
						colBall = &ball;
						colShape = &shape;
						colShapeProp = &shapeGroup.property;
					}
				}
			}
		};

		// checks
		for (Ball& ball : balls) {
			if (!ball.enabled || ball.pocketIndex != -1) {
				continue;
			}

			// Balls in air check additional collisions
			if (ball.position.z != scalar_t(0) || ball.velocity.z != scalar_t(0)) {
				for (auto& shapeGroup : scene) {
					stepCollisionDetection(ball, shapeGroup);
				}
			}

			// Check cushion collisions
			for (const ColliderShapeGroup& shapeGroup : cushions) {
				stepCollisionDetection(ball, shapeGroup);
			}

			// Check ball collisions
			for (Ball& ball2 : balls) {
				if (&ball2 == &ball || !ball2.enabled) {
					continue;
				}

				dt1 = CalcCollisionTime(ball, ball2); // dt1 should be negative
				if (dt1 < dtmin && dt1 > -dt) {
					if (!BallCollided(ball2, ball, -dt)) {
						// dont strobe apart
						dtmin = dt1;
						colBall = &ball2;
						colBall2 = &ball;
						colShape = nullptr;
					}
				}
			}
		}

		if (colShape) {
			MoveBalls(dtmin);
			BallInteraction(*colBall, *colShape, *colShapeProp);

			// TODO: record ball-shape collision

			return StepSimulationEuler(-dtmin, depth + 1);

		} else if (colBall && colBall2) {
			MoveBalls(dtmin);
			BallInteraction(*colBall, *colBall2);

			// TODO: record ball-ball collision

			return StepSimulationEuler(-dtmin, depth + 1);
		}
	}

	int Engine::StepSimulation(scalar_t dt)
	{
		int balls_moving = 0, onPocketIndex;
		Vector accel, waccel, uspeed, uspeed_eff, uspeed2, uspeed_eff2, fricaccel, fricmom, rollmom, totmom;
		Vector bp_pos, bp_normal;
		scalar_t uspeed_eff_par, uspeed_eff2_par, roll_mom_r;
		scalar_t ball_radius, dist;

		// timestep with actual speeds, omegas,...
		StepSimulationEuler(dt, 0);

		// Calc new accelerations and speeds
		for (Ball& ball : balls) {
			if (!ball.enabled) {
				continue;
			}

			// check if balls still moving
			if (ball.velocity.Abs() != scalar_t(0) || ball.angularVelocity.Abs() != scalar_t(0)) {
				balls_moving = 1;
			}

			ball_radius = ball.diameter / scalar_t(2);

			if (ball.pocketIndex == -1) {
				// Ball is not inside pocket

				onPocketIndex = CanPocketBall(ball);
				if (onPocketIndex == -1) {
					// calc accel 3D
					accel = Vector::ZERO; // init acceleration

					// absolute and relative perimeter speed
					uspeed = PerimeterSpeed(ball);
					uspeed_eff = uspeed + ball.velocity;

					// only if ball not flying do sliding/rolling
					if (ball.position.z == scalar_t(0)) {
						// if sliding
						if (uspeed_eff.Abs() > SlideThreshSpeed) {
							// acc caused by friction
							fricaccel = uspeed_eff.Unit() * (-MuSlide * Gravity);
							accel += fricaccel;

							// angular acc caused by friction
							fricmom = fricaccel.Cross(Vector {scalar_t(0), scalar_t(0), -ball_radius}) * ball.mass;
							waccel = fricmom * (scalar_t(-1) / ball.massMom);

							// perform accel
							ball.angularVelocity += waccel * dt;
							ball.velocity += accel * dt;
							uspeed2 = PerimeterSpeed(ball);
							uspeed_eff2 = uspeed2 + ball.velocity;

							// if uspeed_eff passes 0
							uspeed_eff_par = uspeed_eff.Mul(uspeed_eff - uspeed_eff2);
							uspeed_eff2_par = uspeed_eff2.Mul(uspeed_eff - uspeed_eff2);

							if (Vector::ZERO.NDist(uspeed_eff, uspeed_eff2) <= SlideThreshSpeed &&
								((uspeed_eff_par > scalar_t(0) && uspeed_eff2_par < scalar_t(0)) || (uspeed_eff2_par > scalar_t(0) && uspeed_eff_par < scalar_t(0)))
							) {
								// make rolling if uspeed_eff passed 0
								ball.velocity = ball.angularVelocity.Cross(Vector {scalar_t(0), scalar_t(0), ball_radius});
							}

						} else {
							// rolling forces
							fricmom = Vector::ZERO;

							// moment of rotation around ballspot
							if (abs(Vector {scalar_t(0), scalar_t(0), scalar_t(1)}.Mul(ball.angularVelocity)) > OmegaMin) {
								fricmom += Vector {scalar_t(0), scalar_t(0), ball.angularVelocity.z}.Unit() * (MuSlide * ball.mass * Gravity * SpotR);
							}

							// wirkabstand von rollwid.-kraft
							// What does the above comment means ??
							roll_mom_r = MuRoll * ball.massMom / ball.mass / ball.diameter;
							rollmom = Vector {scalar_t(0), scalar_t(0), ball.mass * Gravity * roll_mom_r}.Cross(ball.velocity.Unit());

							totmom = fricmom + rollmom;
							waccel = totmom * (scalar_t(-1) / ball.massMom);

							ball.angularVelocity += waccel * dt;

							// align v with w to assure rolling
							ball.velocity = ball.angularVelocity.Cross(Vector {scalar_t(0), scalar_t(0), ball_radius});
						}

						if (ball.angularVelocity.Abs() < OmegaMin && ball.velocity.Abs() < SlideThreshSpeed) {
							ball.velocity = Vector::ZERO;
							ball.angularVelocity = Vector::ZERO;
						}
					}

					// Snap the ball to the table plane
					dist = (Gravity * dt) * ball.diameter;
					if (abs(ball.position.z) < dist && abs(ball.velocity.z) <= (Gravity * dt)) {
						ball.position.z = scalar_t(0);
						ball.velocity.z = scalar_t(0);
					}

					// Ball collided with the table plane
					if (ball.position.z < scalar_t(0) && ball.velocity.z < scalar_t(0)) {
						// IMPROVE: use better algorithm to find the exact position (not just z) the ball should be placed.
						ball.position.z = scalar_t(0);
						BallInteraction(ball);
					}

				} else {
					// The ball is inside a pocket area
					PocketHole& pocket = pockets[onPocketIndex];

					// Check ball collision with the pocket edge
					if (ball.position.z > -ball_radius) {
						bp_pos = ball.position - pocket.position;
						bp_pos.z = scalar_t(0);

						dist = bp_pos.Abs();
						if (dist > pocket.radius - ball_radius) {
							// Get the possible point of collision on the pocket edge
							bp_pos = pocket.position + bp_pos.Unit() * pocket.radius;
							bp_pos.z = -ball_radius;

							// Check if the ball is intersecting the collision point
							if ((ball.position - bp_pos).Abs() < ball_radius) {
								bp_normal = (ball.position - bp_pos).Unit();
								ball.velocity -= ball.velocity.Proj(bp_normal);
								ball.position = bp_pos + bp_normal * ball_radius;

								// Not affecting the angular velocity though...
							}
						}

					} else {
						// Ball is completely inside the pocket now.
						ball.pocketIndex = onPocketIndex;
					}
				}

			} else {
				// Simpler physics for balls inside pockets
				if (ball.pocketIndex != -1) {
					PocketHole& pocket = pockets[ball.pocketIndex];

					// Bounce the ball inside the pocket
					Vector dr = Vector {ball.position.x, ball.position.y, scalar_t(0)} - Vector {pocket.position.x, pocket.position.y, scalar_t(0)};
					if (dr.Abs() > (pocket.radius - ball_radius)) {
						scalar_t z = ball.velocity.z;
						ball.velocity -= ball.velocity.Proj(ball.position - pocket.position) * scalar_t(2);
						ball.velocity.z = z;

						dr = dr.Unit();
						z = ball.position.z;
						ball.position = pocket.position + (dr * pocket.radius) - (dr * ball_radius);
						ball.position.z = z;
					}

					// Clamp balls inside the pocket
					if (ball.position.z < -pocket.depth) {
						ball.position.z = -pocket.depth;
						ball.velocity = Vector::ZERO;
						ball.angularVelocity = Vector::ZERO;
						ball.enabled = false;

						// TODO: register ball pocket
					}
				}
			}

			// GRAVITY
			{
				// Air resistance
				scalar_t v = ball.velocity.Abs();
				if (v != scalar_t(0)) {
					scalar_t dv = v * AirResistance;
					if (dv > v) {
						ball.velocity = Vector::ZERO;
					} else {
						ball.velocity -= ball.velocity * (dv / v);
					}
				}

				//if ((ball.position.z != 0 && ball.velocity.z != 0) || ball.pocketIndex != -1 || onPocketIndex != -1) {
				if (ball.position.z != scalar_t(0) || ball.pocketIndex != -1 || onPocketIndex != -1) {
					ball.velocity.z += -Gravity * dt;
				}
			}

#if 1
			// Correct intersecting balls
			for (Ball& ball2 : balls) {
				if (&ball2 == &ball || !ball2.enabled) {
					continue;
				}

				dist = ((ball.diameter + ball2.diameter) / scalar_t(2)) - (ball.position - ball2.position).Abs();
				if (dist >= scalar_t(0)) {
					ball.position += (ball.position - ball2.position).Unit() * (dist + MinColDist);
				}
			}
#endif
		}

		return balls_moving;
	}
}
