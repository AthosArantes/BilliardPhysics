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
#include <algorithm>
#include <cassert>

namespace BilliardPhysics
{
	constexpr scalar_t SQRTM1 = scalar_t(1000000);

	static scalar_t plane_dist(Vector r, Vector rp, Vector n)
	{
		return (r - rp).Dot(n);
	}

	// ==========================================================================================
	Ball::Ball() :
		position(Vector::ZERO()),
		velocity(Vector::ZERO()),
		angularVelocity(Vector::ZERO()),
		inPocket(nullptr),
		enabled(true)
	{
		Define(scalar_t(285) / scalar_t(10000), scalar_t(17) / scalar_t(100));
	}

	void Ball::Define(scalar_t radius_, scalar_t mass_)
	{
		radius = radius_;
		mass = mass_;

		diameter = radius * scalar_t(2);
		massMom = scalar_t(2) / scalar_t(5) * mass * (radius * radius);
		invMass = mass > scalar_t(0) ? scalar_t(1) / mass : scalar_t(0);
		invMassMom = massMom > scalar_t(0) ? scalar_t(1) / massMom : scalar_t(0);
	}

	void Ball::ApplyRotation(scalar_t dt)
	{
	}

	void Ball::OnCollided(const Ball* ball, scalar_t dt)
	{
	}

	void Ball::OnCollided(const Collider::Shape* shape, const Collider* collider, scalar_t dt)
	{
	}

	void Ball::OnPocketed(const Pocket* pocket)
	{
	}

	void Ball::OnLoss()
	{
	}

	// ==========================================================================================
	void Collider::Update()
	{
		bbox = BoundingBox {};

		for (Shape& shape : shapes) {
			switch (shape.type) {
				case Shape::Type::Triangle:
				{
					Vector(&v)[3] = shape.triangle.position;
					for (int i = 0; i < 3; ++i) {
						shape.bbox.Merge(BoundingBox {v[i], v[i]});
					}
					break;
				}
				case Shape::Type::Line:
				{
					Vector(&v)[2] = shape.line.position;
					for (int i = 0; i < 2; ++i) {
						shape.bbox.Merge(BoundingBox {v[i], v[i]});
					}
					break;
				}
				case Shape::Type::Point:
				{
					Vector& v = shape.point.position;
					shape.bbox = BoundingBox {v, v};
					break;
				}
				case Shape::Type::Cylinder:
				{
					Vector& v = shape.cylinder.position;
					scalar_t r = shape.cylinder.radius;
					scalar_t h = shape.cylinder.height;
					shape.bbox = BoundingBox {
						Vector {v.x - r, v.y - r, v.z},
						Vector {v.x + r, v.y + r, v.z + h}
					};
					break;
				}
			}

			bbox.Merge(shape.bbox);
		}
	}

	// ==========================================================================================
	Engine::Engine()
	{
		Gravity = scalar_t(981) / scalar_t(100);
		MuRoll = scalar_t(8) / scalar_t(1000);
		MuSlide = scalar_t(2) / scalar_t(10);
		MuBall = scalar_t(5) / scalar_t(100);

		SlideThreshSpeed = scalar_t(5) / scalar_t(1000);
		SlideThreshSpeed *= SlideThreshSpeed;
		SpotR = scalar_t(10);

		ContactThreshold = scalar_t(1) / SQRTM1;
		SlateBound = Vector {scalar_t(2), scalar_t(1), scalar_t(-1)};

		fieldProperty.mu = scalar_t(20) / scalar_t(100);
		fieldProperty.loss0 = scalar_t(60) / scalar_t(100);
		fieldProperty.loss_max = scalar_t(90) / scalar_t(100);
		fieldProperty.loss_wspeed = scalar_t(15) / scalar_t(10);
	}

	Engine::~Engine() = default;

	Pocket* Engine::GetPocketBallInside(const Ball* ball) const
	{
		for (Pocket* pocket : pockets) {
			Vector dr = ball->position - pocket->position;
			dr.z = scalar_t(0);
			if (dr.LengthSqr() < pocket->radius * pocket->radius) {
				return pocket;
			}
		}
		return nullptr;
	}

	bool Engine::BallCollided(const Ball* ball, const Collider::Shape* shape, scalar_t dt) const
	{
		switch (shape->type) {
			case Collider::Shape::Type::Triangle:
			{
				const Collider::Shape::Triangle& triangle = shape->triangle;
				return (ball->velocity.Dot(triangle.normal) > scalar_t(0));
			}
			case Collider::Shape::Type::Line:
			{
				const Collider::Shape::Line& line = shape->line;
				const Vector(&v)[2] = line.position;

				Vector ballpos = ball->position + (ball->velocity * dt);
				return (ball->velocity.Dot((ballpos - v[0]).NComp(v[1] - v[0])) > scalar_t(0));
			}
			case Collider::Shape::Type::Point:
			{
				const Collider::Shape::Point& point = shape->point;

				Vector ballpos = ball->position + (ball->velocity * dt);
				return ((ballpos - point.position).Dot(ball->velocity) > scalar_t(0));
			}
			case Collider::Shape::Type::Cylinder:
			{
				const Collider::Shape::Cylinder& cylinder = shape->cylinder;

				Vector ballpos = ball->position + (ball->velocity * dt);
				Vector n = ballpos - cylinder.position;
				n.z = scalar_t(0);
				Vector r = n.Unit() * cylinder.radius;
				r.z = ballpos.z;

				return ((ballpos - r).Dot(ball->velocity) > scalar_t(0));
			}
		}
		return true;
	}

	bool Engine::BallCollided(const Ball* b1, const Ball* b2, scalar_t dt) const
	{
		Vector b1pos = b1->position + (b1->velocity * dt);
		Vector b2pos = b2->position + (b2->velocity * dt);
		return (b2pos - b1pos).Dot(b2->velocity - b1->velocity) > scalar_t(0);
	}

	bool Engine::IsBallInColliderRange(const Ball* ball, const Collider::Shape* shape) const
	{
		switch (shape->type) {
			case Collider::Shape::Type::Triangle:
			{
				const Collider::Shape::Triangle& triangle = shape->triangle;
				const Vector(&v)[3] = triangle.position;

				Vector dr1 = v[1] - v[0];
				Vector dr2 = v[2] - v[1];
				Vector dr3 = v[0] - v[2];
				Vector n = dr1.Cross(dr2).Unit();
				return (
					plane_dist(ball->position, v[0], n.Cross(dr1).Unit()) >= scalar_t(0) &&
					plane_dist(ball->position, v[1], n.Cross(dr2).Unit()) >= scalar_t(0) &&
					plane_dist(ball->position, v[2], n.Cross(dr3).Unit()) >= scalar_t(0)
				);
			}
			case Collider::Shape::Type::Line:
			{
				const Collider::Shape::Line& line = shape->line;
				const Vector(&v)[2] = line.position;

				Vector r = ball->position - v[0];
				Vector dr = v[1] - v[0];
				scalar_t dra = dr.Length();
				scalar_t d = r.Dot(dr) / dra;
				return (d >= scalar_t(0) && d < dra);
			}
			case Collider::Shape::Type::Cylinder:
			{
				const Collider::Shape::Cylinder& cylinder = shape->cylinder;
				return (ball->position.z >= cylinder.position.z && ball->position.z <= cylinder.position.z + cylinder.height);
			}
		}
		return true;
	}

	scalar_t Engine::CalcCollisionTime(const Ball* ball, const Collider::Shape* shape) const
	{
		if (!IsBallInColliderRange(ball, shape)) {
			return SQRTM1;
		}

		switch (shape->type) {
			case Collider::Shape::Type::Triangle:
			{
				Vector dr = ball->position - shape->triangle.position[0];
				scalar_t h = dr.Dot(shape->triangle.normal) - ball->radius;
				scalar_t vn = ball->velocity.Dot(shape->triangle.normal);

				return (vn == scalar_t(0)) ? SQRTM1 : -h / vn;
			}
			case Collider::Shape::Type::Line:
			{
				// del all comps par to cylinder
				Vector dr = shape->line.position[1] - shape->line.position[0];
				Vector r = ball->position - shape->line.position[0];
				r -= r.Proj(dr);
				Vector v = ball->velocity;
				v -= v.Proj(dr);

				scalar_t vls = v.LengthSqr();
				scalar_t ph = v.Dot(r) / vls;
				scalar_t q = (r.LengthSqr() - (ball->radius * ball->radius)) / vls;

				if (ph * ph > q) {
					scalar_t t1 = -ph + sqrt(ph * ph - q);
					scalar_t t2 = -ph - sqrt(ph * ph - q);

					// solve |r+vt|=d/2
					return (t1 < t2) ? t1 : t2;
				}

				break;
			}
			case Collider::Shape::Type::Point:
			{
				scalar_t vls = ball->velocity.LengthSqr();
				Vector r = ball->position - shape->point.position;
				scalar_t ph = ball->velocity.Dot(r) / vls;
				scalar_t q = (r.LengthSqr() - (ball->radius * ball->radius)) / vls;

				if (ph * ph > q) {
					scalar_t t1 = -ph + sqrt(ph * ph - q);
					scalar_t t2 = -ph - sqrt(ph * ph - q);
					return (t1 < t2) ? t1 : t2;
				}

				break;
			}
			case Collider::Shape::Type::Cylinder:
			{
				scalar_t vls = ball->velocity.LengthSqr();

				Vector dr = ball->position - shape->cylinder.position;
				dr.z = scalar_t(0);

				scalar_t rs = ball->radius + shape->cylinder.radius;
				rs *= rs;

				scalar_t ph = (ball->velocity.Dot(dr)) / vls;
				scalar_t q = (dr.LengthSqr() - rs) / vls;

				scalar_t t1, t2;
				if (ph * ph > q) {
					t1 = -ph + sqrt(ph * ph - q);
					t2 = -ph - sqrt(ph * ph - q);
					return (t1 < t2) ? t1 : t2;
				}

				break;
			}
		}

		return SQRTM1;
	}

	scalar_t Engine::CalcCollisionTime(const Ball* b1, const Ball* b2) const
	{
		Vector dv = b1->velocity - b2->velocity;
		scalar_t vs = dv.LengthSqr();
		if (vs == scalar_t(0)) {
			return SQRTM1;
		}

		Vector dr = b1->position - b2->position;
		scalar_t rs = dr.LengthSqr();
		scalar_t ds = b1->radius + b2->radius;
		ds *= ds;

		scalar_t p = dv.Dot(dr) / vs;
		scalar_t q = (rs - ds) / vs;
		q = (p * p > q) ? sqrt(p * p - q) : SQRTM1;
		scalar_t t1 = -p + q;
		scalar_t t2 = -p - q;

		return (t1 < t2) ? t1 : t2;
	}

	void Engine::MoveBalls(scalar_t dt)
	{
		for (Ball* ball : balls) {
			if (!ball->enabled) {
				continue;
			}

			// Pre-translate bbox
			ball->bbox = BoundingBox {ball->position, ball->position};

			// Translate ball
			Vector dx = ball->velocity * dt;
			ball->position += dx;

			// Post-translate bbox
			ball->bbox.Merge(BoundingBox {ball->position, ball->position});
			ball->bbox.Dilate(ball->radius);
		}
	}

	void Engine::CollectColliders(const Ball* ball)
	{
		colliders.clear();

		// Check cushion colliders
		for (const Collider* collider : fieldColliders) {
			if (collider->bbox.Intersects(ball->bbox)) {
				colliders.push_back(collider);
			}
		}

		bool onGround = ((ball->position.z - ball->radius) == scalar_t(0));

		// Balls in air check additional colliders
		if (!onGround || ball->velocity.z != scalar_t(0)) {
			for (const Collider* collider : envColliders) {
				if (collider->bbox.Intersects(ball->bbox)) {
					colliders.push_back(collider);
				}
			}
		}
	}

	void Engine::BallInteraction(Ball* ball, const Collider::Shape* shape, const Collider* collider)
	{
		Vector hit_normal;

		if (shape) {
			switch (shape->type) {
				case Collider::Shape::Type::Triangle:
				{
					const Collider::Shape::Triangle& triangle = shape->triangle;
					hit_normal = triangle.normal;
					break;
				}
				case Collider::Shape::Type::Line:
				{
					const Collider::Shape::Line& line = shape->line;
					Vector dr = line.position[1] - line.position[0];
					hit_normal = ball->position - line.position[0];
					hit_normal = (hit_normal - hit_normal.Proj(dr)).Unit();
					break;
				}
				case Collider::Shape::Type::Point:
				{
					const Collider::Shape::Point& point = shape->point;
					hit_normal = (ball->position - point.position).Unit();
					break;
				}
				case Collider::Shape::Type::Cylinder:
				{
					const Collider::Shape::Cylinder& cylinder = shape->cylinder;
					Vector dr = ball->position - cylinder.position;
					dr.z = scalar_t(0);
					hit_normal = dr.Unit();
					break;
				}
			}

		} else {
			// If no shape then the collision is with the table slate
			hit_normal = Vector {scalar_t(0), scalar_t(0), scalar_t(1)};
		}

		const Collider::Property* prop = (collider && collider->property) ? collider->property : &fieldProperty;
		Vector vn = ball->velocity.Proj(hit_normal);
		Vector vp = ball->velocity - vn;

		// normal component
		scalar_t loss = prop->loss0 + (prop->loss_max - prop->loss0) * (scalar_t(1) - exp(-vn.Length() / prop->loss_wspeed));
		Vector dv = vn * -(scalar_t(1) + sqrt(scalar_t(1) - loss));
		ball->velocity += dv;

		// parallel component
		Vector ps = ball->PerimeterSpeed(hit_normal);
		dv = (ps + vp).Unit() * (-dv.Length() * prop->mu);
		Vector dw = (dv * (ball->mass / scalar_t(2) / ball->massMom)).Cross(hit_normal * ball->diameter);
		Vector dw2 = dw + ball->angularVelocity.Proj(dw);

		if (dw2.Dot(ball->angularVelocity) < scalar_t(0)) {
			dw = dw - dw2;
			dv = dv.Unit() * (dw.Length() * scalar_t(2) * ball->massMom / ball->mass / ball->diameter);
		}

		ball->angularVelocity += dw;
		ball->velocity += dv;

		// optionally disable z velocity if no jump shots desired.
		// maybe some angular momentum loss has to be implemented here
	}

	void Engine::BallInteraction(Ball* b1, Ball* b2)
	{
		const scalar_t half = scalar_t(1) / scalar_t(2);

		Vector dvec = b1->position - b2->position;
		Vector duvec = dvec.Unit();

		// balls in coord system of b1
		// stoss
		Vector b2v = b2->velocity - b1->velocity;
		Vector b1v = Vector::ZERO();

		Vector dvn = duvec * duvec.Dot(b2v);
		Vector dvp = b2v - dvn;

		scalar_t invSumMass = ((b1->mass + b2->mass) > scalar_t(0)) ? scalar_t(1) / (b1->mass + b2->mass) : scalar_t(0);

		b2v = b2v - dvn;
		b2v += dvn * (b2->mass - b1->mass) * invSumMass;
		b1v = dvn * scalar_t(2) * b2->mass * invSumMass; // (momentum transfer)/m1
		b2->velocity = b1->velocity + b2v;
		b1->velocity = b1->velocity + b1v;

		// angular momentum transfer
		Vector dpn = b1v * b1->mass; // momentum transfer from ball2 to ball1
		Vector perimeter_speed_b1 = b1->angularVelocity.Cross(duvec * -b1->radius);
		Vector perimeter_speed_b2 = b2->angularVelocity.Cross(duvec * b2->radius);
		Vector fric_dir = ((perimeter_speed_b2 - perimeter_speed_b1) + dvp).Unit();
		Vector dpp = fric_dir * (-dpn.Length() * MuBall); // dp parallel of ball2

		Vector dw2 = dpp.Cross(duvec) * b2->diameter * b2->invMassMom;
		Vector dw1 = dpp.Cross(duvec) * b1->diameter * b1->invMassMom;
		Vector dw2max = (b2->angularVelocity - b1->angularVelocity).Proj(dw2) * half;
		Vector dw1max = (b1->angularVelocity - b2->angularVelocity).Proj(dw2) * half;
		if (dw1.Length() > dw1max.Length() || dw2.Length() > dw2max.Length()) {
			// Correct momentum transfer to max
			dpp *= dw2max.Length() / dw2.Length();
			// Correct amg mom transfer to max
			dw2 = dw2max;
			dw1 = dw1max;
		}
		b1->angularVelocity = b1->angularVelocity - dw1;
		b2->angularVelocity = b2->angularVelocity - dw2;

		// parallel momentum transfer due to friction between balls
		Vector dv1 = dpp * -b1->mass;
		Vector dv2 = dpp * b2->mass;
		dv1.z = scalar_t(0);
		dv2.z = scalar_t(0);
		b1->velocity += dv1;
		b2->velocity += dv2;
	}

	void Engine::ApplyContactThreshold(Ball* ball)
	{
		for (const Ball* ball2 : balls) {
			if (ball2 == ball || !ball2->enabled) {
				continue;
			}

			Vector dr = ball->position - ball2->position;
			scalar_t dist = ball->radius + ball2->radius;
			if (dr.LengthSqr() <= dist * dist) {
				ball->position = ball2->position + dr.Unit() * (dist + ContactThreshold);
			}
		}
	}

	void Engine::StepSimulationEuler(scalar_t dt, int depth)
	{
		scalar_t dtmin = scalar_t(0);

		MoveBalls(dt);

		// checks
		for (Ball* ball : balls) {
			if (!ball->enabled || ball->IsPocketed()) {
				continue;
			}

			CollectColliders(ball);

			// Check slate collision
			if (ball->position.z < ball->radius && GetPocketBallInside(ball) == nullptr && IsBallInsidePlayfield(ball)) {
				scalar_t ct = (ball->position.z - ball->radius) / ball->velocity.z;
				ct = (scalar_t(1) - ct) * -dt;
				assert(ct >= -dt);

				// Clear collisions vector if this collision happened earlier
				if (ct < dtmin) {
					dtmin = ct;
					collisions.clear();
				}

				// Store this collision
				auto it = std::find_if(collisions.cbegin(), collisions.cend(), [ball](const Collision& c)
				{
					return c.ball1 == ball && c.ball2 == nullptr && c.shape == nullptr;
				});
				if (it == collisions.cend()) {
					collisions.emplace_back(Collision {ball, nullptr, nullptr, nullptr});
				}
			}

			// Check shape collisions
			for (const Collider* collider : colliders) {
				for (const Collider::Shape& shape : collider->shapes) {
					scalar_t ct = CalcCollisionTime(ball, &shape);
					if (ct >= -dt && ct <= dtmin && !BallCollided(ball, &shape, -dt)) {
						// dont strobe apart

						// Clear collisions vector if this collision happened earlier
						if (ct < dtmin) {
							dtmin = ct;
							collisions.clear();
						}

						// Store this collision
						auto it = std::find_if(collisions.cbegin(), collisions.cend(), [ball, shape](const Collision& c)
						{
							return c.ball1 == ball && c.shape == &shape;
						});
						if (it == collisions.cend()) {
							collisions.emplace_back(Collision {ball, nullptr, &shape, collider});
						}
					}
				}
			}

			// Check ball collisions
			for (Ball* ball2 : balls) {
				if (ball2 == ball || !ball2->enabled || ball2->IsPocketed() || !ball->bbox.Intersects(ball2->bbox)) {
					continue;
				}

				scalar_t ct = CalcCollisionTime(ball, ball2);
				if (ct >= -dt && ct <= dtmin && !BallCollided(ball2, ball, -dt)) {
					// dont strobe apart

					// Clear collisions vector if this collision happened earlier
					if (ct < dtmin) {
						dtmin = ct;
						collisions.clear();
					}

					// Store this collision
					auto it = std::find_if(collisions.cbegin(), collisions.cend(), [ball](const Collision& c)
					{
						return c.ball1 == ball || c.ball2 == ball;
					});
					if (it == collisions.cend()) {
						collisions.emplace_back(Collision {ball, ball2, nullptr, nullptr});
					}
				}
			}
		}

		// No collision occurred
		if (collisions.empty()) {
			return;
		}

		MoveBalls(dtmin);

		for (Collision& c : collisions) {
			if (c.shape) {
				BallInteraction(c.ball1, c.shape, c.collider);
				c.ball1->OnCollided(c.shape, c.collider, -dtmin);

			} else if (c.ball2) {
				BallInteraction(c.ball1, c.ball2);
				c.ball1->OnCollided(c.ball2, -dtmin);

			} else {
				// Assure the ball is exactly "on table"
				c.ball1->position.z = c.ball1->radius;

				BallInteraction(c.ball1, nullptr, nullptr);

				// Prevent ball from bouncing forever...
				if (c.ball1->velocity.z < Gravity * dt * scalar_t(2)) {
					c.ball1->velocity.z = scalar_t(0);
				}

				c.ball1->OnCollided(nullptr, nullptr, -dtmin);
			}
		}
		collisions.clear();

		StepSimulationEuler(-dtmin, ++depth);
	}

	int Engine::StepSimulation(scalar_t dt)
	{
		int balls_moving = 0;

		// Timestep with actual speeds, omegas, interactions, etc
		StepSimulationEuler(dt, 0);

		// Calc new accelerations and speeds
		for (Ball* ball : balls) {
			if (!ball->enabled) {
				continue;
			}

			// Check out-of-play balls
			if (ball->position.z < SlateBound.z) {
				ball->enabled = false;
				ball->OnLoss();
				continue;
			}

			// Check if ball still moving
			bool moving = (ball->velocity + ball->angularVelocity != Vector::ZERO());
			bool onGround = ((ball->position.z - ball->radius) == scalar_t(0));
			Pocket* pocket = ball->inPocket;

			if (moving) {
				++balls_moving;

				ball->ApplyRotation(dt);

				if (!pocket) {
					// Ball is not fully inside a pocket, but may be falling into one.
					pocket = GetPocketBallInside(ball);

					if (!pocket) {
						// absolute and relative perimeter speed
						Vector uspeed = ball->PerimeterSpeed();
						Vector uspeed_eff = uspeed + ball->velocity;

						// only if ball not flying do sliding/rolling
						if (onGround) {
							// if sliding
							if (uspeed_eff.LengthSqr() > SlideThreshSpeed) {
								// acc caused by friction
								Vector fricaccel = uspeed_eff.Unit() * -(MuSlide * Gravity);

								// angular acc caused by friction
								Vector fricmom = fricaccel.Cross(Vector {scalar_t(0), scalar_t(0), -ball->radius}) * ball->mass;
								Vector waccel = fricmom * -ball->invMassMom;

								// perform accel
								ball->angularVelocity += waccel * dt;
								ball->velocity += fricaccel * dt;
								Vector uspeed2 = ball->PerimeterSpeed();
								Vector uspeed_eff2 = uspeed2 + ball->velocity;

								// if uspeed_eff passes 0
								scalar_t uspeed_eff_par = uspeed_eff.Dot(uspeed_eff - uspeed_eff2);
								scalar_t uspeed_eff2_par = uspeed_eff2.Dot(uspeed_eff - uspeed_eff2);

								if (Vector::ZERO().NDist2(uspeed_eff, uspeed_eff2) <= SlideThreshSpeed &&
									((uspeed_eff_par > scalar_t(0) && uspeed_eff2_par < scalar_t(0)) || (uspeed_eff2_par > scalar_t(0) && uspeed_eff_par < scalar_t(0)))
								) {
									// make rolling if uspeed_eff passed 0
									ball->velocity = ball->angularVelocity.Cross(Vector {scalar_t(0), scalar_t(0), ball->radius});
								}

							} else { // if rolling
								Vector waccel;
								{
									scalar_t roll_mom_r = MuRoll * ball->massMom / ball->mass / ball->diameter;
									Vector rollmom = Vector {scalar_t(0), scalar_t(0), ball->mass * Gravity * roll_mom_r}.Cross(ball->velocity.Unit());
									waccel = rollmom * -ball->invMassMom;

									// Spin deceleration rate
									waccel.z = SpotR * scalar_t(int(ball->angularVelocity.z < scalar_t(0)) - int(scalar_t(0) < ball->angularVelocity.z));

									waccel *= dt;
									if (waccel.z * waccel.z > ball->angularVelocity.z * ball->angularVelocity.z) {
										waccel.z = -ball->angularVelocity.z;
									}
								}
								scalar_t ws = ball->angularVelocity.LengthSqr();

								ball->angularVelocity += waccel;

								// Rolling on a flat surface should never cause angular velocity to increase,
								// so if waccel increased the angular velocity (due to gravity force being used), it's time to zero.
								if (ball->angularVelocity.LengthSqr() > ws) {
									ball->angularVelocity = Vector::ZERO();
								}

								// Align velocity with angularVelocity to assure rolling
								ball->velocity = ball->angularVelocity.Cross(Vector {scalar_t(0), scalar_t(0), ball->radius});
							}
						}

					} else {
						if (ball->position.z <= scalar_t(0)) {
							// Ball is completely inside the pocket now.
							ball->inPocket = pocket;

						} else if (ball->position.z < ball->radius) {
							// Check ball collision with the pocket edge
							Vector dr = ball->position - pocket->position;
							dr.z = scalar_t(0);
							dr = dr.Unit() * pocket->radius;

							if (dr != Vector::ZERO()) {
								dr += pocket->position; // dr is a point in pocket edge

								// Check if the ball is intersecting the collision point
								if ((ball->position - dr).LengthSqr() < ball->radius * ball->radius) {
									Vector normal = (ball->position - dr).Unit();
									ball->velocity -= ball->velocity.Proj(normal);
									ball->position = dr + normal * ball->radius;
								}
							}
						}
					}

				} else {
					// Bounce the ball inside the pocket
					Vector dr = ball->position - pocket->position;
					dr.z = scalar_t(0);

					scalar_t rs = (pocket->radius - ball->radius);
					rs *= rs;

					if (dr.LengthSqr() > rs) {
						scalar_t z = ball->velocity.z;
						ball->velocity -= ball->velocity.Proj(ball->position - pocket->position) * scalar_t(2);
						ball->velocity.z = z;

						dr = dr.Unit();
						z = ball->position.z;
						ball->position = pocket->position + (dr * pocket->radius) - (dr * ball->radius);
						ball->position.z = z;
					}

					// Clamp balls inside the pocket
					scalar_t minz = -pocket->depth + ball->radius;
					if (ball->position.z < minz) {
						ball->position.z = minz;
						ball->velocity = Vector::ZERO();
						ball->angularVelocity = Vector::ZERO();
						ball->enabled = false;
						ball->OnPocketed(pocket);
					}
				}
			}

			// Gravity
			if (!onGround || pocket != nullptr || !IsBallInsidePlayfield(ball)) {
				ball->velocity.z += -Gravity * dt;
			}

			// Correct intersecting balls
			ApplyContactThreshold(ball);
		}

		return balls_moving;
	}
}
