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

#define BILLIARDPHYSICS_ADV_EDGE_COLLISION

namespace BilliardPhysics
{
	static const scalar_t SQRTM1 {100000};

	static scalar_t plane_dist(Vector r, Vector rp, Vector n)
	{
		return (r - rp).Dot(n);
	}

	// ==========================================================================================
	Ball::Ball() :
		inPocket(nullptr),
		position(Vector::ZERO),
		velocity(Vector::ZERO),
		angularVelocity(Vector::ZERO),
		enabled(true)
	{
		Define(scalar_t(285) / scalar_t(10000), scalar_t(17) / scalar_t(100));
	}

	Ball::~Ball() = default;

	void Ball::Define(scalar_t radius_, scalar_t mass_)
	{
		radius = radius_;
		mass = mass_;
		massMom = scalar_t(2) * mass * radius * radius / scalar_t(3);
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

	// ==========================================================================================
	Pocket::~Pocket() = default;

	// ==========================================================================================
	Collider::~Collider() = default;

	void Collider::Update()
	{
		bbox = {
			Vector {scalar_t_max_value(), scalar_t_max_value(), scalar_t_max_value()},
			Vector {scalar_t_min_value(), scalar_t_min_value(), scalar_t_min_value()}
		};

		for (Shape& shape : shapes) {
			switch (shape.type) {
				case Shape::Type::Triangle:
					bbox.min.x = std::min({bbox.min.x, shape.r1.x, shape.r2.x, shape.r3.x});
					bbox.min.y = std::min({bbox.min.y, shape.r1.y, shape.r2.y, shape.r3.y});
					bbox.min.z = std::min({bbox.min.z, shape.r1.z, shape.r2.z, shape.r3.z});
					bbox.max.x = std::max({bbox.max.x, shape.r1.x, shape.r2.x, shape.r3.x});
					bbox.max.y = std::max({bbox.max.y, shape.r1.y, shape.r2.y, shape.r3.y});
					bbox.max.z = std::max({bbox.max.z, shape.r1.z, shape.r2.z, shape.r3.z});
					break;
				case Shape::Type::Line:
					bbox.min.x = std::min({bbox.min.x, shape.r1.x, shape.r2.x});
					bbox.min.y = std::min({bbox.min.y, shape.r1.y, shape.r2.y});
					bbox.min.z = std::min({bbox.min.z, shape.r1.z, shape.r2.z});
					bbox.max.x = std::max({bbox.max.x, shape.r1.x, shape.r2.x});
					bbox.max.y = std::max({bbox.max.y, shape.r1.y, shape.r2.y});
					bbox.max.z = std::max({bbox.max.z, shape.r1.z, shape.r2.z});
					break;
				case Shape::Type::Point:
					bbox.min.x = std::min(bbox.min.x, shape.r1.x);
					bbox.min.y = std::min(bbox.min.y, shape.r1.y);
					bbox.min.z = std::min(bbox.min.z, shape.r1.z);
					bbox.max.x = std::max(bbox.max.x, shape.r1.x);
					bbox.max.y = std::max(bbox.max.y, shape.r1.y);
					bbox.max.z = std::max(bbox.max.z, shape.r1.z);
					break;
			}
		}
	}

	// ==========================================================================================
	Engine::Engine()
	{
		Gravity = scalar_t(981) / scalar_t(100);
		MuRoll = scalar_t(1) / scalar_t(100);
		MuSlide = scalar_t(2) / scalar_t(10);
		MuBall = scalar_t(5) / scalar_t(100);

		SlideThreshSpeed = scalar_t(5) / scalar_t(1000);
		SpotR = sqrt(scalar_t(5));

		//OmegaMin = scalar_t(1) / scalar_t(10);
		//AirResistance = scalar_t(1) / scalar_t(1000);

		ContactThreshold = scalar_t_epsilon();

		fieldProperty.mu = scalar_t(2) / scalar_t(10);
		fieldProperty.loss0 = scalar_t(6) / scalar_t(10);
		fieldProperty.loss_max = scalar_t(80) / scalar_t(100);
		fieldProperty.loss_wspeed = scalar_t(2);
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

	bool Engine::IsBallInColliderRange(const Ball* ball, const Collider::Shape* shape) const
	{
		switch (shape->type) {
			case Collider::Shape::Type::Triangle:
			{
				Vector dr1 = shape->r2 - shape->r1;
				Vector dr2 = shape->r3 - shape->r2;
				Vector dr3 = shape->r1 - shape->r3;
				Vector n = dr1.Cross(dr2).Unit();
				return (
					plane_dist(ball->position, shape->r1, n.Cross(dr1).Unit()) >= scalar_t(0) &&
					plane_dist(ball->position, shape->r2, n.Cross(dr2).Unit()) >= scalar_t(0) &&
					plane_dist(ball->position, shape->r3, n.Cross(dr3).Unit()) >= scalar_t(0)
				);
			}
			case Collider::Shape::Type::Line:
			{
				Vector r = ball->position - shape->r1;
				Vector dr = shape->r2 - shape->r1;
				scalar_t dra = dr.Length();
				scalar_t d = r.Dot(dr) / dra;
				return (d >= scalar_t(0) && d < dra);
			}
			case Collider::Shape::Type::Point:
			{
				return true;
			}
		}
		return true;
	}

	scalar_t Engine::BallDistance(const Ball* ball, const Collider::Shape* shape) const
	{
		if (IsBallInColliderRange(ball, shape)) {
			switch (shape->type) {
				case Collider::Shape::Type::Triangle:
				{
					return (ball->position - shape->r1).Dot(shape->normal);
				}
				case Collider::Shape::Type::Line:
				{
					Vector r = ball->position - shape->r1;
					Vector dr = shape->r2 - shape->r1;
					return (r - r.Proj(dr)).Length();
				}
				case Collider::Shape::Type::Point:
				{
					return (ball->position - shape->r1).Length();
				}
			}
		}
		return scalar_t_infinity(); //old return = -1.0E20
	}

	bool Engine::BallCollided(const Ball* ball, const Collider::Shape* shape, scalar_t dt) const
	{
		switch (shape->type) {
			case Collider::Shape::Type::Triangle:
			{
				return (ball->velocity.Dot(shape->normal) > scalar_t(0));
			}
			case Collider::Shape::Type::Line:
			{
				Vector ballpos = ball->position + (ball->velocity * dt);
				return (ball->velocity.Dot((ballpos - shape->r1).NComp(shape->r2 - shape->r1)) > scalar_t(0));
			}
			case Collider::Shape::Type::Point:
			{
				Vector ballpos = ball->position + (ball->velocity * dt);
				return ((ballpos - shape->r1).Dot(ball->velocity) > scalar_t(0));
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

	Vector Engine::PerimeterSpeed(const Ball* ball) const
	{
		return ball->angularVelocity.Cross(Vector {scalar_t(0), scalar_t(0), -ball->radius});
	}

	Vector Engine::PerimeterSpeedNormal(const Ball* ball, const Vector& normal) const
	{
		return ball->angularVelocity.Cross(normal * -ball->radius);
	}

	scalar_t Engine::CalcCollisionTime(const Ball* b1, const Ball* b2) const
	{
		Vector dv = b1->velocity - b2->velocity;
		Vector dr = b1->position - b2->position;
		scalar_t vs = dv.x * dv.x + dv.y * dv.y + dv.z * dv.z;
		if (vs == scalar_t(0)) {
			return scalar_t_infinity();
		}

		scalar_t rs = dr.x * dr.x + dr.y * dr.y + dr.z * dr.z;
		scalar_t ds = b1->radius + b2->radius;
		ds *= ds;

		scalar_t p = (dv.x * dr.x + dv.y * dr.y + dv.z * dr.z) / vs;
		scalar_t q = (rs - ds) / vs;
		q = (p * p > q) ? sqrt(p * p - q) : SQRTM1;
		scalar_t t1 = -p + q;
		scalar_t t2 = -p - q;

		return (t1 < t2) ? t1 : t2;
	}

	scalar_t Engine::CalcCollisionTime(const Ball* ball, const Collider::Shape* shape) const
	{
		constexpr scalar_t inf = scalar_t_infinity();
		scalar_t rval {};

		switch (shape->type) {
			case Collider::Shape::Type::Triangle:
			{
				Vector dr = ball->position - shape->r1;
				scalar_t h = dr.Dot(shape->normal) - ball->radius;
				scalar_t vn = ball->velocity.Dot(shape->normal);
				if (vn == scalar_t(0)) {
					rval = inf;
					break;
				}
				rval = -h / vn;
				break;
			}
			case Collider::Shape::Type::Line:
			{
				// del all comps par to cylinder
				Vector dr = shape->r2 - shape->r1;
				Vector r = ball->position - shape->r1;
				r -= r.Proj(dr);
				Vector v = ball->velocity;
				v -= v.Proj(dr);

				scalar_t vls = v.LengthSqr();
				if (vls == scalar_t(0)) {
					rval = inf;
					break;
				}

				scalar_t ph = v.Dot(r) / vls;
				scalar_t q = (r.LengthSqr() - (ball->radius * ball->radius)) / vls;

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
			case Collider::Shape::Type::Point:
			{
				scalar_t vls = ball->velocity.LengthSqr();
				if (vls == scalar_t(0)) {
					rval = inf;
					break;
				}

				Vector r = ball->position - shape->r1;
				scalar_t ph = ball->velocity.Dot(r) / vls;
				scalar_t q = (r.LengthSqr() - (ball->radius * ball->radius)) / vls;

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
		for (Ball* ball : balls) {
			if (!ball->enabled || (ball->velocity.LengthSqr() + ball->angularVelocity.LengthSqr()) == scalar_t(0)) {
				continue;
			}

			// Pre-translate bbox
			Vector vr {ball->radius, ball->radius, ball->radius};
			ball->bbox.min = ball->position - vr;
			ball->bbox.max = ball->position + vr;

			// Translate ball
			Vector dx = ball->velocity * dt;
			ball->position += dx;

			// Post-translate bbox
			ball->bbox.Merge(BoundingBox {
				ball->position - vr,
				ball->position + vr
			});

			// Rotate ball
			ball->ApplyRotation(dt);
		}
	}

	void Engine::CollectColliders(const Ball* ball)
	{
		colliders.clear();

		// Check cushion collisions
		for (const Collider* collider : fieldColliders) {
			if (collider->bbox.Intersects(ball->bbox)) {
				colliders.push_back(collider);
			}
		}

		scalar_t ground = ball->position.z - ball->radius;

		// Balls in air check additional collisions
		if (ground != scalar_t(0) || ball->velocity.z != scalar_t(0)) {
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

		switch (shape->type) {
			case Collider::Shape::Type::Triangle:
			{
				hit_normal = shape->normal;
				break;
			}
			case Collider::Shape::Type::Line:
			{
				Vector dr = shape->r2 - shape->r1;
				hit_normal = ball->position - shape->r1;
				hit_normal = (hit_normal - hit_normal.Proj(dr)).Unit();
				break;
			}
			case Collider::Shape::Type::Point:
			{
				hit_normal = (ball->position - shape->r1).Unit();
				break;
			}
		}

		const Collider::Property* prop = (collider && collider->property) ? collider->property : &fieldProperty;
		Vector vn = ball->velocity.Proj(hit_normal);
		Vector vp = ball->velocity - vn;

		// normal component
		scalar_t loss = prop->loss0 + (prop->loss_max - prop->loss0) * (scalar_t(1) - exp(-vn.Length() / prop->loss_wspeed));
		Vector dv = vn * -(scalar_t(1) + sqrt(scalar_t(1) - loss));
		ball->velocity += dv;

		// parallel component
		scalar_t diameter = ball->radius * scalar_t(2);
		Vector ps = PerimeterSpeedNormal(ball, hit_normal);
		dv = (ps + vp).Unit() * (-dv.Length() * prop->mu);
		Vector dw = (dv * (ball->mass / scalar_t(2) / ball->massMom)).Cross(hit_normal * diameter);
		Vector dw2 = dw + ball->angularVelocity.Proj(dw);

		if (dw2.Dot(ball->angularVelocity) < scalar_t(0)) {
			dw = dw - dw2;
			dv = dv.Unit() * (dw.Length() * scalar_t(2) * ball->massMom / ball->mass / diameter);
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
		Vector b1v = Vector::ZERO;

		Vector dvn = duvec * duvec.Dot(b2v);
		Vector dvp = b2v - dvn;

		b2v = b2v - dvn;
		b2v += dvn * ((b2->mass - b1->mass) / (b1->mass + b2->mass));
		b1v = dvn * (scalar_t(2) * b2->mass / (b1->mass + b2->mass)); // (momentum transfer)/m1
		b2->velocity = b1->velocity + b2v;
		b1->velocity = b1->velocity + b1v;

		// angular momentum transfer
		Vector dpn = b1v * b1->mass; // momentum transfer from ball2 to ball1
		Vector perimeter_speed_b1 = b1->angularVelocity.Cross(duvec * -b1->radius);
		Vector perimeter_speed_b2 = b2->angularVelocity.Cross(duvec * b2->radius);
		Vector fric_dir = ((perimeter_speed_b2 - perimeter_speed_b1) + dvp).Unit();
		Vector dpp = fric_dir * (-dpn.Length() * MuBall); // dp parallel of ball2

		Vector dw2 = dpp.Cross(duvec) * (b2->radius * scalar_t(2) / b2->massMom);
		Vector dw1 = dpp.Cross(duvec) * (b1->radius * scalar_t(2) / b1->massMom);
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

	// Ball interaction with table playfield
	void Engine::BallInteraction(Ball* ball)
	{
		const Collider::Property& prop = fieldProperty;

		Vector hit_normal {scalar_t(0), scalar_t(0), scalar_t(1)};
		Vector vn = ball->velocity.Proj(hit_normal);
		Vector vp = ball->velocity - vn;

		// normal component
		scalar_t loss = prop.loss0 + (prop.loss_max - prop.loss0) * (scalar_t(1) - exp(-vn.Length() / prop.loss_wspeed));
		Vector dv = vn * -(scalar_t(1) + sqrt(scalar_t(1) - loss));
		ball->velocity += dv;

		// parallel component
		scalar_t diameter = ball->radius * scalar_t(2);
		Vector ps = PerimeterSpeedNormal(ball, hit_normal);
		dv = (ps + vp).Unit() * (-dv.Length() * prop.mu);
		Vector dw = (dv * (ball->mass / scalar_t(2) / ball->massMom)).Cross(hit_normal * diameter);
		Vector dw2 = dw + ball->angularVelocity.Proj(dw);

		if (dw2.Dot(ball->angularVelocity) < scalar_t(0)) {
			dw = dw - dw2;
			dv = dv.Unit() * (dw.Length() * scalar_t(2) * ball->massMom / ball->mass / diameter);
		}

		ball->angularVelocity += dw;
		ball->velocity += dv;
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
				dist -= dr.Length();
				ball->position += dr.Unit() * (dist + ContactThreshold);
			}
		}
	}

	void Engine::StepSimulationEuler(scalar_t dt, int depth)
	{
		scalar_t dt1 {};
		scalar_t dtmin {};

		Ball* colBall = nullptr;
		Ball* colBall2 = nullptr;
		const Collider* col = nullptr;
		const Collider::Shape* colShape = nullptr;

#ifdef BILLIARDPHYSICS_ADV_EDGE_COLLISION
		const Pocket* colPocket = nullptr;

		Collider::Shape edgePoint;
		edgePoint.type = Collider::Shape::Type::Point;
#endif

		MoveBalls(dt);

		// checks
		for (Ball* ball : balls) {
			if (!ball->enabled || ball->IsPocketed()) {
				continue;
			}

			CollectColliders(ball);

			// Check shape collisions
			for (const Collider* collider : colliders) {
				for (const Collider::Shape& shape : collider->shapes) {
					dt1 = CalcCollisionTime(ball, &shape);
					if (dt1 < dtmin && dt1 > -dt && !BallCollided(ball, &shape, -dt)) {
						// dont strobe apart
						dtmin = dt1;
						colBall = ball;
						col = collider;
						colShape = &shape;
					}
				}
			}

			// Check ball collisions
			for (Ball* ball2 : balls) {
				if (ball2 == ball || !ball2->enabled) {
					continue;
				}

				dt1 = CalcCollisionTime(ball, ball2);
				if (dt1 < dtmin && dt1 > -dt && !BallCollided(ball2, ball, -dt)) {
					// dont strobe apart
					dtmin = dt1;
					colBall = ball2;
					colBall2 = ball;
					col = nullptr;
				}
			}

#ifdef BILLIARDPHYSICS_ADV_EDGE_COLLISION
			// Check pocket edge collision
			const Pocket* pocket = GetPocketBallInside(ball);
			if (pocket) {
				Vector dr = ball->position - pocket->position;
				dr.z = scalar_t(0);

				// Get a collision point at the pocket edge
				edgePoint.r1 = pocket->position + dr.Unit() * pocket->radius;
				edgePoint.r1.z = pocket->position.z;

				dt1 = CalcCollisionTime(ball, &edgePoint);
				if (dt1 < dtmin && dt1 > -dt && !BallCollided(ball, &edgePoint, -dt)) {
					// dont strobe apart
					dtmin = dt1;
					colBall = ball;
					colBall2 = nullptr;
					col = nullptr;
					colPocket = pocket;
				}
			}
#endif

		}

		if (col) {
			MoveBalls(dtmin);
			BallInteraction(colBall, colShape, col);
			colBall->OnCollided(colShape, col, -dtmin);

			StepSimulationEuler(-dtmin, depth + 1);

		} else if (colBall && colBall2) {
			MoveBalls(dtmin);
			BallInteraction(colBall, colBall2);
			colBall->OnCollided(colBall2, -dtmin);

			StepSimulationEuler(-dtmin, depth + 1);
		}
#ifdef BILLIARDPHYSICS_ADV_EDGE_COLLISION
		else if (colBall && colPocket) {
			MoveBalls(dtmin);
			BallInteraction(colBall, &edgePoint, nullptr);
			colBall->OnCollided(&edgePoint, nullptr, -dtmin);

			StepSimulationEuler(-dtmin, depth + 1);
		}
#endif

	}

	int Engine::StepSimulation(scalar_t dt)
	{
		int balls_moving = 0;

		// timestep with actual speeds, omegas,...
		StepSimulationEuler(dt, 0);

		// Calc new accelerations and speeds
		for (Ball* ball : balls) {
			if (!ball->enabled) {
				continue;
			}

			// check if balls still moving
			if (ball->velocity.LengthSqr() != scalar_t(0) || ball->angularVelocity.LengthSqr() != scalar_t(0)) {
				++balls_moving;
			}

			scalar_t ground = ball->position.z - ball->radius;
			Pocket* pocket = ball->inPocket;

			if (!pocket) {
				// Ball is not fully inside a pocket, but may be hovering one.
				pocket = GetPocketBallInside(ball);

				if (pocket == nullptr) {
					// absolute and relative perimeter speed
					Vector uspeed = PerimeterSpeed(ball);
					Vector uspeed_eff = uspeed + ball->velocity;

					// only if ball not flying do sliding/rolling
					if (ground == scalar_t(0)) {
						// if sliding
						if (uspeed_eff.Length() > SlideThreshSpeed) {
							// acc caused by friction
							Vector fricaccel = uspeed_eff.Unit() * (-MuSlide * Gravity);

							// angular acc caused by friction
							Vector fricmom = fricaccel.Cross(Vector {scalar_t(0), scalar_t(0), -ball->radius}) * ball->mass;
							Vector waccel = fricmom * -ball->invMassMom;

							// perform accel
							ball->angularVelocity += waccel * dt;
							ball->velocity += fricaccel * dt;
							Vector uspeed2 = PerimeterSpeed(ball);
							Vector uspeed_eff2 = uspeed2 + ball->velocity;

							// if uspeed_eff passes 0
							scalar_t uspeed_eff_par = uspeed_eff.Dot(uspeed_eff - uspeed_eff2);
							scalar_t uspeed_eff2_par = uspeed_eff2.Dot(uspeed_eff - uspeed_eff2);

							if (Vector::ZERO.NDist(uspeed_eff, uspeed_eff2) <= SlideThreshSpeed &&
								((uspeed_eff_par > scalar_t(0) && uspeed_eff2_par < scalar_t(0)) || (uspeed_eff2_par > scalar_t(0) && uspeed_eff_par < scalar_t(0)))
							) {
								// make rolling if uspeed_eff passed 0
								ball->velocity = ball->angularVelocity.Cross(Vector {scalar_t(0), scalar_t(0), ball->radius});
							}

						} else { // if rolling
							scalar_t diameter = ball->radius * scalar_t(2);

							Vector waccel;
							{
								scalar_t roll_mom_r = MuRoll * ball->massMom / ball->mass / diameter;
								Vector rollmom = Vector {scalar_t(0), scalar_t(0), ball->mass * Gravity * roll_mom_r}.Cross(ball->velocity.Unit());
								waccel = rollmom * -ball->invMassMom;
							}

							scalar_t wls = ball->angularVelocity.LengthSqr();
							scalar_t vls = ball->velocity.LengthSqr();

							// Independently update spin component
							{
								scalar_t alpha = SpotR * ball->radius * Gravity / diameter;
								alpha *= dt;
								if (abs(ball->angularVelocity.z) < alpha) {
									ball->angularVelocity.z = scalar_t(0);
								} else {
									ball->angularVelocity.z += (ball->angularVelocity.z > scalar_t(0)) ? -alpha : alpha;
								}
							}

							ball->angularVelocity += waccel * dt;

							// align velocity with angularVelocity to assure rolling
							ball->velocity = ball->angularVelocity.Cross(Vector {scalar_t(0), scalar_t(0), ball->radius});

							// Check for low velocities
							if (ball->angularVelocity.LengthSqr() > wls && ball->velocity.LengthSqr() > vls) {
								ball->velocity = Vector::ZERO;
								ball->angularVelocity = Vector::ZERO;
							}
						}
					}

					// Ball collided with the table plane
					if (ground < scalar_t(0)) {
						BallInteraction(ball);

						if (ball->velocity.z < Gravity * dt * scalar_t(2)) {
							ball->velocity.z = scalar_t(0);
						}

						// Reposition the ball to be on the table plane.
						// IMPROVE: use better algorithm to find the exact position (not just z) the ball should be placed.
						ball->position.z = ball->radius;
						ground = scalar_t(0);
					}

				} else {
#ifndef BILLIARDPHYSICS_ADV_EDGE_COLLISION
					// Check ball collision with the pocket edge
					if (ground > -ball->radius) {
						Vector dr = ball->position - pocket->position;
						dr.z = scalar_t(0);

						scalar_t rs = (pocket->radius - ball->radius);
						rs *= rs;
						if (dr.LengthSqr() > rs) {
							// Get a collision point at the pocket edge
							dr = pocket->position + dr.Unit() * pocket->radius;
							dr.z = pocket->position.z;

							// Check if the ball is intersecting the collision point
							if ((ball->position - dr).LengthSqr() < ball->radius * ball->radius) {
								Vector normal = (ball->position - dr).Unit();
								ball->velocity -= ball->velocity.Proj(normal);
								ball->position = dr + normal * ball->radius;
							}
						}

					} else {
						// Ball is completely inside the pocket now.
						ball->inPocket = pocket;
					}
#else
					if (ball->position.z < scalar_t(0)) {
						// Ball is completely inside the pocket now.
						ball->inPocket = pocket;
					}
#endif
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
					ball->velocity = Vector::ZERO;
					ball->angularVelocity = Vector::ZERO;
					ball->enabled = false;
					ball->OnPocketed(pocket);
				}
			}

			// Gravity
			if (ground != scalar_t(0) || pocket != nullptr) {
				ball->velocity.z += -Gravity * dt;
			}

			// Correct intersecting balls
			ApplyContactThreshold(ball);
		}

		return balls_moving;
	}
}
