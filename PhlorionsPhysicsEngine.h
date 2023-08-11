#pragma once
#include <iostream>
#include <vector>
#include <math.h>
#include <limits.h>
#include <assert.h>

namespace phlorion
{
	// Solve a quadratic equation
	static std::vector<float> quadratic(const float a, const float b, const float c)
	{
		float D = b * b - 4 * a * c;
		if (D > 0)
		{
			std::vector<float> res = { (-1 * b + sqrtf(D)) / (2 * a), (-1 * b - sqrtf(D)) / (2 * a) };
			return res;
		}
		else if (D == 0)
		{
			std::vector<float> res = { (-1 * b) / (2 * a) };
			return res;
		}
	}

	/**
	* Holds a vector in 3 dimensions.
	* Four data members are allocated to ensure alignment in an array
	*/
	class Vector3
	{
	private:
		// Padding to ensure 4-word alignment
		float pad;

	public:
		float x;
		float y;
		float z;

		// default set all components to 0
		Vector3() : x(0), y(0), z(0) {};
		// The explicit constructor creates a vector with the given components
		Vector3(const float x, const float y, const float z) : x(x), y(y), z(z) {};
		// Deconstructor
		~Vector3() {};

		// Zero all the components of the vector
		void clear()
		{
			x = y = z = 0;
		}

		// invert the vector
		void invert()
		{
			x = -x;
			y = -y;
			z = -z;
		}

		// find the magnitude of the vector
		float magnitude()
		{
			return sqrtf(x * x + y * y + z * z);
		}

		// find the squared magnitude
		float sq_magnitude()
		{
			return x * x + y * y + z * z;
		}

		// Turns a non-zero vector into a vector of unit length
		void normalize()
		{
			float l = magnitude();
			if (l > 0)
			{
				*this = *this * (float)(1 / l);
			}
		}

		// Adds the given vector to this, scaled by the given amount
		void addScaledVector(const Vector3& vec, float scale)
		{
			x = x + scale * vec.x;
			y = y + scale * vec.y;
			z = z + scale * vec.z;
		}

		// Calculates and returns the scalar product of this vector with the given vector
		float scalarProduct(const Vector3& vec) const
		{
			return x * vec.x + y * vec.y + z * vec.z;
		}

		// Calculates and returns the vector product of this vector with the given vector
		Vector3 vectorProduct(const Vector3& vec) const
		{
			return Vector3(y * vec.z - z * vec.y, z * vec.x - x * vec.z, x * vec.y - y * vec.x);
		}

		/*
		-- Vector multiplication with constant --
		Operator overload *
		We want to multiply a vector with a constant
		Returns a new Vector3 multiplied by the constant
		*/
		Vector3 operator*(const float& value) const
		{
			return Vector3(x*value, y*value, z*value);
		}

		/*
		-- Vector addition --
		Operator overload +
		We want to add two vectors
		Returns a new Vector3, the result of the sum between vector a and b
		*/
		Vector3 operator+(const Vector3& vec) const
		{
			return Vector3(x + vec.x, y + vec.y, z + vec.z);
		}

		/*
		-- Vector subtraction --
		Operator overload -
		We want to subtract two vectors
		Returns a new Vector3, the result of the subtraction between vector a and b
		*/
		Vector3 operator-(const Vector3& vec) const
		{
			return Vector3(x - vec.x, y - vec.y, z - vec.z);
		}

	};

	/**
	* A particle is the simplest object that can be simulated in the
	* physics system.
	*/
	class Particle
	{
	protected:
		// Holds the linear position of the particle in world space
		Vector3 position;
		//  Holds the linear velocity of the particle in world space
		Vector3 velocity;
		//  Holds the acceleration of the particle. This value can be used to set acceleration due to gravity(its primary use) or any other constant acceleration
		Vector3 acceleration;
		Vector3 last_acceleration;

		// Holds the accumulated force to be applied at the next simulation iteration only. This value is zeroed at each integration step
		Vector3 forceAccum;

		// Holds the amount of damping applied to linear motion. Damping is required to remove energy added through numerical instability in the integrator
		float damping;

		/**
		* Holds the inverse of the mass of the particle. It is more useful to hold the inverse mass because integration is simpler and because in real-time
		* simulation it is more useful to have objects with infinite mass (immovable) than zero mass (completely unstable in numerical simulation).
		*/
		float inverseMass;

	public:
		/**
		* Default Constructor
		* Sets position = velocity = acceleration = forceAccum = 0, damping ~ 1, inverseMass = 1
		*/
		Particle() : position(Vector3()), velocity(Vector3()), acceleration(Vector3()), forceAccum(Vector3()), damping(0.991f), inverseMass(1.f) {};
		/**
		* Explicit Constructor
		* Sets position, velocity, acceleration, forceAccum, damping, mass to given components (if mass != 0)
		*/
		Particle(const Vector3& pos, const Vector3& vel, const Vector3& acc, const Vector3& force, const float damp, const float mass) :
			position(pos), velocity(vel), acceleration(acc), forceAccum(force), damping(damp) { assert(mass != 0); inverseMass = ((float)1.f) / mass; };
		// Deconstructor
		~Particle() {};

		/**
		* Sets the mass of the particle.
		*
		* @param mass The new mass of the body. This may not be zero. Small masses can produce unstable rigid bodies under simulation.
		*
		* @warning This invalidates internal data for the particle. Either an integration function, or the calculateInternals
		* function should be called before trying to get any settings from the particle.
		*/
		void setMass(const float mass)
		{
			assert(mass != 0);
			inverseMass = ((float)1.f) / mass;
		}

		/**
		* Gets the mass of the particle.
		*
		* @return The current mass of the particle.
		*/
		float getMass() const
		{
			if (inverseMass == 0) {
				return LONG_MAX;
			}
			else {
				return ((float)1.f) / inverseMass;
			}
		}

		/**
		* Sets the inverse mass of the particle.
		*
		* @param inverseMass The new inverse mass of the body. This may be zero, for a body with infinite mass (i.e. unmovable).
		*
		* @warning This invalidates internal data for the particle. Either an integration function, or the calculateInternals
		* function should be called before trying to get any settings from the particle.
		*/
		void setInverseMass(const float invMass)
		{
			inverseMass = invMass;
		}

		/**
		* Gets the inverse mass of the particle.
		*
		* @return The current inverse mass of the particle.
		*/
		float getInverseMass() const
		{
			return inverseMass;
		}

		/**
		 * Returns true if the mass of the particle is not-infinite.
		 */
		bool hasFiniteMass() const
		{
			return inverseMass >= 0.0f;
		}

		/**
		* Sets the damping of the particle.
		*/
		void setDamping(const float damp)
		{
			damping = damp;
		}

		/**
		* Gets the current damping value.
		*/
		float getDamping() const
		{
			return damping;
		}

		/**
		* Sets the position of the particle.
		*
		* @param position The new position of the particle.
		*/
		void setPosition(const Vector3& pos)
		{
			position = pos;
		}
		/**
		* Sets the position of the particle by component.
		*
		* @param x The x coordinate of the new position of the rigid body.
		* @param y The y coordinate of the new position of the rigid body.
		* @param z The z coordinate of the new position of the rigid body.
		*/
		void setPosition(const float x, const float y, const float z)
		{
			position.x = x;
			position.y = y;
			position.z = z;
		}

		/**
		* Gets the position of the particle.
		*
		* @return The position of the particle.
		*/
		Vector3 getPosition() const
		{
			return position;
		}

		/**
		* Fills the given vector with the position of the particle.
		*
		* @param position A pointer to a vector into which to write the position.
		*/
		void getPosition(Vector3* position) const
		{
			*position = Particle::position;
		}

		/**
		* Sets the velocity of the particle.
		*
		* @param velocity The new velocity of the particle.
		*/
		void setVelocity(const Vector3& vel)
		{
			velocity = vel;
		}
		/**
		* Sets the velocity of the particle by component.
		*
		* @param x The x coordinate of the new velocity of the rigid body.
		* @param y The y coordinate of the new velocity of the rigid body.
		* @param z The z coordinate of the new velocity of the rigid body.
		*/
		void setVelocity(const float x, const float y, const float z)
		{
			velocity.x = x;
			velocity.y = y;
			velocity.z = z;
		}

		/**
		* Gets the velocity of the particle.
		*
		* @return The velocity of the particle. The velocity is given in world local space.
		*/
		Vector3 getVelocity() const
		{
			return velocity;
		}

		/**
		* Fills the given vector with the velocity of the particle.
		*
		* @param position A pointer to a vector into which to write the position.
		*/
		void getVelocity(Vector3* velocity) const
		{
			*velocity = Particle::velocity;
		}

		/**
		* Sets the constant acceleration of the particle.
		*
		* @param acceleration The new acceleration of the particle.
		*/
		void setAcceleration(const Vector3& acc)
		{
			acceleration = acc;
		}
		/**
		* Sets the constant acceleration of the particle by component.
		*
		* @param x The x coordinate of the new acceleration of the rigid body.
		* @param y The y coordinate of the new acceleration of the rigid body.
		* @param z The z coordinate of the new acceleration of the rigid body.
		*/
		void setAcceleration(const float x, const float y, const float z)
		{
			acceleration.x = x;
			acceleration.y = y;
			acceleration.z = z;
		}

		/**
		* Gets the acceleration of the particle.
		*
		* @return The acceleration of the particle. The acceleration is given in world local space.
		*/
		Vector3 getAcceleration() const
		{
			return acceleration;
		}

		/**
		* Fills the given vector with the acceleration of the particle.
		*
		* @param position A pointer to a vector into which to write the position.
		*/
		void getAcceleration(Vector3* acceleration) const
		{
			*acceleration = Particle::acceleration;
		}

		Vector3 getLastAcceleration() const
		{
			return last_acceleration;
		}

		/**
		* Clears the forces applied to the particle. This will be called automatically after each integration step.
		*/
		void clearAccumulator()
		{
			forceAccum.clear();
		}

		/**
		* Adds the given force to the particle, to be applied at the next iteration only.
		*
		* @param force The force to apply.
		*/
		void addForce(const Vector3& force)
		{
			forceAccum = forceAccum + force;
		}

		/** 
		* Integrates the particle forward in time by the given amount. This function uses a Newton - Euler integration method, which is a
		* linear approximation of the correct integral. For this reason it may be inaccurate in some cases.
		*/
		void integrate(float t)
		{
			assert(t > 0.f);

			// Update linear position
			position.addScaledVector(velocity, t);

			// Work out the acceleration from the force
			Vector3 resultingAcc = acceleration;
			resultingAcc.addScaledVector(forceAccum, inverseMass);
			last_acceleration = resultingAcc;

			// Update linear velocity from the acceleration
			velocity.addScaledVector(resultingAcc, t);

			// Impose drag
			velocity = velocity * powf(damping, t);

			// Clear the forces.
			clearAccumulator();
		}
	};

	/**
	* Holds all the force generators and the particles they apply to.
	*/
	class ParticleForceRegistry
	{
	protected:
		/**
		* Keeps track of one force generator and the particle applies to.
		*/
		struct ParticleForceRegistration
		{
			Particle* particle;
			ParticleForceGenerator* fg;
		};

		/**
		* Holds the list of registrations.
		*/
		typedef std::vector<ParticleForceRegistration> Registry;
		Registry registrations;

	public:
		/**
		* Registers the given force generator to apply to the given particle.
		*/
		void add(Particle* particle, ParticleForceGenerator* fg)
		{
			ParticleForceRegistration PFR;
			PFR.particle = particle; PFR.fg = fg;
			registrations.push_back(PFR);
		}

		/**
		* Removes the given registered pair from the registry. If the pair is not registered, this method will have no effect.
		*/
		void remove(Particle* particle, ParticleForceGenerator* fg)
		{
			Registry::iterator i = registrations.begin();
			for (; i != registrations.end(); i++)
			{
				if (particle == i->particle && fg == i->fg)
				{
					registrations.erase(i);
				}
			}
		}

		/**
		* Clears all registrations from the registry. This will not delete the particles or the force generators themselves, just the records of their connection.
		*/
		void clear()
		{
			registrations.clear();
		}

		/**
		* Calls all the force generators to update the forces of their corresponding particles.
		*/
		void updateForces(float t)
		{
			Registry::iterator i = registrations.begin();
			for (; i != registrations.end(); i++)
			{
				i->fg->updateForce(i->particle, t);
			}
		}
	};

	/**
	* A force generator can be asked to add a force to one or more particles.
	*/
	class ParticleForceGenerator
	{
	public:
		/**
		* Overload this in implementations of the interface to calculate and update the force applied to the given particle.
		*/
		virtual void updateForce(Particle* particle, float t) = 0;
	};

	/**
	* A force generator that applies a gravitational force. One instance can be used for multiple particles.
	*/
	class ParticleGravity : public ParticleForceGenerator
	{
		/** Holds the acceleration due to gravity. */
		Vector3 gravity;

	public:
		/** Creates the generator with the given acceleration. */
		ParticleGravity(const Vector3& gravity) : gravity(gravity) {}

		/** Applies the gravitational force to the given particle. */
		virtual void updateForce(Particle* particle, float t)
		{
			// Check that we do not have infinite mass.
			if (!particle->hasFiniteMass()) return;

			// Apply the mass-scaled force to the particle.
			particle->addForce(gravity * particle->getMass());
		}
	};

	/**
	* A force generator that applies a drag force. One instance can be used for multiple particles.
	*/
	class ParticleDrag : public ParticleForceGenerator
	{
		/** Holds the velocity drag coefficient. */
		float k1;
		/** Holds the velocity squared drag coefficient. */
		float k2;

	public:
		/** Creates the generator with the given coefficients. */
		ParticleDrag(float k1, float k2) : k1(k1), k2(k2) {}

		/** Applies the drag force to the given particle. */
		virtual void updateForce(Particle* particle, float t)
		{
			Vector3 force;
			particle->getVelocity(&force);

			// Calculate the total drag coefficient
			float dragCoeff = force.magnitude();
			dragCoeff = k1 * dragCoeff + k2 * dragCoeff * dragCoeff;

			// Calculate the final force and apply it
			force.normalize();
			force = force * -dragCoeff;
			particle->addForce(force);
		}
	};

	/**
	* A force generator that applies a spring force.
	*/
	class ParticleSpring : public ParticleForceGenerator
	{
		/** The particle at the other end of the spring. */
		Particle* other;

		/** Holds the spring constant. */
		float springConstant;

		/** Holds the rest length of the spring. */
		float restLength;

	public:
		/** Creates a new spring with the given parameters. */
		ParticleSpring(Particle* other, float springConstant, float restLength) : 
			other(other), springConstant(springConstant), restLength(restLength) {}

		/**Applies the spring force to the given particle. */
		virtual void update(Particle* particle, float t)
		{
			// Calculate the vector of the spring.
			Vector3 force;
			particle->getPosition(&force);
			force = force - other->getPosition();

			// Calculate the magnitude of the force.
			float magnitude = force.magnitude();
			magnitude = fabsf(magnitude - restLength);
			magnitude = magnitude * springConstant;

			// Calculate the final force and apply it.
			force.normalize();
			force = force * -magnitude;
			particle->addForce(force);
		}
	};

	/**
	* A force generator that applies a spring force, where one end is attached to a fixed point in space.
	*/
	class ParticleAnchoredSpring : public ParticleForceGenerator
	{
		/** The location of the anchored end of the spring. */
		Vector3* anchor;

		/** Holds the spring constant. */
		float springConstant;

		/** Holds the rest length of the spring. */
		float restLength;

	public:
		/** Creates a new spring with the given parameters. */
		ParticleAnchoredSpring(Vector3* anchor, float springConstant, float restLength) :
			anchor(anchor), springConstant(springConstant), restLength(restLength) {}

		/** Applies the spring force to the given particle. */
		virtual void updateForce(Particle* particle, float t)
		{
			// Calculate the vector of the spring.
			Vector3 force;
			particle->getPosition(&force);
			force = force - *anchor;

			// Calculate the magnitude of the force.
			float magnitude = force.magnitude();
			magnitude = fabsf(magnitude - restLength);
			magnitude = magnitude * springConstant;

			// Calculate the final force and apply it.
			force.normalize();
			force = force * -magnitude;
			particle->addForce(force);
		}
	};

	/**
	* A force generator that applies a spring force only when extended.
	*/
	class ParticleBungee : public ParticleForceGenerator
	{
		/** The particle at the other end of the spring. */
		Particle* other;

		/** Holds the spring constant. */
		float springConstant;

		/**  Holds the length of the bungee at the point it begins to generate a force. */
		float restLength;

	public:
		/** Creates a new bungee with the given parameters. */
		ParticleBungee(Particle* other, float springConstant, float restLength) :
			other(other), springConstant(springConstant), restLength(restLength) {}

		/** Applies the spring force to the given particle. */
		virtual void updateForce(Particle* particle, float t)
		{
			// Calculate the vector of the spring.
			Vector3 force;
			particle->getPosition(&force);
			force = force - other->getPosition();

			// Check if the bungee is compressed.
			float magnitude = force.magnitude();
			if (magnitude <= restLength) return;

			// Calculate the magnitude of the force.
			magnitude = springConstant * (restLength - magnitude);

			// Calculate the final force and apply it.
			force.normalize();
			force = force * -magnitude;
			particle->addForce(force);

		}
	};

}
