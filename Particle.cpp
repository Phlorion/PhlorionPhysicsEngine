#include "Particle.h"

GraphicsParticle::GraphicsParticle() :
	Particle()
{
	radius = 0.5f;
	position.x = 6.f;
	position.y = 2.f;
	time_alive = 0.f;
	br.fill_color[0] = 1.f; br.fill_color[1] = 1.f; br.fill_color[2] = 1.f;

	// set trajectory buffer to positions that cannot exist
	for (int i = 0; i < tr_len; i++)
	{
		tr_buffer[i] = phlorion::Vector3(-1, -1, 0);
	}
}

GraphicsParticle::GraphicsParticle(const phlorion::Vector3& pos, const phlorion::Vector3& vel, const phlorion::Vector3& acc, const phlorion::Vector3& force,
	const float damp, const float mass, const float radius, const graphics::Brush br) :
	Particle(pos, vel, acc, force, damp, mass)
{
	this->radius = radius;
	time_alive = 0.f;
	this->br = br;

	// set trajectory buffer to positions that cannot exist
	for (int i = 0; i < tr_len; i++)
	{
		tr_buffer[i] = phlorion::Vector3(-1, -1, 0);
	}
}

GraphicsParticle::~GraphicsParticle() 
{
	delete[] tr_buffer;
}

void GraphicsParticle::setRadius(const float rad)
{
	radius = rad;
}

float GraphicsParticle::getRadius() const
{
	return radius;
}

void GraphicsParticle::update(float ms)
{
	time_alive += ms;

	/**
	* Trajectory calculation
	* Stores previous positions to draw the trajectory of the particle
	*/
	// when trajectory buffer is empty this while sets the first positions of the trajectory
	if (can_draw_trajectory)
	{
		int i = 0;
		while (i < tr_len)
		{
			if (tr_buffer[i].x == -1)
			{
				tr_buffer[i] = position;
				break;
			}
			i++;
		}
		// when the trajectory buffer gets filled, replace the new position with the index position
		if (i == tr_len)
		{
			tr_buffer[index] = position;
			can_draw_trajectory = true;
			index++;
		}
		// reset index
		if (index == tr_len) index = 0;
	}

	/**
	* Basic collision with wall
	*/
	if (position.y >= CANVAS_HEIGHT - radius && getVelocity().y > 0)
	{
		float prev_pos_y = getPosition().y;
		//position.y = CANVAS_HEIGHT - radius;
		
		velocity.y = getVelocity().y * -1;

		/**
		* Debug on collision with ground particle can go below the ground (due to update rate). Return the particle's position at exactly ground level
		* and calculate all differences
		* 
		* std::cout << "Collision position correction: " << prev_pos_y << " ";
		* std::cout << position.y << std::endl;
		* std::cout << "Time alive: " << time_alive << std::endl;
		* 
		*/
		std::cout << "Velocity on impact: " << getVelocity().y << std::endl;
	}
	else if (position.y <= radius && getVelocity().y < 0)
	{
		velocity.y = getVelocity().y * -1;
		std::cout << "Velocity on impact: " << getVelocity().y << std::endl;
	}
	else if (position.x >= CANVAS_WIDTH - radius && getVelocity().x > 0)
	{
		velocity.x = getVelocity().x * -1;
		std::cout << "Velocity on impact: " << getVelocity().x << std::endl;
	}
	else if (position.x <= radius && getVelocity().x < 0)
	{
		velocity.x = getVelocity().x * -1;
		std::cout << "Velocity on impact: " << getVelocity().x << std::endl;
	}

	// integration step, update particles position, velocity, acceleration, etc.
	if (ms != 0) integrate(ms / 1000);
}

void GraphicsParticle::draw()
{
	br.fill_color[0] = br.fill_color[1] = br.fill_color[2] = 1.f;
	graphics::drawDisk(position.x, position.y, radius, br);

	// draw trajectory
	if (can_draw_trajectory)
	{
		for (int i = 0; i < tr_len; i++)
		{
			// if we are at the first posiiton of the buffer, draw a line between this and the last position of the buffer
			if (i == 0 && i - 1 != index - 1 && tr_buffer[i].x != -1)
				graphics::drawLine(tr_buffer[tr_len - 1].x, tr_buffer[tr_len - 1].y, tr_buffer[i].x, tr_buffer[i].y, br);
			// else check if the previous position is the current index
			else if (i - 1 != index - 1 && tr_buffer[i].x != -1)
				graphics::drawLine(tr_buffer[i - 1].x, tr_buffer[i - 1].y, tr_buffer[i].x, tr_buffer[i].y, br);
		}
	}
}
