#pragma once
#include <iostream>
#include "graphics.h"
#include "config.h"
#include "PhlorionsPhysicsEngine.h"

/**
* Creates a Particle based on the particle of the physics engine.
* We modify the base particle so it can be displayed on the graphics library.
*/
class GraphicsParticle : public phlorion::Particle
{
private:
	// Radius of the particle
	float radius;
	// Time alive (ms)
	float time_alive;
	// Brush we will use to draw it on the screen
	graphics::Brush br;

	// Trajectory buffer
	const int tr_len = 200;
	phlorion::Vector3* tr_buffer = new phlorion::Vector3[tr_len];
	bool can_draw_trajectory = false;
	/**
	* This index is responsible for the new position to be added to the trajectory buffer when it is full
	* 
	* It starts from 0 meaning when the buffer gets filled for the first time the new position will replace the first position of the buffer etc.
	* This method gives the trajectory a nice effect, removing the oldest lines.
	*/
	int index = 0;
public:
	/**
	* Default Constructor
	*/
	GraphicsParticle();
	/**
	* Explicit Constructor
	* 
	*/
	GraphicsParticle(const phlorion::Vector3& pos, const phlorion::Vector3& vel, const phlorion::Vector3& acc, const phlorion::Vector3& force, 
		const float damp, const float mass, const float radius, const graphics::Brush br);
	// Deconstructor
	~GraphicsParticle();

	/**
	* Set the radius of the particle to be drawn.
	* 
	* @param The new radius of the particle.
	*/
	void setRadius(const float rad);
	
	/**
	* Get the radius of the particle.
	*/
	float getRadius() const;

	/**
	* Update function
	* Updates the data each frame
	* 
	* @param Time between each frame
	*/
	void update(float ms);

	/**
	* Draw function
	*/
	void draw();
};
