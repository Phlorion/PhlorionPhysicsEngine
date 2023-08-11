#pragma once
#include <iostream>
#include "config.h"
#include "graphics.h"
#include "Particle.h"


/**
* Info Box 
* Provides info for a particle (position, velocity, acceleration, etc.)
*/
class Info
{
private:
	float x;
	float y;
	float size;

	GraphicsParticle* p;

public:
	Info(GraphicsParticle* particle);
	~Info();

	void update(float ms);
	void draw();
};