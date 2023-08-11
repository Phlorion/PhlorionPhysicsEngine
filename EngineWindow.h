#pragma once
#include <iostream>
#include "graphics.h"
#include "config.h"
#include "UI.h"
#include "Particle.h"


class EngineWindow
{
private:
	Info* infoBox;
	graphics::MouseState mouse;

	GraphicsParticle* p;

public:
	EngineWindow();
	~EngineWindow();

	void update(float ms);
	void draw();
	void init();
};