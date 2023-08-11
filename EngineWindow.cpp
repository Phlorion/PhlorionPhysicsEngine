#include "EngineWindow.h"

EngineWindow::EngineWindow()
{
	p = new GraphicsParticle();

	infoBox = new Info(p);
}

EngineWindow::~EngineWindow()
{
	delete infoBox;
	delete p;
}

void EngineWindow::update(float ms)
{
	infoBox->update(ms);

	/**
	* Debug for window to canvas units
	* 
	* graphics::getMouseState(mouse);
	* std::cout << graphics::windowToCanvasX(mouse.cur_pos_x) << " " << graphics::windowToCanvasY(mouse.cur_pos_y) << std::endl;
	*/

	if (p)
	{
		phlorion::Vector3 gravity = phlorion::Vector3();
		gravity.y = 10.f;
		p->addForce(gravity);
		p->update(ms);
	}
}

void EngineWindow::draw()
{
	// canvas bountries
	graphics::Brush br;
	br.fill_color[0] = br.fill_color[1] = br.fill_color[2] = 1.f;
	graphics::drawLine(0, 0, CANVAS_WIDTH, 0, br);
	graphics::drawLine(0.01, 0, 0.01, CANVAS_HEIGHT, br);
	graphics::drawLine(0, CANVAS_HEIGHT - 0.01, CANVAS_WIDTH, CANVAS_HEIGHT - 0.01, br);
	graphics::drawLine(CANVAS_WIDTH, 0, CANVAS_WIDTH, CANVAS_HEIGHT, br);

	infoBox->draw();
	p->draw();
}

void EngineWindow::init()
{
	graphics::setFont("assets\\big_noodle_titling.ttf");
}
