#include "UI.h"

Info::Info(GraphicsParticle* particle)
{
	x = 10.f;
	y = 20.f;
	size = 0.15;

	p = particle;
}

Info::~Info()
{

}

void Info::update(float ms)
{
}

void Info::draw()
{
	graphics::Brush br;
	br.fill_color[0] = br.fill_color[1] = br.fill_color[2] = 1.f;
	graphics::drawText(graphics::windowToCanvasX(x), graphics::windowToCanvasY(y), size, "Position[x, y]: [" + std::to_string(p->getPosition().x) + ", " + std::to_string(p->getPosition().y) + "]", br);
	graphics::drawText(graphics::windowToCanvasX(x), graphics::windowToCanvasY(y + 20), size, "Velocity[x, y]: [" + std::to_string(p->getVelocity().x) + ", " + std::to_string(p->getVelocity().y) + "]", br);
	graphics::drawText(graphics::windowToCanvasX(x), graphics::windowToCanvasY(y + 40), size, "Acceleration[x, y]: [" + std::to_string(p->getLastAcceleration().x) + ", " + std::to_string(p->getLastAcceleration().y) + "]", br);
}
