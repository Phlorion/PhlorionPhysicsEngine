#include "config.h"
#include "EngineWindow.h"
#include "graphics.h"

void update(float ms)
{
	ms = graphics::getDeltaTime();
	EngineWindow* engine_window = reinterpret_cast<EngineWindow*>(graphics::getUserData());
	engine_window->update(ms);
}

void draw()
{
	EngineWindow* engine_window = reinterpret_cast<EngineWindow*>(graphics::getUserData());
	engine_window->draw();
}

int main()
{
	EngineWindow engine_window;

	graphics::createWindow(SCR_WIDTH, SCR_HEIGHT, "Physics Engine");

	graphics::setUserData(&engine_window);

	graphics::setDrawFunction(draw);
	graphics::setUpdateFunction(update);

	graphics::setCanvasSize(CANVAS_WIDTH, CANVAS_HEIGHT);
	graphics::setCanvasScaleMode(graphics::CANVAS_SCALE_FIT);

	engine_window.init();

	graphics::startMessageLoop();

	return 0;
}