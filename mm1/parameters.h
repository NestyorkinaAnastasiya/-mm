#pragma once
#pragma once
#include <string>
namespace parameters
{
	bool useLU = false;
	int solver = 0;
	// ѕараметры регул€ризации
	double alpha = 0;
	double betta = 0;

	/* ѕараметры решател€ */
	int maxCountOfIterations;
	double epsilon;
	int gmresM;

	struct Parameters
	{
		Parameters();
		~Parameters();
	};
}