#pragma once
#pragma once
#include <string>
namespace parameters
{
	bool useLU = false;
	int solver = 0;
	// ��������� �������������
	double alpha = 0;
	double betta = 0;

	/* ��������� �������� */
	int maxCountOfIterations;
	double epsilon;
	int gmresM;

	struct Parameters
	{
		Parameters();
		~Parameters();
	};
}