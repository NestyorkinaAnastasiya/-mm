#pragma once
#include "integration.h"

namespace tests
{
	struct Tests
	{
		int test;
		double Lambda();
		double Sigma();
		double Ug(double x, double y);
		double Fi(double x, double y);
		Tests();
	};
}