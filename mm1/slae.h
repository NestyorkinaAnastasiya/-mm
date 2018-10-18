#pragma once
#include "integration.h"
#include "parameters.h"
using namespace matrix;
using namespace basis;
using namespace integration;
using namespace parameters;

namespace slae
{
	class SLAE : private GaussIntegration
	{
		//локальные матрицы и векторы
		void CalculateLocals(int elementNumber);
		void CalculateA(int elementNumber, const Point* points, int kPoints);
		void CalculateAlphaMatrix(int elementNumber);
		void CalculateBettaMatrix(int elementNumber);
		void CalculateF(int elementNumber, const Point* points, int kPoints);
		/*
		//Вектор праввой части для первого краевого 
		array<double, 2> g;
		//Нахождение правой части для 1ого краевого условия
		void Calculate_g(int formNumber, int orientation, int elNumber);
		//Вычисление 1ого краевого условия для одного узла
		void CalculateBoundaries1ForNode(int node, double gi, double weight);
		//Учёт первого краевого условия
		void CalculateBoundaries1(int number);*/
	public:
		SLAE();
		//Размерность задачи
		int n;
		//Глобальная матрица
		Matrix A;
		//Глобальный вектор правой части
		vector <double> F;
		vector <double> q;
		SplineBasis splBasis;
		BilinearBasis blBasis;
		void GenerateSLAE();
		double GetSolutionInThePoint(double x, double y, int elementNumber);
		double GetSolutionInThePointWithSpline(double x, double y, int elementNumber);
		~SLAE() {};
	};
}