#pragma once
#include "slae.h"
#include "bmp24_lib.h"
#include <set>
using namespace slae;
namespace solver
{
	// ѕеречисление Ч это пользовательский тип, состо€щий из набора
	// целочисленных констант, называемых перечислител€ми.
	enum Derivative { dx, dy };
	class Solver
	{
		SLAE slae;

		vector<double> qSpline;

		//¬ектор нев€зки
		vector <double> r;
		//¬ектор спуска
		vector <double> z;

		//¬ычисление нормы вектора
		double Norm(const vector<double>& x);
		//—кал€рное произведение векторов
		double Scalar(const vector<double>& x, const vector<double>& y);
		double RelDiscrepancy();
		void LULOS();

		void linear(solver::Solver& s);

		//расчЄт решени€ и производных в точке
		double get_solution_in_point_bilinear(double x, double y, vector<double> qi);
		double get_solution_in_point_with_spline(double x, double y, vector<double> qi);
		double get_solution_in_point_bilinear(double x,	double y, int elementNumber, vector<double> qi);
		double get_solution_in_point_with_spline(double x, double y, int elementNumber, vector<double> qi);

		double get_derivative_in_point_bilinear(double x, double y,	vector<double> qi, Derivative d);
		double get_derivative_in_point_with_spline(double x, double y, vector<double> qi, Derivative d);
		double get_derivative_in_point_bilinear(double x, double y,	int elementNumber,	vector<double> qi,	Derivative d);
		double get_derivative_in_point_with_spline(double x, double y, int elementNumber, vector<double> qi, Derivative d);

	public:

		Solver();
		~Solver();

		void Solve();
		void GetSolutions();

		void Draw(unsigned int width, unsigned int height, unsigned int countOfIsolines);
		void Draw(unsigned int width, unsigned int height, unsigned int num_isolines, std::set<int> &nd);
	};
}