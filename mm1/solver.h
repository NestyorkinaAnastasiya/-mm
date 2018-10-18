#pragma once
#include "slae.h"
#include "bmp24_lib.h"
#include <set>
using namespace slae;
namespace solver
{
	// ������������ � ��� ���������������� ���, ��������� �� ������
	// ������������� ��������, ���������� ���������������.
	enum Derivative { dx, dy };
	class Solver
	{
		SLAE slae;

		vector<double> qSpline;

		//������ �������
		vector <double> r;
		//������ ������
		vector <double> z;

		//���������� ����� �������
		double Norm(const vector<double>& x);
		//��������� ������������ ��������
		double Scalar(const vector<double>& x, const vector<double>& y);
		double RelDiscrepancy();
		void LULOS();

		void linear(solver::Solver& s);

		//������ ������� � ����������� � �����
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