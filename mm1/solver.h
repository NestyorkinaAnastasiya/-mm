#pragma once
#include "integration.h"
#include "slae.h"
#include "parameters.h"
using namespace integration;
using namespace basis;
using namespace slae;
using namespace parameters;
namespace solver
{
	class Solver : private GaussIntegration
	{
		SLAE slae;

		vector<double> qSpline;
		vector<double> qSolution;
		BilinearBasis blBasis;
		SplineBasis splBasis;
		//��������� ������� � �������
		void CalculateLocals(int elementNumber);
		void CalculateA(int elementNumber, const Point* points, int kPoints);
		void CalculateAlphaMatrix(int elementNumber);
		void CalculateBettaMatrix(int elementNumber);
		void CalculateF(int elementNumber, const Point* points, int kPoints);


		void linear(solver::Solver& s);

		//������ ������� � ����������� � �����
		double get_solution_in_point_bilinear(double x, double y, myvector::MyVector qi);
		double get_solution_in_point_with_spline(double x, double y, myvector::MyVector qi);
		double get_solution_in_point_bilinear(double x,
			double y,
			int element_number,
			myvector::MyVector qi);
		double get_solution_in_point_with_spline(double x,
			double y,
			int element_number,
			myvector::MyVector qi);

		double get_derivative_in_point_bilinear(double x,
			double y,
			myvector::MyVector qi,
			Derivative d_var);
		double get_derivative_in_point_with_spline(double x,
			double y,
			myvector::MyVector qi,
			Derivative d_var);
		double get_derivative_in_point_bilinear(double x,
			double y,
			int element_number,
			myvector::MyVector qi,
			Derivative d_var);
		double get_derivative_in_point_with_spline(double x,
			double y,
			int element_number,
			myvector::MyVector qi,
			Derivative d_var);

		//��������� ������� � �������
		void calculate_locals(int element_number);
		void calculate_A(int element_number, const point::Point* points, int k_points);
		void calculate_alpha_matrix(int element_number);
		void calculate_betta_matrix(int element_number);
		void calculate_F(int element_number, const point::Point* points, int k_points);

	public:

		Solver();
		~Solver();

		void Solve();
		void get_solutions(std::ofstream& solution_f_out,
			std::ofstream& info_f_out,
			std::ifstream& points);

		void draw(unsigned int width, unsigned int height, unsigned int num_isolines);
		void draw(unsigned int width, unsigned int height, unsigned int num_isolines, std::set<int> &nd);
	};
}