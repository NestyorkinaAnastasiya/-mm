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
		//��������� ������� � �������
		void CalculateLocals(int elementNumber);
		void CalculateA(int elementNumber, const Point* points, int kPoints);
		void CalculateAlphaMatrix(int elementNumber);
		void CalculateBettaMatrix(int elementNumber);
		void CalculateF(int elementNumber, const Point* points, int kPoints);
		/*
		//������ ������� ����� ��� ������� �������� 
		array<double, 2> g;
		//���������� ������ ����� ��� 1��� �������� �������
		void Calculate_g(int formNumber, int orientation, int elNumber);
		//���������� 1��� �������� ������� ��� ������ ����
		void CalculateBoundaries1ForNode(int node, double gi, double weight);
		//���� ������� �������� �������
		void CalculateBoundaries1(int number);*/
	public:
		SLAE();
		//����������� ������
		int n;
		//���������� �������
		Matrix A;
		//���������� ������ ������ �����
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