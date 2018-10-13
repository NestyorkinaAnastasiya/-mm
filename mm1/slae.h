#include "integration.h"
using namespace matrix;
using namespace basis;
using namespace integration;

namespace slae
{
	class SLAE : private GaussIntegration
	{
		//����������� ������
		int n;
		//������������ ���������� �������� � ��������
		int maxiter = 10000;
		//�������� ������� ����
		const double eps = 1e-10;
		BilinearBasis basis4;
		BiquadraticBasis basis9;
		//���������� �������
		Matrix A;
		//��������� �������
		//������� ��������
		array<array<double, 4>, 4> G;
		//������� �����
		array<array<double, 4>, 4> M;
		//��������� ������ ������ �����
		array <double, 4> locF;
		//���������� ������ ������ �����
		vector <double> F;
		//������ ������������� �������
		vector <double> u;

		double dphix(int i, int elementNumber, double ksi, double etta);
		double dphiy(int i, int elementNumber, double ksi, double etta);
		double GetAbsJ(int elementNumber, double ksi, double etta);
		double GetAbsBQJ(int elementNumber, double ksi, double etta);

		//������ ��������� ������ ��������
		void CalculateG(int elementNumber);
		//������ ��������� ������ ����
		void CalculateM(int elementNumber);
		//������ ��������� ������ ������
		void CalculateLocalF(int elementNumber);
		//������ ��������� ������(��������) � ���������� � ����������
		void CalculateLocals(int elementNumber);

		//������ ������� ����� ��� ������� �������� 
		array<double, 2> g;
		//���������� ������ ����� ��� 1��� �������� �������
		void Calculate_g(int formNumber, int orientation, int elNumber);
		//���������� 1��� �������� ������� ��� ������ ����
		void CalculateBoundaries1ForNode(int node, double gi, double weight);
		//���� ������� �������� �������
		void CalculateBoundaries1(int number);

		
		//������ �������
		vector <double> r;
		//������ ������
		vector <double> z;

		//���������� ����� �������
		double Norm(const vector<double>& x);
		//��������� ������������ ��������
		double Scalar(const vector<double>& x, const vector<double>& y);

		//��������� ���� �� i-�� �������� �� �������
		void GenerateSLAE();
		
		double Rel_Discrepancy();
		//�������� ��� � LU-�������������
		void LULOS();
	public:
		SLAE();
		void Solve();
		~SLAE() {};
	};
}