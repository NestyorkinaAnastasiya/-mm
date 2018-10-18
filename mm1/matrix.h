#pragma once
#pragma once
#include <algorithm>
#include "basis.h"

namespace matrix
{
	class Matrix
	{
		// ������� 
		// ����������� �������
		int n;
		// ������ �������� ggl, ggu, jg
		int size;

		// ��������� ������ �����(��������)
		vector <int> ig;
		// ������ �������� ��������������� ��������� 
		vector <int> jg;
		// ������������ �������� �������
		vector <double> di;
		// ��������������� �������� ������� ������������ �������
		vector <double> ggl;
		// ��������������� �������� �������� ������������ �������
		vector <double> ggu;

		//���������� ������� � �������������
		vector <double> L;
		vector <double> D;
		vector <double> U;

	public:
		// ��������� �������� �������
		void CreatePortret();
		// ������������� ������� ����� ��������� ��������
		void Initialize(int size1, int size2);
		// ��������� ������� �� ������
		vector <double> operator*(vector<double> x);

		void Clear();

		// ��������� ����������������� ������� �� ������
		void MultiplyATx(vector<double> a, vector<double>& result);
		// �������� ������� � �������
		void AddElement(int i, int j, double element);
		// �������� ������� � �������
		void ChangeElement(int i, int j, double element);
		// �������� ��������� � ������ (��� ��������)
		void Sum(vector<double>& result);
		//LU-������������
		void LU();
		//��������������� ������� ��� ��������
		void LYF(const vector<double>& C, vector<double>& yl);
		void UXY(const vector<double>& C, vector<double>& yu);
	};
}