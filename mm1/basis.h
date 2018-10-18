#pragma once
#include "grid.h"

using namespace grd;
namespace basis
{
	// ����� ���������� �������
	const int nLFunc = 4;
	// ����� �������� �������
	const int nLFunc1D = 2;

	struct BilinearBasis
	{
		//��������� �� ������� ���������� �������� ������� � �����
		array <function <double(double, double)>, nLFunc> phi;
		//��������� �� ������� ���������� d/dksi �������� ������� � �����
		array <function <double(double, double)>, nLFunc> dphiksi;
		//��������� �� ������� ���������� d/detta �������� ������� � �����
		array <function <double(double, double)>, nLFunc> dphietta;

		//�����������(���������� �������� �������)
		BilinearBasis();
	};

	const int nBQFunc = 9;
	const int nBQFunc1D = 3;

	struct BiquadraticBasis
	{
		//��������� �� ������� ���������� �������� ������� � �����
		array<function<double(double, double)>, nBQFunc> phi;
		//��������� �� ������� ���������� d/dksi �������� ������� � �����
		array<function<double(double, double)>, nBQFunc> dphiksi;
		//��������� �� ������� ���������� d/detta �������� ������� � �����
		array<function<double(double, double)>, nBQFunc> dphietta;
		BiquadraticBasis();
	};

	const int nSplineFunc = 16; //����� �������� ������� �������
	const int nSplineFunc1D = 4;

	class SplineBasis
	{
		std::pair<double, double> GetKsiEtta(double x, double y, int element_number, double &hx, double &hy);
		std::pair<double, double> GetKsiEtta(double x, double y, int element_number);

	public:
		//������� ���������� �� ������� � ����� ������-���������, �� �����������, ������ � ������
		array <function<double(double, double, double, double)>, nSplineFunc> phi;
		array <function<double(double, double, double, double)>, nSplineFunc> dphiksi;
		array <function<double(double, double, double, double)>, nSplineFunc> dphietta;
		array <function<double(double, double, double, double)>, nSplineFunc> d2phiksi;
		array <function<double(double, double, double, double)>, nSplineFunc> d2phietta;
		SplineBasis();

		//������� ���������� �� ������� � ����� (x,y), �� �����������
		double phi_i(int i, double x, double y, int element_number);
		double dphi_ix(int i, double x, double y, int element_number);
		double dphi_iy(int i, double x, double y, int element_number);

	};
}