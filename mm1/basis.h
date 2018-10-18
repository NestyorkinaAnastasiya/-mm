#pragma once
#include "grid.h"

using namespace grd;
namespace basis
{
	// число билинейных функций
	const int nLFunc = 4;
	// число линейных функций
	const int nLFunc1D = 2;

	struct BilinearBasis
	{
		//указатели на функции вычисления базисных функций в точке
		array <function <double(double, double)>, nLFunc> phi;
		//указатели на функции вычисления d/dksi базисных функций в точке
		array <function <double(double, double)>, nLFunc> dphiksi;
		//указатели на функции вычисления d/detta базисных функций в точке
		array <function <double(double, double)>, nLFunc> dphietta;

		//конструктор(вычисление базисных функций)
		BilinearBasis();
	};

	const int nBQFunc = 9;
	const int nBQFunc1D = 3;

	struct BiquadraticBasis
	{
		//Указатели на функции вычисления базисных функций в точке
		array<function<double(double, double)>, nBQFunc> phi;
		//Указатели на функции вычисления d/dksi базисных функций в точке
		array<function<double(double, double)>, nBQFunc> dphiksi;
		//Указатели на функции вычисления d/detta базисных функций в точке
		array<function<double(double, double)>, nBQFunc> dphietta;
		BiquadraticBasis();
	};

	const int nSplineFunc = 16; //число базисных функций сплайна
	const int nSplineFunc1D = 4;

	class SplineBasis
	{
		std::pair<double, double> GetKsiEtta(double x, double y, int element_number, double &hx, double &hy);
		std::pair<double, double> GetKsiEtta(double x, double y, int element_number);

	public:
		//функции вычисления бф сплайна в точке мастер-координат, их производных, первых и вторых
		array <function<double(double, double, double, double)>, nSplineFunc> phi;
		array <function<double(double, double, double, double)>, nSplineFunc> dphiksi;
		array <function<double(double, double, double, double)>, nSplineFunc> dphietta;
		array <function<double(double, double, double, double)>, nSplineFunc> d2phiksi;
		array <function<double(double, double, double, double)>, nSplineFunc> d2phietta;
		SplineBasis();

		//функции вычисления бф сплайна в точке (x,y), их производных
		double phi_i(int i, double x, double y, int element_number);
		double dphi_ix(int i, double x, double y, int element_number);
		double dphi_iy(int i, double x, double y, int element_number);

	};
}