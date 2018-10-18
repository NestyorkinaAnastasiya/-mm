#include "basis.h"

namespace basis
{
	BilinearBasis::BilinearBasis()
	{
		//указатели на функции вычисления одномерных базисных функций в точке
		array <function <double(double)>, nLFunc1D> _phi_;
		//указатели на функции вычисления d/dksi одномерных базисных функций в точке
		array <function <double(double)>, nLFunc1D> dphi_ksi;

		_phi_[0] = [](double ksi) { return 1 - ksi; };
		_phi_[1] = [](double ksi) { return ksi; };
		dphi_ksi[0] = [](double ksi) { return -1; };
		dphi_ksi[1] = [](double ksi) { return  1; };

		//указатели на функции вычисления базисных функций в точке
		phi[0] = [_phi_](double ksi, double etta) { return _phi_[0](ksi) * _phi_[0](etta); };
		phi[1] = [_phi_](double ksi, double etta) { return _phi_[1](ksi) * _phi_[0](etta); };
		phi[2] = [_phi_](double ksi, double etta) { return _phi_[0](ksi) * _phi_[1](etta); };
		phi[3] = [_phi_](double ksi, double etta) { return _phi_[1](ksi) * _phi_[1](etta); };
		//указатели на функции вычисления d/dksi базисных функций в точке
		dphiksi[0] = [_phi_, dphi_ksi](double ksi, double etta) { return dphi_ksi[0](ksi) * _phi_[0](etta); };
		dphiksi[1] = [_phi_, dphi_ksi](double ksi, double etta) { return dphi_ksi[1](ksi) * _phi_[0](etta); };
		dphiksi[2] = [_phi_, dphi_ksi](double ksi, double etta) { return dphi_ksi[0](ksi) * _phi_[1](etta); };
		dphiksi[3] = [_phi_, dphi_ksi](double ksi, double etta) { return dphi_ksi[1](ksi) * _phi_[1](etta); };
		//указатели на функции вычисления d/detta базисных функций в точке
		dphietta[0] = [_phi_, dphi_ksi](double ksi, double etta) { return _phi_[0](ksi) * dphi_ksi[0](etta); };
		dphietta[1] = [_phi_, dphi_ksi](double ksi, double etta) { return _phi_[1](ksi) * dphi_ksi[0](etta); };
		dphietta[2] = [_phi_, dphi_ksi](double ksi, double etta) { return _phi_[0](ksi) * dphi_ksi[1](etta); };
		dphietta[3] = [_phi_, dphi_ksi](double ksi, double etta) { return _phi_[1](ksi) * dphi_ksi[1](etta); };
	}

	BiquadraticBasis::BiquadraticBasis()
	{
		array <function<double(double)>, nBQFunc1D> phi_;
		phi_[0] = [](double ksi) { return 2 * (ksi - 0.5) * (ksi - 1); };
		phi_[1] = [](double ksi) { return -4 * ksi * (ksi - 1); };
		phi_[2] = [](double ksi) { return 2 * ksi * (ksi - 0.5); };

		array <function<double(double)>, nBQFunc1D> dphi_ksi;
		dphi_ksi[0] = [](double ksi) { return 4 * ksi - 3; };
		dphi_ksi[1] = [](double ksi) { return  -8 * ksi + 4; };
		dphi_ksi[2] = [](double ksi) { return 4 * ksi - 1; };

		phi[0] = [phi_](double ksi, double etta) { return phi_[0](ksi) * phi_[0](etta); };
		phi[1] = [phi_](double ksi, double etta) { return phi_[1](ksi) * phi_[0](etta); };
		phi[2] = [phi_](double ksi, double etta) { return phi_[2](ksi) * phi_[0](etta); };
		phi[3] = [phi_](double ksi, double etta) { return phi_[0](ksi) * phi_[1](etta); };
		phi[4] = [phi_](double ksi, double etta) { return phi_[1](ksi) * phi_[1](etta); };
		phi[5] = [phi_](double ksi, double etta) { return phi_[2](ksi) * phi_[1](etta); };
		phi[6] = [phi_](double ksi, double etta) { return phi_[0](ksi) * phi_[2](etta); };
		phi[7] = [phi_](double ksi, double etta) { return phi_[1](ksi) * phi_[2](etta); };
		phi[8] = [phi_](double ksi, double etta) { return phi_[2](ksi) * phi_[2](etta); };

		dphiksi[0] = [phi_, dphi_ksi](double ksi, double etta) { return dphi_ksi[0](ksi) * phi_[0](etta); };
		dphiksi[1] = [phi_, dphi_ksi](double ksi, double etta) { return dphi_ksi[1](ksi) * phi_[0](etta); };
		dphiksi[2] = [phi_, dphi_ksi](double ksi, double etta) { return dphi_ksi[2](ksi) * phi_[0](etta); };
		dphiksi[3] = [phi_, dphi_ksi](double ksi, double etta) { return dphi_ksi[0](ksi) * phi_[1](etta); };
		dphiksi[4] = [phi_, dphi_ksi](double ksi, double etta) { return dphi_ksi[1](ksi) * phi_[1](etta); };
		dphiksi[5] = [phi_, dphi_ksi](double ksi, double etta) { return dphi_ksi[2](ksi) * phi_[1](etta); };
		dphiksi[6] = [phi_, dphi_ksi](double ksi, double etta) { return dphi_ksi[0](ksi) * phi_[2](etta); };
		dphiksi[7] = [phi_, dphi_ksi](double ksi, double etta) { return dphi_ksi[1](ksi) * phi_[2](etta); };
		dphiksi[8] = [phi_, dphi_ksi](double ksi, double etta) { return dphi_ksi[2](ksi) * phi_[2](etta); };

		dphietta[0] = [phi_, dphi_ksi](double ksi, double etta) { return phi_[0](ksi) * dphi_ksi[0](etta); };
		dphietta[1] = [phi_, dphi_ksi](double ksi, double etta) { return phi_[1](ksi) * dphi_ksi[0](etta); };
		dphietta[2] = [phi_, dphi_ksi](double ksi, double etta) { return phi_[2](ksi) * dphi_ksi[0](etta); };
		dphietta[3] = [phi_, dphi_ksi](double ksi, double etta) { return phi_[0](ksi) * dphi_ksi[1](etta); };
		dphietta[4] = [phi_, dphi_ksi](double ksi, double etta) { return phi_[1](ksi) * dphi_ksi[1](etta); };
		dphietta[5] = [phi_, dphi_ksi](double ksi, double etta) { return phi_[2](ksi) * dphi_ksi[1](etta); };
		dphietta[6] = [phi_, dphi_ksi](double ksi, double etta) { return phi_[0](ksi) * dphi_ksi[2](etta); };
		dphietta[7] = [phi_, dphi_ksi](double ksi, double etta) { return phi_[1](ksi) * dphi_ksi[2](etta); };
		dphietta[8] = [phi_, dphi_ksi](double ksi, double etta) { return phi_[2](ksi) * dphi_ksi[2](etta); };
	}

	SplineBasis::SplineBasis()
	{
		array <function<double(double, double)>, nSplineFunc1D> phi_;
		phi_[0] = [](double ksi, double h) { return 1 - 3 * ksi*ksi + 2 * ksi*ksi*ksi; };
		phi_[1] = [](double ksi, double h) { return h * (ksi - 2 * ksi*ksi + ksi * ksi*ksi); };
		phi_[2] = [](double ksi, double h) { return 3 * ksi*ksi - 2 * ksi*ksi; };
		phi_[3] = [](double ksi, double h) { return h * (ksi*ksi*ksi - ksi * ksi); };

		array <function<double(double, double)>, nSplineFunc1D> dphi_ksi;
		dphi_ksi[0] = [](double ksi, double h) { return 6 * (ksi*ksi - ksi); };
		dphi_ksi[1] = [](double ksi, double h) { return h * (1 - 4 * ksi + 3 * ksi*ksi); };
		dphi_ksi[2] = [](double ksi, double h) { return 6 * (ksi - ksi * ksi); };
		dphi_ksi[3] = [](double ksi, double h) { return h * (3 * ksi*ksi - 2 * ksi); };

		array <function<double(double, double)>, nSplineFunc1D> d2phi_ksi;
		d2phi_ksi[0] = [](double ksi, double h) { return 12 * ksi - 6; };
		d2phi_ksi[1] = [](double ksi, double h) { return h * (6 * ksi - 4); };
		d2phi_ksi[2] = [](double ksi, double h) { return 6 - 12 * ksi; };
		d2phi_ksi[3] = [](double ksi, double h) { return h * (6 * ksi - 2); };


		phi[0] = [phi_](double ksi, double etta, double hx, double hy) { return phi_[0](ksi, hx) * phi_[0](etta, hy); };
		phi[1] = [phi_](double ksi, double etta, double hx, double hy) { return phi_[1](ksi, hx) * phi_[0](etta, hy); };
		phi[2] = [phi_](double ksi, double etta, double hx, double hy) { return phi_[0](ksi, hx) * phi_[1](etta, hy); };
		phi[3] = [phi_](double ksi, double etta, double hx, double hy) { return phi_[1](ksi, hx) * phi_[1](etta, hy); };
		phi[4] = [phi_](double ksi, double etta, double hx, double hy) { return phi_[2](ksi, hx) * phi_[0](etta, hy); };
		phi[5] = [phi_](double ksi, double etta, double hx, double hy) { return phi_[3](ksi, hx) * phi_[0](etta, hy); };
		phi[6] = [phi_](double ksi, double etta, double hx, double hy) { return phi_[2](ksi, hx) * phi_[1](etta, hy); };
		phi[7] = [phi_](double ksi, double etta, double hx, double hy) { return phi_[3](ksi, hx) * phi_[1](etta, hy); };
		phi[8] = [phi_](double ksi, double etta, double hx, double hy) { return phi_[0](ksi, hx) * phi_[2](etta, hy); };
		phi[9] = [phi_](double ksi, double etta, double hx, double hy) { return phi_[1](ksi, hx) * phi_[2](etta, hy); };
		phi[10] = [phi_](double ksi, double etta, double hx, double hy) { return phi_[0](ksi, hx) * phi_[3](etta, hy); };
		phi[11] = [phi_](double ksi, double etta, double hx, double hy) { return phi_[1](ksi, hx) * phi_[3](etta, hy); };
		phi[12] = [phi_](double ksi, double etta, double hx, double hy) { return phi_[2](ksi, hx) * phi_[2](etta, hy); };
		phi[13] = [phi_](double ksi, double etta, double hx, double hy) { return phi_[3](ksi, hx) * phi_[2](etta, hy); };
		phi[14] = [phi_](double ksi, double etta, double hx, double hy) { return phi_[2](ksi, hx) * phi_[3](etta, hy); };
		phi[15] = [phi_](double ksi, double etta, double hx, double hy) { return phi_[3](ksi, hx) * phi_[3](etta, hy); };

		dphiksi[0] = [phi_, dphi_ksi](double ksi, double etta, double hx, double hy) { return dphi_ksi[0](ksi, hx) * phi_[0](etta, hy); };
		dphiksi[1] = [phi_, dphi_ksi](double ksi, double etta, double hx, double hy) { return dphi_ksi[1](ksi, hx) * phi_[0](etta, hy); };
		dphiksi[2] = [phi_, dphi_ksi](double ksi, double etta, double hx, double hy) { return dphi_ksi[0](ksi, hx) * phi_[1](etta, hy); };
		dphiksi[3] = [phi_, dphi_ksi](double ksi, double etta, double hx, double hy) { return dphi_ksi[1](ksi, hx) * phi_[1](etta, hy); };
		dphiksi[4] = [phi_, dphi_ksi](double ksi, double etta, double hx, double hy) { return dphi_ksi[2](ksi, hx) * phi_[0](etta, hy); };
		dphiksi[5] = [phi_, dphi_ksi](double ksi, double etta, double hx, double hy) { return dphi_ksi[3](ksi, hx) * phi_[0](etta, hy); };
		dphiksi[6] = [phi_, dphi_ksi](double ksi, double etta, double hx, double hy) { return dphi_ksi[2](ksi, hx) * phi_[1](etta, hy); };
		dphiksi[7] = [phi_, dphi_ksi](double ksi, double etta, double hx, double hy) { return dphi_ksi[3](ksi, hx) * phi_[1](etta, hy); };
		dphiksi[8] = [phi_, dphi_ksi](double ksi, double etta, double hx, double hy) { return dphi_ksi[0](ksi, hx) * phi_[2](etta, hy); };
		dphiksi[9] = [phi_, dphi_ksi](double ksi, double etta, double hx, double hy) { return dphi_ksi[1](ksi, hx) * phi_[2](etta, hy); };
		dphiksi[10] = [phi_, dphi_ksi](double ksi, double etta, double hx, double hy) { return dphi_ksi[0](ksi, hx) * phi_[3](etta, hy); };
		dphiksi[11] = [phi_, dphi_ksi](double ksi, double etta, double hx, double hy) { return dphi_ksi[1](ksi, hx) * phi_[3](etta, hy); };
		dphiksi[12] = [phi_, dphi_ksi](double ksi, double etta, double hx, double hy) { return dphi_ksi[2](ksi, hx) * phi_[2](etta, hy); };
		dphiksi[13] = [phi_, dphi_ksi](double ksi, double etta, double hx, double hy) { return dphi_ksi[3](ksi, hx) * phi_[2](etta, hy); };
		dphiksi[14] = [phi_, dphi_ksi](double ksi, double etta, double hx, double hy) { return dphi_ksi[2](ksi, hx) * phi_[3](etta, hy); };
		dphiksi[15] = [phi_, dphi_ksi](double ksi, double etta, double hx, double hy) { return dphi_ksi[3](ksi, hx) * phi_[3](etta, hy); };

		dphietta[0] = [phi_, dphi_ksi](double ksi, double etta, double hx, double hy) { return phi_[0](ksi, hx) * dphi_ksi[0](etta, hy); };
		dphietta[1] = [phi_, dphi_ksi](double ksi, double etta, double hx, double hy) { return phi_[1](ksi, hx) * dphi_ksi[0](etta, hy); };
		dphietta[2] = [phi_, dphi_ksi](double ksi, double etta, double hx, double hy) { return phi_[0](ksi, hx) * dphi_ksi[1](etta, hy); };
		dphietta[3] = [phi_, dphi_ksi](double ksi, double etta, double hx, double hy) { return phi_[1](ksi, hx) * dphi_ksi[1](etta, hy); };
		dphietta[4] = [phi_, dphi_ksi](double ksi, double etta, double hx, double hy) { return phi_[2](ksi, hx) * dphi_ksi[0](etta, hy); };
		dphietta[5] = [phi_, dphi_ksi](double ksi, double etta, double hx, double hy) { return phi_[3](ksi, hx) * dphi_ksi[0](etta, hy); };
		dphietta[6] = [phi_, dphi_ksi](double ksi, double etta, double hx, double hy) { return phi_[2](ksi, hx) * dphi_ksi[1](etta, hy); };
		dphietta[7] = [phi_, dphi_ksi](double ksi, double etta, double hx, double hy) { return phi_[3](ksi, hx) * dphi_ksi[1](etta, hy); };
		dphietta[8] = [phi_, dphi_ksi](double ksi, double etta, double hx, double hy) { return phi_[0](ksi, hx) * dphi_ksi[2](etta, hy); };
		dphietta[9] = [phi_, dphi_ksi](double ksi, double etta, double hx, double hy) { return phi_[1](ksi, hx) * dphi_ksi[2](etta, hy); };
		dphietta[10] = [phi_, dphi_ksi](double ksi, double etta, double hx, double hy) { return phi_[0](ksi, hx) * dphi_ksi[3](etta, hy); };
		dphietta[11] = [phi_, dphi_ksi](double ksi, double etta, double hx, double hy) { return phi_[1](ksi, hx) * dphi_ksi[3](etta, hy); };
		dphietta[12] = [phi_, dphi_ksi](double ksi, double etta, double hx, double hy) { return phi_[2](ksi, hx) * dphi_ksi[2](etta, hy); };
		dphietta[13] = [phi_, dphi_ksi](double ksi, double etta, double hx, double hy) { return phi_[3](ksi, hx) * dphi_ksi[2](etta, hy); };
		dphietta[14] = [phi_, dphi_ksi](double ksi, double etta, double hx, double hy) { return phi_[2](ksi, hx) * dphi_ksi[3](etta, hy); };
		dphietta[15] = [phi_, dphi_ksi](double ksi, double etta, double hx, double hy) { return phi_[3](ksi, hx) * dphi_ksi[3](etta, hy); };

		d2phiksi[0] = [phi_, d2phi_ksi](double ksi, double etta, double hx, double hy) { return d2phi_ksi[0](ksi, hx) * phi_[0](etta, hy); };
		d2phiksi[1] = [phi_, d2phi_ksi](double ksi, double etta, double hx, double hy) { return d2phi_ksi[1](ksi, hx) * phi_[0](etta, hy); };
		d2phiksi[2] = [phi_, d2phi_ksi](double ksi, double etta, double hx, double hy) { return d2phi_ksi[0](ksi, hx) * phi_[1](etta, hy); };
		d2phiksi[3] = [phi_, d2phi_ksi](double ksi, double etta, double hx, double hy) { return d2phi_ksi[1](ksi, hx) * phi_[1](etta, hy); };
		d2phiksi[4] = [phi_, d2phi_ksi](double ksi, double etta, double hx, double hy) { return d2phi_ksi[2](ksi, hx) * phi_[0](etta, hy); };
		d2phiksi[5] = [phi_, d2phi_ksi](double ksi, double etta, double hx, double hy) { return d2phi_ksi[3](ksi, hx) * phi_[0](etta, hy); };
		d2phiksi[6] = [phi_, d2phi_ksi](double ksi, double etta, double hx, double hy) { return d2phi_ksi[2](ksi, hx) * phi_[1](etta, hy); };
		d2phiksi[7] = [phi_, d2phi_ksi](double ksi, double etta, double hx, double hy) { return d2phi_ksi[3](ksi, hx) * phi_[1](etta, hy); };
		d2phiksi[8] = [phi_, d2phi_ksi](double ksi, double etta, double hx, double hy) { return d2phi_ksi[0](ksi, hx) * phi_[2](etta, hy); };
		d2phiksi[9] = [phi_, d2phi_ksi](double ksi, double etta, double hx, double hy) { return d2phi_ksi[1](ksi, hx) * phi_[2](etta, hy); };
		d2phiksi[10] = [phi_, d2phi_ksi](double ksi, double etta, double hx, double hy) { return d2phi_ksi[0](ksi, hx) * phi_[3](etta, hy); };
		d2phiksi[11] = [phi_, d2phi_ksi](double ksi, double etta, double hx, double hy) { return d2phi_ksi[1](ksi, hx) * phi_[3](etta, hy); };
		d2phiksi[12] = [phi_, d2phi_ksi](double ksi, double etta, double hx, double hy) { return d2phi_ksi[2](ksi, hx) * phi_[2](etta, hy); };
		d2phiksi[13] = [phi_, d2phi_ksi](double ksi, double etta, double hx, double hy) { return d2phi_ksi[3](ksi, hx) * phi_[2](etta, hy); };
		d2phiksi[14] = [phi_, d2phi_ksi](double ksi, double etta, double hx, double hy) { return d2phi_ksi[2](ksi, hx) * phi_[3](etta, hy); };
		d2phiksi[15] = [phi_, d2phi_ksi](double ksi, double etta, double hx, double hy) { return d2phi_ksi[3](ksi, hx) * phi_[3](etta, hy); };

		d2phietta[0] = [phi_, d2phi_ksi](double ksi, double etta, double hx, double hy) { return phi_[0](ksi, hx) * d2phi_ksi[0](etta, hy); };
		d2phietta[1] = [phi_, d2phi_ksi](double ksi, double etta, double hx, double hy) { return phi_[1](ksi, hx) * d2phi_ksi[0](etta, hy); };
		d2phietta[2] = [phi_, d2phi_ksi](double ksi, double etta, double hx, double hy) { return phi_[0](ksi, hx) * d2phi_ksi[1](etta, hy); };
		d2phietta[3] = [phi_, d2phi_ksi](double ksi, double etta, double hx, double hy) { return phi_[1](ksi, hx) * d2phi_ksi[1](etta, hy); };
		d2phietta[4] = [phi_, d2phi_ksi](double ksi, double etta, double hx, double hy) { return phi_[2](ksi, hx) * d2phi_ksi[0](etta, hy); };
		d2phietta[5] = [phi_, d2phi_ksi](double ksi, double etta, double hx, double hy) { return phi_[3](ksi, hx) * d2phi_ksi[0](etta, hy); };
		d2phietta[6] = [phi_, d2phi_ksi](double ksi, double etta, double hx, double hy) { return phi_[2](ksi, hx) * d2phi_ksi[1](etta, hy); };
		d2phietta[7] = [phi_, d2phi_ksi](double ksi, double etta, double hx, double hy) { return phi_[3](ksi, hx) * d2phi_ksi[1](etta, hy); };
		d2phietta[8] = [phi_, d2phi_ksi](double ksi, double etta, double hx, double hy) { return phi_[0](ksi, hx) * d2phi_ksi[2](etta, hy); };
		d2phietta[9] = [phi_, d2phi_ksi](double ksi, double etta, double hx, double hy) { return phi_[1](ksi, hx) * d2phi_ksi[2](etta, hy); };
		d2phietta[10] = [phi_, d2phi_ksi](double ksi, double etta, double hx, double hy) { return phi_[0](ksi, hx) * d2phi_ksi[3](etta, hy); };
		d2phietta[11] = [phi_, d2phi_ksi](double ksi, double etta, double hx, double hy) { return phi_[1](ksi, hx) * d2phi_ksi[3](etta, hy); };
		d2phietta[12] = [phi_, d2phi_ksi](double ksi, double etta, double hx, double hy) { return phi_[2](ksi, hx) * d2phi_ksi[2](etta, hy); };
		d2phietta[13] = [phi_, d2phi_ksi](double ksi, double etta, double hx, double hy) { return phi_[3](ksi, hx) * d2phi_ksi[2](etta, hy); };
		d2phietta[14] = [phi_, d2phi_ksi](double ksi, double etta, double hx, double hy) { return phi_[2](ksi, hx) * d2phi_ksi[3](etta, hy); };
		d2phietta[15] = [phi_, d2phi_ksi](double ksi, double etta, double hx, double hy) { return phi_[3](ksi, hx) * d2phi_ksi[3](etta, hy); };

	}

	pair<double, double> SplineBasis::GetKsiEtta(double x, double y, int element_number)
	{
		double x_left = grid.nodes[grid.elements[element_number].nodes[0]].x;
		double x_right = grid.nodes[grid.elements[element_number].nodes[1]].x;
		double y_low = grid.nodes[grid.elements[element_number].nodes[0]].y;
		double y_up = grid.nodes[grid.elements[element_number].nodes[3]].y;
		double hx = x_right - x_left, hy = y_up - y_low;
		double ksi = (x - x_left) / hx,
			etta = (y - y_low) / hy;

		return pair<double, double>(ksi, etta);
	}

	pair<double, double> SplineBasis::GetKsiEtta(double x, double y, int element_number, double &hx, double &hy)
	{
		double x_left = grid.nodes[grid.elements[element_number].nodes[0]].x;
		double x_right = grid.nodes[grid.elements[element_number].nodes[1]].x;
		double y_low = grid.nodes[grid.elements[element_number].nodes[0]].y;
		double y_up = grid.nodes[grid.elements[element_number].nodes[3]].y;
		hx = x_right - x_left;
		hy = y_up - y_low;
		double ksi = (x - x_left) / hx,
			etta = (y - y_low) / hy;

		return pair<double, double>(ksi, etta);
	}

	double SplineBasis::phi_i(int i, double x, double y, int element_number)
	{
		double hx, hy;
		pair<double, double> ksi_etta = GetKsiEtta(x, y, element_number, hx, hy);
		return phi[i](ksi_etta.first, ksi_etta.second, hx, hy);
	}

	double SplineBasis::dphi_ix(int i, double x, double y, int element_number)
	{
		double hx, hy;
		pair<double, double> ksi_etta = GetKsiEtta(x, y, element_number, hx, hy);
		return dphiksi[i](ksi_etta.first, ksi_etta.second, hx, hy) / hx;
	}

	double SplineBasis::dphi_iy(int i, double x, double y, int element_number)
	{
		double hx, hy;
		pair<double, double> ksi_etta = GetKsiEtta(x, y, element_number, hx, hy);
		return dphietta[i](ksi_etta.first, ksi_etta.second, hx, hy) / hy;
	}
}