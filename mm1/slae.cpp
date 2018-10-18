#include "slae.h"
namespace slae
{

	SLAE::SLAE()
	{
		//Размерность задачи соответствует общему числу базисных функций
		n = grid.nodes.size();
		F.resize(n);
		q.resize(n);
		//Генерация портрета матрицы и её инициализация
		A.CreatePortret();

	}
	void SLAE::GenerateSLAE()
	{
		int size = grid.elements.size();
		for (int i = 0; i < size; i++)
			CalculateLocals(i);
	}
	void SLAE::CalculateLocals(int elementNumber)
	{
		Element element = grid.elements[elementNumber];

		// Локальные точки по которым считается сплайн		
		//const int local_nodes_num = 4;
		//const Point local_nodes[local_nodes_num] =
		//{
		//	Point(0.0, 0.0),
		//	Point(1.0, 0.0),
		//	Point(0.0, 1.0),
		//	Point(1.0, 1.0)
		//};

		const size_t localNodesNum = 16;
		const Point localNodes[localNodesNum] =
		{
			Point(0.0, 0.0),
			Point(1.0 / 3.0, 0.0),
			Point(2.0 / 3.0, 0.0),
			Point(1.0, 0.0),
			Point(0.0, 1.0 / 3.0),
			Point(1.0 / 3.0, 1.0 / 3.0),
			Point(2.0 / 3.0, 1.0 / 3.0),
			Point(1.0, 1.0 / 3.0),
			Point(0.0, 2.0 / 3.0),
			Point(1.0 / 3.0, 2.0 / 3.0),
			Point(2.0 / 3.0, 2.0 / 3.0),
			Point(1.0, 2.0 / 3.0),
			Point(0.0, 1.0),
			Point(1.0 / 3.0, 1.0),
			Point(2.0 / 3.0, 1.0),
			Point(1.0, 1.0)
		};

		CalculateA(elementNumber, localNodes, localNodesNum);
		CalculateAlphaMatrix(elementNumber);
		CalculateBettaMatrix(elementNumber);
		CalculateF(elementNumber, localNodes, localNodesNum);
	}
	void SLAE::CalculateA(int elementNumber, const Point * points, int kPoints)
	{
		Element element = grid.elements[elementNumber];

		vector <vector<double>> AA;
		AA.resize(nSplineFunc);
		for (int i = 0; i < nSplineFunc; i++)
			AA[i].resize(nSplineFunc);

		double hx = grid.hx(elementNumber);
		double hy = grid.hy(elementNumber);

		for (int i = 0; i < nSplineFunc; i++) {
			for (int j = i; j < nSplineFunc; j++) {
				// Цикл по точкам
				for (int l = 0; l < kPoints; l++) {
					double p_ksi = points[l].x,
						p_etta = points[l].y;
					AA[i][j] += splBasis.phi[i](p_ksi, p_etta, hx, hy) * splBasis.phi[j](p_ksi, p_etta, hx, hy);
				}
			}
		}

		for (int i = 1; i < nSplineFunc; i++)
			for (int j = 0; j < i; j++)
				AA[i][j] = AA[j][i];

		for (int i = 0; i < nSplineFunc; i++) {
			int id_i = element.dof[i];
			for (int j = 0; j < nSplineFunc; j++) {
				int id_j = element.dof[j];
				A.AddElement(id_i, id_j, AA[i][j]);
			}
		}
	}
	void SLAE::CalculateAlphaMatrix(int elementNumber)
	{
		vector <vector<double>> AM;
		AM.resize(nSplineFunc);
		for (int i = 0; i < nSplineFunc; i++)
			AM[i].resize(nSplineFunc);
		Element element = grid.elements[elementNumber];

		double hx = grid.hx(elementNumber);
		double hy = grid.hy(elementNumber);
		double hx2 = hx * hx;
		double hy2 = hy * hy;

		double jacobian = hx * hy / 4.0;

		for (int i = 0; i < nSplineFunc; i++) {
			for (int j = i; j < nSplineFunc; j++) {
				double am1 = 0, am2 = 0;
				for (int k = 0; k < 25; k++) {
					double p_ksi = 0.5 + 0.5 * gaussPoints[0][k],
						p_etta = 0.5 + 0.5 * gaussPoints[0][k];
					am1 += gaussWeights[k] *
						splBasis.dphiksi[i](p_ksi, p_etta, hx, hy) * splBasis.dphiksi[j](p_ksi, p_etta, hx, hy);
					am2 += gaussWeights[k] *
						splBasis.dphietta[i](p_ksi, p_etta, hx, hy) * splBasis.dphietta[j](p_ksi, p_etta, hx, hy);
				}
				AM[i][j] = (am1 / hx2 + am2 / hy2) * jacobian * alpha;
			}
		}

		for (int i = 1; i < 4; i++)
			for (int j = 0; j < i; j++)
				AM[i][j] = AM[j][i];

		for (int i = 0; i < nSplineFunc; i++) {
			int id_i = element.dof[i];
			for (int j = 0; j < nSplineFunc; j++) {
				int id_j = element.dof[j];
				A.AddElement(id_i, id_j, AM[i][j]);
			}
		}
	}
	void SLAE::CalculateBettaMatrix(int elementNumber)
	{
		vector <vector<double>> BM;
		BM.resize(nSplineFunc);
		for (int i = 0; i < nSplineFunc; i++)
			BM[i].resize(nSplineFunc);
		Element element = grid.elements[elementNumber];

		double hx = grid.hx(elementNumber);
		double hy = grid.hy(elementNumber);
		double hx2 = hx * hx;
		double hy2 = hy * hy;

		double jacobian = hx * hy / 4.0;

		for (int i = 0; i < nSplineFunc; i++) {
			for (int j = i; j < nSplineFunc; j++) {
				double bm = 0;
				for (int k = 0; k < 25; k++) {
					double p_ksi = 0.5 + 0.5 * gaussPoints[0][k],
						p_etta = 0.5 + 0.5 * gaussPoints[0][k];
					bm += gaussWeights[k] *
						(splBasis.d2phiksi[i](p_ksi, p_etta, hx, hy) / hx + splBasis.d2phietta[i](p_ksi, p_etta, hx, hy) / hy) *
						(splBasis.d2phiksi[j](p_ksi, p_etta, hx, hy) / hx + splBasis.d2phietta[j](p_ksi, p_etta, hx, hy) / hy);
				}
				BM[i][j] = bm * jacobian * betta;
			}
		}

		for (int i = 1; i < 4; i++)
			for (int j = 0; j < i; j++)
				BM[i][j] = BM[j][i];

		for (int i = 0; i < nSplineFunc; i++) {
			int id_i = element.dof[i];
			for (int j = 0; j < nSplineFunc; j++) {
				int id_j = element.dof[j];
				A.AddElement(id_i, id_j, BM[i][j]);
			}
		}
	}
	void SLAE::CalculateF(int elementNumber, const Point * points, int kPoints)
	{
		Element element = grid.elements[elementNumber];

		vector<double> F;
		F.resize(nSplineFunc);

		double x0 = grid.nodes[element.nodes[0]].x;
		double y0 = grid.nodes[element.nodes[0]].y;

		double hx = grid.hx(elementNumber);
		double hy = grid.hy(elementNumber);

		for (int i = 0; i < nSplineFunc; i++) {
			// Цикл по точкам
			for (int l = 0; l < kPoints; l++) {
				double p_ksi = points[l].x,
					p_etta = points[l].y;
				double f = GetSolutionInThePoint(hx * p_ksi + x0, hy * p_etta + y0, elementNumber);
				F[i] += splBasis.phi[i](p_ksi, p_etta, hx, hy) * f;
			}
		}

		for (int i = 0; i < nSplineFunc; i++) {
			int id_i = element.dof[i];
			F[id_i] += F[i];
		}
	}
	double SLAE::GetSolutionInThePoint(double x, double y, int elementNumber)
	{
		vector <int> indexes;
		indexes.resize(4);

		vector <double> qLocal;
		qLocal.resize(4);

		//собираем глобальные номера с элемента
		for (int j = 0; j < 4; j++)
			indexes[j] = grid.elements[elementNumber].nodes[j];

		//собираем локальный набор весов
		for (int j = 0; j < 4; j++)
			qLocal[j] = q[indexes[j]];

		double solutionInThePoint = 0;
		for (int j = 0; j < 4; j++)
			solutionInThePoint += qLocal[j] * splBasis.phi_i(j, x, y, elementNumber);

		//return x*x*x + y*y*y;
		return solutionInThePoint;//??
	}
	double SLAE::GetSolutionInThePointWithSpline(double x, double y, int elementNumber)
	{
		vector <int> indexes;
		indexes.resize(nSplineFunc);

		vector <double> qLocal;
		qLocal.resize(nSplineFunc);

		//собираем глобальные номера с элемента
		for (int j = 0; j < nSplineFunc; j++)
			indexes[j] = grid.elements[elementNumber].dof[j];

		//собираем локальный набор весов
		for (int j = 0; j < nSplineFunc; j++)
			qLocal[j] = q[indexes[j]];

		double solutionInThePoint = 0;
		for (int j = 0; j < nSplineFunc; j++)
			solutionInThePoint += qLocal[j] * splBasis.phi_i(j, x, y, elementNumber);

		return solutionInThePoint; // кажется здесь я заебалась
	}
	/* какая-то хрень с производными*/
	/*double SLAE::dphix(int i, int elementNumber, double ksi, double etta)
	{
		double yy[4];
		for (int j = 0; j < 4; j++)
			yy[j] = grid.nodes[grid.elements[elementNumber].nodes[j]].y;

		double b3 = yy[2] - yy[0];
		double b4 = yy[1] - yy[0];
		double b6 = yy[0] - yy[1] - yy[2] + yy[3];

		double J_inv = 1 / GetAbsJ(elementNumber, ksi, etta);
		double dphiix = J_inv * (basis4.dphiksi[i](ksi, etta) * (b6 * ksi + b3)
			- basis4.dphietta[i](ksi, etta) * (b6 * etta + b4));
		return dphiix;
	}

	double SLAE::dphiy(int i, int elementNumber, double ksi, double etta)
	{
		double xx[4];
		for (int i = 0; i < 4; i++)
			xx[i] = grid.nodes[grid.elements[elementNumber].nodes[i]].x;

		double b1 = xx[2] - xx[0];
		double b2 = xx[1] - xx[0];
		double b5 = xx[0] - xx[1] - xx[2] + xx[3];

		double J_inv = 1 / GetAbsJ(elementNumber, ksi, etta);
		double dphiiy = J_inv * (basis4.dphietta[i](ksi, etta) * (b5 * etta + b2)
			- basis4.dphiksi[i](ksi, etta) * (b5 * ksi + b1));
		return dphiiy;
	}

	double SLAE::GetAbsJ(int elementNumber, double ksi, double etta)
	{

		Element element = grid.elements[elementNumber];

		double xx[4], yy[4], J[4] = { 0,0,0,0 };
		for (int i = 0; i < 4; i++)
		{
			xx[i] = grid.nodes[element.nodes[i]].x;
			yy[i] = grid.nodes[element.nodes[i]].y;
		}

		for (int i = 0; i < 4; i++)
			J[0] += xx[i] * basis4.dphiksi[i](ksi, etta);
		for (int i = 0; i < 4; i++)
			J[1] += yy[i] * basis4.dphiksi[i](ksi, etta);
		for (int i = 0; i < 4; i++)
			J[2] += xx[i] * basis4.dphietta[i](ksi, etta);
		for (int i = 0; i < 4; i++)
			J[3] += yy[i] * basis4.dphietta[i](ksi, etta);

		double jacobian = J[0] * J[3] - J[1] * J[2];
		return fabs(jacobian);
	}

	double SLAE::GetAbsBQJ(int elementNumber, double ksi, double etta)
	{
		Element element = grid.elements[elementNumber];

		double xx[9], yy[9], J[4] = { 0,0,0,0 };
		for (int i = 0; i < 9; i++)
		{
			xx[i] = grid.dofs[element.dof[i]].x;
			yy[i] = grid.dofs[element.dof[i]].y;
		}

		for (int i = 0; i < 9; i++)
			J[0] += xx[i] * basis9.dphiksi[i](ksi, etta);
		for (int i = 0; i < 9; i++)
			J[1] += yy[i] * basis9.dphiksi[i](ksi, etta);
		for (int i = 0; i < 9; i++)
			J[2] += xx[i] * basis9.dphietta[i](ksi, etta);
		for (int i = 0; i < 9; i++)
			J[3] += yy[i] * basis9.dphietta[i](ksi, etta);

		double jacobian = J[0] * J[3] - J[1] * J[2];
		return fabs(jacobian);
	}*/
	/* сборка старых лок. матриц*/
	/*
	//Сборка локальных матриц жёсткости
	void SLAE::CalculateG(int elementNumber)
	{
		Element element = grid.elements[elementNumber];
		double ksi, etta, lambda, jacobian, gij, g1, g2;

		for (int i = 0; i < 4; i++)
		{
			for (int j = i; j < 4; j++)
			{
				gij = 0;
				for (int k = 0; k < 25; k++)
				{
					ksi = 0.5 + 0.5 * gaussPoints[0][k]; etta = 0.5 + 0.5 * gaussPoints[1][k];
					lambda = tests.Lambda();
					jacobian = GetAbsBQJ(elementNumber, ksi, etta);
					g1 = dphix(i, elementNumber, ksi, etta) * dphix(j, elementNumber, ksi, etta);
					g2 = dphiy(i, elementNumber, ksi, etta) * dphiy(j, elementNumber, ksi, etta);
					gij += gaussWeights[k] * (g1 + g2) * jacobian * lambda;
				}
				G[i][j] = gij * 0.25;
			}
		}
		//матрица симметричная, заполняем нижний треугольник
		for (int i = 1; i < 4; i++)
			for (int j = 0; j < i; j++)
				G[i][j] = G[j][i];
	}

	//Сборка локальных матриц масс
	void SLAE::CalculateM(int elementNumber)
	{
		Element element = grid.elements[elementNumber];
		double  m, ksi, etta, sigma, jacobian;
		for (int i = 0; i < 4; i++)
		{
			for (int j = i; j < 4; j++)
			{
				m = 0;
				for (int k = 0; k < 25; k++)
				{
					ksi = 0.5 + 0.5 * gaussPoints[0][k]; etta = 0.5 + 0.5 * gaussPoints[1][k];
					sigma = tests.Sigma();
					jacobian = GetAbsBQJ(elementNumber, ksi, etta);
					m += jacobian * sigma * gaussWeights[k]
						* basis4.phi[i](ksi, etta) * basis4.phi[j](ksi, etta);
				}
				M[i][j] = m * 0.25;
			}
		}
		//матрица симметричная, заполняем нижний треугольник
		for (int i = 1; i < 4; i++)
			for (int j = 0; j < i; j++)
				M[i][j] = M[j][i];
	}

	//Сборка локальных правых частей
	void SLAE::CalculateLocalF(int elementNumber)
	{
		Element element = grid.elements[elementNumber];
		double ksi, etta, xx[9], yy[9], x_, y_, jacobian;
		for (int i = 0; i < 9; i++)
		{
			xx[i] = grid.dofs[element.dof[i]].x;
			yy[i] = grid.dofs[element.dof[i]].y;
		}
		//интегрирование(Гаусс 5)
		for (int i = 0; i < 4; i++)
		{
			locF[i] = 0;
			for (int k = 0; k < 25; k++)
			{
				ksi = 0.5 + 0.5 * gaussPoints[0][k]; etta = 0.5 + 0.5 * gaussPoints[1][k];
				x_ = 0, y_ = 0;
				for (int j = 0; j < 9; j++)
				{
					x_ += xx[j] * basis9.phi[j](ksi, etta);
					y_ += yy[j] * basis9.phi[j](ksi, etta);
				}
				jacobian = GetAbsBQJ(elementNumber, ksi, etta);
				locF[i] += jacobian * tests.Fi(x_, y_) * gaussWeights[k] * basis4.phi[i](ksi, etta);
			}
			locF[i] *= 0.25;
		}
	}

	//Сборка локальных матриц(векторов) и добавление в глобальные
	void SLAE::CalculateLocals(int elementNumber)
	{
		Element element = grid.elements[elementNumber];
		int ki, kj;

		//вычисление локальных матриц
		CalculateG(elementNumber);
		CalculateM(elementNumber);
		CalculateLocalF(elementNumber);

		for (int i = 0; i < 4; i++)
		{
			ki = element.nodes[i];
			for (int j = 0; j < 4; j++)
			{
				kj = element.nodes[j];
				//добавка в глобальную матрицу А
				AddElementToGlobalMatrix(A, ki, kj, G[i][j] + M[i][j]);
			}
			//добавка в глобальную правую часть
			F[ki] += locF[i];
		}
	}*/

	/*
	//Нахождение правой части для 1ого краевого условия
	void SLAE::Calculate_g(int formNumber, int orientation, int elNumber)
	{
		Element element = grid.elements[elNumber];
		double x1 = 0, x2 = 0, y1 = 0, y2 = 0;

		switch (orientation)
		{
			//левое ребро
		case 0:
		{
			x1 = grid.nodes[element.nodes[0]].x;
			x2 = grid.nodes[element.nodes[2]].x;
			y1 = grid.nodes[element.nodes[0]].y;
			y2 = grid.nodes[element.nodes[2]].y;

		}
		break;
		//правое ребро
		case 1:
		{
			x1 = grid.nodes[element.nodes[1]].x;
			x2 = grid.nodes[element.nodes[3]].x;
			y1 = grid.nodes[element.nodes[1]].y;
			y2 = grid.nodes[element.nodes[3]].y;
		}
		break;
		//нижнее ребро
		case 2:
		{
			x1 = grid.nodes[element.nodes[0]].x;
			x2 = grid.nodes[element.nodes[1]].x;
			y1 = grid.nodes[element.nodes[0]].y;
			y2 = grid.nodes[element.nodes[1]].y;
		}
		break;
		//верхнее ребро
		case 3:
		{
			x1 = grid.nodes[element.nodes[2]].x;
			x2 = grid.nodes[element.nodes[3]].x;
			y1 = grid.nodes[element.nodes[2]].y;
			y2 = grid.nodes[element.nodes[3]].y;
		}
		break;
		default:; break;
		}
		g[0] = tests.Ug(x1, y1);
		g[1] = tests.Ug(x2, y2);
	}

	//Вычисление 1ого краевого условия для одного узла
	void SLAE::CalculateBoundaries1ForNode(int node, double gi, double weight)
	{
		int id;
		F[node] = gi;
		A.di[node] = weight;

		for (int j = 0; j < n; j++)
			if (node < j)
			{
				bool flag = false;
				for (id = A.ig[j]; !flag && id <= A.ig[j + 1] - 1; id++)
					if (A.jg[id] == node) flag = true;
				if (flag) A.ggu[id - 1] = 0.0;
			}
			else
			{
				bool flag = false;
				for (id = A.ig[node]; !flag && id <= A.ig[node + 1] - 1; id++)
					if (A.jg[id] == j) flag = true;
				if (flag) A.ggl[id - 1] = 0.0;
			}
	}

	//Учёт первого краевого условия
	void SLAE::CalculateBoundaries1(int number)
	{
		Element element = grid.elements[grid.ku[0][number].elem];

		if (grid.ku[0][number].edges[0] == 1)
		{
			int indexes[2] = { element.nodes[0], element.nodes[2] };
			Calculate_g(grid.ku[0][number].formNumber[0], 0, grid.ku[0][number].elem);
			for (int i = 0; i < 2; i++)
				CalculateBoundaries1ForNode(indexes[i], g[i], 1);
		}
		if (grid.ku[0][number].edges[1] == 1)
		{
			int indexes[2] = { element.nodes[1], element.nodes[3] };
			Calculate_g(grid.ku[0][number].formNumber[1], 1, grid.ku[0][number].elem);
			for (int i = 0; i < 2; i++)
				CalculateBoundaries1ForNode(indexes[i], g[i], 1);
		}
		if (grid.ku[0][number].edges[2] == 1)
		{
			int indexes[2] = { element.nodes[0], element.nodes[1] };
			Calculate_g(grid.ku[0][number].formNumber[2], 2, grid.ku[0][number].elem);
			for (int i = 0; i < 2; i++)
				CalculateBoundaries1ForNode(indexes[i], g[i], 1);
		}
		if (grid.ku[0][number].edges[3] == 1)
		{
			int indexes[2] = { element.nodes[2], element.nodes[3] };
			Calculate_g(grid.ku[0][number].formNumber[3], 3, grid.ku[0][number].elem);
			for (int i = 0; i < 2; i++)
				CalculateBoundaries1ForNode(indexes[i], g[i], 1);
		}
	}*/

}