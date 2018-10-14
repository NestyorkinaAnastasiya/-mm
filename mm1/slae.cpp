#include "slae.h"

namespace slae
{
	SLAE::SLAE()
	{
		//Размерность задачи соответствует общему числу базисных функций
		n = grid.nodes.size();

		F.resize(n);
		u.resize(n);
		r.resize(n);
		z.resize(n);
		//Генерация портрета матрицы и её инициализация
		A.CreatePortret();

	}

	double SLAE::dphix(int i, int elementNumber, double ksi, double etta)
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
	}
	
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
				locF[i] += jacobian*tests.Fi(x_, y_) * gaussWeights[k] * basis4.phi[i](ksi, etta);
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
	}

	//Генерация СЛАУ
	void SLAE::GenerateSLAE()
	{
		for (int i = 0; i < n; i++)
			F[i] = 0;
		A.Clear();

		//Высчитывание локальных матриц(векторов) и добавление в глобальные
		for (int i = 0; i < grid.elements.size(); i++)
			CalculateLocals(i);

		//Учёт краевых условий
		for (int i = 0; i < grid.ku[0].size(); i++)
			CalculateBoundaries1(i);
	}

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
	}

	//Вычисление нормы вектора
	double SLAE::Norm(const vector<double> &x)
	{
		double norm = 0;
		int size = x.size();

		for (int i = 0; i < size; i++)
			norm += x[i] * x[i];

		return sqrt(norm);
	}

	//Скалярное произведение векторов
	double SLAE::Scalar(const vector<double> &x, const vector<double> &y)
	{
		double sum = 0;
		int size = x.size();
		for (int i = 0; i < size; i++)
			sum += x[i] * y[i];

		return sum;
	}

	double SLAE::Rel_Discrepancy()
	{
		double dis1, dis2;
		dis1 = Scalar(r, r);
		dis2 = Scalar(F, F);
		double dis = dis1 / dis2;
		return sqrt(dis);
	}

	void SLAE::LULOS()
	{
		double a, b, pp, dis, rr;
		int i, k;
		vector <double> Ax(n), C(n), p(n);

		A.LU();
		//Ax0
		Ax = A * u;
		//f-Ax0
		for (i = 0; i < n; i++)
			r[i] = F[i] - Ax[i];

		//r0=L^(-1)(f-Ax0)
		A.LYF(r, r);

		//z0=U^(-1)r0->r0=Uz0
		A.UXY(r, z);

		//p0=L^(-1)Az0
		Ax = A * z;//Az0
		A.LYF(Ax, p);

		rr = Scalar(r, r);
		dis = Scalar(r, r) / rr;
		dis = sqrt(dis);
		k = 0;

		for (k = 0; dis > eps && k <= maxiter; k++)	{
			//Аk
			pp = Scalar(p, p);
			a = Scalar(p, r) / pp;

			//Xk, Rk
			for (i = 0; i < n; i++)	{
				u[i] = u[i] + a * z[i];
				r[i] = r[i] - a * p[i];
			}

			//UY=rk->Y=U^(-1)rk
			A.UXY(r, C);
			//AU^(-1)rk=Ax
			Ax = A * C;
			//L^(-1)AU^(-1)rk=Y2->L^(-1)B=Y2->LY2=B->Y2=L^(-1)AU^(-1)rk
			A.LYF(Ax, Ax);
			//bk
			b = -Scalar(p, Ax) / pp;

			//zk=U^(-1)rk+bkz[k-1]
			//pk
			for (i = 0; i < n; i++)	{
				z[i] = C[i] + b * z[i];
				p[i] = Ax[i] + b * p[i];
			}
			dis = Scalar(r, r) / rr;
			dis = sqrt(dis);
		}

		dis = Rel_Discrepancy();
		printf("%le\t iter = %d", dis, k);
		getchar();
	}

	void SLAE::Solve()
	{
		FILE *fo;
		GenerateSLAE();
		LULOS();
		fopen_s(&fo, "result.txt", "w");
		fprintf(fo, "x\t\t\ty\t\t\t");
		if (tests.test == 3)
			fprintf(fo, "x*x\t\t\ty*y\t\t\tu\n");
		else fprintf(fo, "u\n");
		for (int i = 0; i < n; i++)
		{
			if (tests.test != 3)
				fprintf(fo, "%lf\t%lf\t%lf\n", grid.nodes[i].x, grid.nodes[i].y, u[i]);
			else
				fprintf(fo, "%lf\t%lf\t%lf\t%lf\t%lf\n", grid.nodes[i].x, grid.nodes[i].y, grid.nodes[i].x * grid.nodes[i].x, grid.nodes[i].y * grid.nodes[i].y, u[i]);

		}
		fprintf(fo, "\n");
		for (int i = 0; i < grid.dofs.size(); i++)
		{		
			fprintf(fo, "%lf\t%lf\n", grid.dofs[i].x, grid.dofs[i].y);
			
		}
		fclose(fo);
	}
}