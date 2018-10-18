#include "solver.h"
namespace solver
{
	Solver::Solver()
	{
		qSpline.resize(slae.n);
		cout << "Building SLAE..." << endl;
		r.resize(slae.n);
		z.resize(slae.n);
	}
	void Solver::Solve()
	{
		slae.A.LU();
		LULOS();
	}

	void Solver::Draw(unsigned int width, unsigned int height, unsigned int countOfIsolines)
	{
		bmp24_file pic(width, height, "plot.bmp");

		// Поиск максимального и минимального значений для отрисовки изолиний
		double maxVal = slae.q[0], minVal = slae.q[0];
		// Поиск максимальной и минимальной координат
		double maxX = grid.nodes[0].x, minX = grid.nodes[0].x, 
			maxY = grid.nodes[0].y, minY = grid.nodes[0].y;
		for (unsigned int i = 1; i < slae.q.size(); i++) {
			if (slae.q[i] > maxVal)	maxVal = slae.q[i];
			if (slae.q[i] < minVal)	minVal = slae.q[i];
			if (grid.nodes[i].x > maxX)	maxX = grid.nodes[i].x;
			if (grid.nodes[i].x < minX)	minX = grid.nodes[i].x;
			if (grid.nodes[i].y > maxY)	maxY = grid.nodes[i].y;
			if (grid.nodes[i].y < minY)	minY = grid.nodes[i].y;
		}

		// Шаги
		double stepX = (maxX - minX) / (double)width, stepY = (maxY - minY) / (double)height;
		cout << "\t> Calculating solution grid..." << endl;
		// Решение по сеткев виде массива
		double ** solutionGrid = new double *[height];

#pragma omp parallel for schedule(dynamic)
		// Цикл по вертикали
		for (int i = 0; i < (int)height; i++) {
			solutionGrid[i] = new double[width];
			// Цикл по горизонтали
			for (unsigned int j = 0; j < width; j++) {
				// Получаем решение в точке
				double x = minX + stepX * (double)j, y = minY + stepY * (double)i;
				int elementNumber = grid.SearchElement(x, y);
				solutionGrid[i][j] = slae.GetSolutionInThePoint(x, y, elementNumber);
			}
		}

		// Белый фон
		for (unsigned int i = 0; i < width; i++) {
			for (unsigned int j = 0; j < height; j++) {
				RGBTRIPLE col = pic.get_pixel(i, j);
				pic.set_pixel(i, j, (BYTE)255, (BYTE)255, (BYTE)255);
			}
		}

		// Сетка
		cout << "\t> Calculating grid lines..." << endl;
		unsigned short colorGrid = 140;

		// Множество координат x
		set<double> gridX;
		for (unsigned int i = 0; i < grid.nodes.size(); i++)
			gridX.insert(grid.nodes[i].x);
		// Множество координат y
		set<double> gridY;
		for (unsigned int i = 0; i < grid.nodes.size(); i++)
			gridY.insert(grid.nodes[i].y);

		// Проходимся по x 
		for (set<double>::iterator i = gridX.begin(); i != gridX.end(); i++) {
			// Проставляем точки по верикали 
			unsigned int coord = (unsigned int)((*i - minX) / stepX);
			for (unsigned int j = 0; j < height; j++)
				pic.set_pixel(coord, j, (BYTE)colorGrid, (BYTE)colorGrid, (BYTE)colorGrid);
		}
		// -//- по y
		for (set<double>::iterator j = gridY.begin(); j != gridY.end(); j++) {
			// Проставляем точки по горизонтали
			unsigned int coord = (unsigned int)((*j - minY) / stepY);
			for (unsigned int i = 0; i < width; i++)
				pic.set_pixel(i, coord, (BYTE)colorGrid, (BYTE)colorGrid, (BYTE)colorGrid);
		}

		// Изолинии
		cout << "\t> Calculating isolines..." << endl;
		unsigned short colorIsolines = 0; //цвет изолиний
		// Шаг, с которым рисовать изолинии
		double isolinesStep = (maxVal - minVal) / (double)(countOfIsolines + 1);

		// Цикл по изолиниям
		for (unsigned int k = 1; k <= countOfIsolines; k++)	{
			// Значение, в котором нужно рисовать изолинию
			double isolineVal = minVal + isolinesStep * (double)k;
			// Цикл по вертикали
			for (unsigned int i = 1; i < height - 1; i++) {
				// Цикл по горизонтали
				for (unsigned int j = 1; j < width - 1; j++) {
					// Если в точке i, j нужно отметить пиксель изолинии, то отмечаем
					if ((solutionGrid[i][j] >= isolineVal && solutionGrid[i + 1][j] < isolineVal) ||
						(solutionGrid[i][j] <= isolineVal && solutionGrid[i + 1][j] > isolineVal) ||
						(solutionGrid[i][j] >= isolineVal && solutionGrid[i][j + 1] < isolineVal) ||
						(solutionGrid[i][j] <= isolineVal && solutionGrid[i][j + 1] > isolineVal)) {
						pic.set_pixel(j, i,	(BYTE)colorIsolines, (BYTE)colorIsolines,	(BYTE)colorIsolines);
					}
				}
			}
		}

		// Вывод изображения в файл
		pic.write();

		for (unsigned int i = 0; i < height; i++)
			delete[] solutionGrid[i];
		delete[] solutionGrid;
	}
	void Solver::Draw(unsigned int width, unsigned int height, unsigned int countOfIsolines, set<int> &nd)
	{
		
		class bmp24_file pic(width, height, "plot_spline.bmp");
		cout << endl << "Drawing picture..." << endl;
		auto iter = nd.begin();
		// Поиск максимального и минимального значений для отрисовки изолиний
		double maxVal = slae.q[0], minVal = slae.q[0];
		// Поиск максимальной и минимальной координат
		double maxX = grid.nodes[0].x, minX = grid.nodes[0].x,
			maxY = grid.nodes[0].y, minY = grid.nodes[0].y;
		unsigned int size = nd.size();
		iter++;
		for (unsigned int i = 1; i < size; i++) {
			if (slae.q[*iter] > maxVal)	maxVal = slae.q[*iter];
			if (slae.q[*iter] < minVal)	minVal = slae.q[*iter];
			if (grid.nodes[i].x > maxX)	maxX = grid.nodes[i].x;
			if (grid.nodes[i].x < minX)	minX = grid.nodes[i].x;
			if (grid.nodes[i].y > maxY)	maxY = grid.nodes[i].y;
			if (grid.nodes[i].y < minY)	minY = grid.nodes[i].y;
			iter++;
		}
		// Шаги
		double stepX = (maxX - minX) / (double)width, stepY = (maxY - minY) / (double)height;
		cout << "\t> Calculating solution grid..." << endl;
		// Решение по сеткев виде массива
		double ** solutionGrid = new double *[height];

#pragma omp parallel for schedule(dynamic)
		// Цикл по вертикали
		for (int i = 0; i < (int)height; i++) {
			solutionGrid[i] = new double[width];
			// Цикл по горизонтали
			for (unsigned int j = 0; j < width; j++) {
				// Получаем решение в точке
				double x = minX + stepX * (double)j, y = minY + stepY * (double)i;
				int elementNumber = grid.SearchElement(x, y);
				solutionGrid[i][j] = slae.GetSolutionInThePointWithSpline(x, y, elementNumber);
			}
		}

		// Белый фон
		for (unsigned int i = 0; i < width; i++) {
			for (unsigned int j = 0; j < height; j++) {
				RGBTRIPLE col = pic.get_pixel(i, j);
				pic.set_pixel(i, j, (BYTE)255, (BYTE)255, (BYTE)255);
			}
		}

		// Сетка
		cout << "\t> Calculating grid lines..." << endl;
		unsigned short colorGrid = 140;

		// Множество координат x
		set<double> gridX;
		for (unsigned int i = 0; i < grid.nodes.size(); i++)
			gridX.insert(grid.nodes[i].x);
		// Множество координат y
		set<double> gridY;
		for (unsigned int i = 0; i < grid.nodes.size(); i++)
			gridY.insert(grid.nodes[i].y);

		// Проходимся по x 
		for (set<double>::iterator i = gridX.begin(); i != gridX.end(); i++) {
			// Проставляем точки по верикали 
			unsigned int coord = (unsigned int)((*i - minX) / stepX);
			for (unsigned int j = 0; j < height; j++)
				pic.set_pixel(coord, j, (BYTE)colorGrid, (BYTE)colorGrid, (BYTE)colorGrid);
		}
		// -//- по y
		for (set<double>::iterator j = gridY.begin(); j != gridY.end(); j++) {
			// Проставляем точки по горизонтали
			unsigned int coord = (unsigned int)((*j - minY) / stepY);
			for (unsigned int i = 0; i < width; i++)
				pic.set_pixel(i, coord, (BYTE)colorGrid, (BYTE)colorGrid, (BYTE)colorGrid);
		}
		// Изолинии
		cout << "\t> Calculating isolines..." << endl;
		unsigned short colorIsolines = 0; //цвет изолиний
		// Шаг, с которым рисовать изолинии
		double isolinesStep = (maxVal - minVal) / (double)(countOfIsolines + 1);

		// Цикл по изолиниям
		for (unsigned int k = 1; k <= countOfIsolines; k++) {
			// Значение, в котором нужно рисовать изолинию
			double isolineVal = minVal + isolinesStep * (double)k;
			// Цикл по вертикали
			for (unsigned int i = 1; i < height - 1; i++) {
				// Цикл по горизонтали
				for (unsigned int j = 1; j < width - 1; j++) {
					// Если в точке i, j нужно отметить пиксель изолинии, то отмечаем
					if ((solutionGrid[i][j] >= isolineVal && solutionGrid[i + 1][j] < isolineVal) ||
						(solutionGrid[i][j] <= isolineVal && solutionGrid[i + 1][j] > isolineVal) ||
						(solutionGrid[i][j] >= isolineVal && solutionGrid[i][j + 1] < isolineVal) ||
						(solutionGrid[i][j] <= isolineVal && solutionGrid[i][j + 1] > isolineVal)) {
						pic.set_pixel(j, i, (BYTE)colorIsolines, (BYTE)colorIsolines, (BYTE)colorIsolines);
					}
				}
			}
		}

		// Вывод изображения в файл
		pic.write();

		for (unsigned int i = 0; i < height; i++)
			delete[] solutionGrid[i];
		delete[] solutionGrid;
	}
	//Вычисление нормы вектора
	double Solver::Norm(const vector<double> &x)
	{
		double norm = 0;
		int size = x.size();

		for (int i = 0; i < size; i++)
			norm += x[i] * x[i];

		return sqrt(norm);
	}
	//Скалярное произведение векторов
	double Solver::Scalar(const vector<double> &x, const vector<double> &y)
	{
		double sum = 0;
		int size = x.size();
		for (int i = 0; i < size; i++)
			sum += x[i] * y[i];

		return sum;
	}
	double Solver::RelDiscrepancy()
	{
		double dis1, dis2;
		dis1 = Scalar(r, r);
		dis2 = Scalar(slae.F, slae.F);
		double dis = dis1 / dis2;
		return sqrt(dis);
	}
	void Solver::LULOS()
	{
		double a, b, pp, dis, rr;
		int i, k;
		vector <double> Ax(slae.n), C(slae.n), p(slae.n);
		//Ax0
		Ax = slae.A * slae.q;
		//f-Ax0
		for (i = 0; i < slae.n; i++)
			r[i] = slae.F[i] - Ax[i];

		//r0=L^(-1)(f-Ax0)
		slae.A.LYF(r, r);

		//z0=U^(-1)r0->r0=Uz0
		slae.A.UXY(r, z);

		//p0=L^(-1)Az0
		Ax = slae.A * z;//Az0
		slae.A.LYF(Ax, p);

		rr = Scalar(r, r);
		dis = Scalar(r, r) / rr;
		dis = sqrt(dis);
		k = 0;

		for (k = 0; dis > epsilon && k <= maxCountOfIterations; k++) {
			//Аk
			pp = Scalar(p, p);
			a = Scalar(p, r) / pp;

			//Xk, Rk
			for (i = 0; i < slae.n; i++) {
				slae.q[i] = slae.q[i] + a * z[i];
				r[i] = r[i] - a * p[i];
			}

			//UY=rk->Y=U^(-1)rk
			slae.A.UXY(r, C);
			//AU^(-1)rk=Ax
			Ax = slae.A * C;
			//L^(-1)AU^(-1)rk=Y2->L^(-1)B=Y2->LY2=B->Y2=L^(-1)AU^(-1)rk
			slae.A.LYF(Ax, Ax);
			//bk
			b = -Scalar(p, Ax) / pp;

			//zk=U^(-1)rk+bkz[k-1]
			//pk
			for (i = 0; i < slae.n; i++) {
				z[i] = C[i] + b * z[i];
				p[i] = Ax[i] + b * p[i];
			}
			dis = Scalar(r, r) / rr;
			dis = sqrt(dis);
		}

		dis = RelDiscrepancy();
		printf("%le\t iter = %d", dis, k);
		getchar();
		
	}
}
