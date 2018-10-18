#include "matrix.h"
namespace matrix
{
	//Генерация портрета матрицы
	void Matrix::CreatePortret()
	{
		std::cout << "Creating matrix-portret..." << endl;
		// Общее число степеней свободы
		int slaeSize = grid.elements.size() * grid.elements[0].nDof;
		int ndof = grid.elements[0].nDof;
		vector <int> elList, unzeroNumbersList;

		vector<vector <int>> list;
		list.reserve(slaeSize);

		for (int i = 0; i < slaeSize; i++) {
			//нашли,в каких элементах присутствует данная i-я степень свободы
			grid.SearchElements(i, elList);
			//собрать список узлов этих элементов, меньших по номеру
			for (unsigned int j = 0; j < elList.size(); j++) {
				for (int k = 0; k < ndof; k++) {
					//если такого номера ещё нет и он меньше по номеру,то добавляем
					if (find(unzeroNumbersList.begin(), unzeroNumbersList.end(), grid.elements[elList[j]].nodes[k])
						== unzeroNumbersList.end() && grid.elements[elList[j]].nodes[k] < i)
						unzeroNumbersList.push_back(grid.elements[elList[j]].nodes[k]);
				}
			}
			sort(unzeroNumbersList.begin(), unzeroNumbersList.end());
			list.push_back(unzeroNumbersList);
			unzeroNumbersList.clear();
			elList.clear();
		}

		//вычисляем размерность ggl,ggu
		int gg_size = 0;
		for (int i = 0; i < slaeSize; i++) {
			if (!list[i].empty())
				gg_size += list[i].size();
		}

		//инициализируем матрицу и собираем портрет
		Initialize(slaeSize, gg_size);

		ig[0] = 0;

		for (int i = 0; i < n; i++) {
			if (!list[i].empty())
				ig[i + 1] = ig[i] + list[i].size();
			else ig[i + 1] = ig[i];
		}

		int k = 0;
		for (int i = 0; i < n; i++) {
			if (!list[i].empty()) {
				for (unsigned int j = 0; j < list[i].size(); j++) {
					jg[k] = list[i][j];
					k++;
				}
				list[i].clear();
			}
		}

		list.clear();
	}

	//Инициализация матрицы после генерации портрета
	void Matrix::Initialize(int size1, int size2)
	{
		n = size1; size = size2;

		ggl.resize(size);
		ggu.resize(size);
		di.resize(n);
		ig.resize(n + 1);
		jg.resize(size);
		L.resize(size);
		D.resize(n);
		U.resize(size);
	}

	vector<double> Matrix::operator*(vector<double> x)
	{
		int i, j, l, ik, iend, k;
		vector <double> result;
		result.resize(x.size());

		for (i = 0; i < n; i++) {
			//начало i-ой строки(столбца)
			l = ig[i];
			//начало (i+1)-ой строки(столбца)
			iend = ig[i + 1];
			//количество элементов в i строке(столбце)
			ik = iend - l;

			result[i] = di[i] * x[i];

			//проходим по всем элементам i строки (столбца)
			for (k = 0; k < ik; k++, l++) {
				j = jg[l];
				result[i] += ggl[l] * x[j];
				result[j] += ggu[l] * x[i];
			}
		}
		return result;
	}

	void Matrix::Clear()
	{
		for (int i = 0; i < n; i++)
			di[i] = 0;
		for (int i = 0; i < ggl.size(); i++)
			ggl[i] = ggu[i] = 0;
	}

	//Умножение транспонированной матрицы на вектор
	void Matrix::MultiplyATx(vector <double> a, vector <double> &result)
	{
		int i, j, l, ik, iend, k;
		for (i = 0; i < n; i++)
		{
			//начало i-ой строки(столбца)
			l = ig[i];
			//начало (i+1)-ой строки(столбца)
			iend = ig[i + 1];
			//количество элементов в i строке(столбце)
			ik = iend - l;

			result[i] = di[i] * a[i];

			//проходим по всем элементам i строки (столбца)
			for (k = 0; k < ik; k++, l++)
			{
				j = jg[l];
				result[i] += ggu[l] * a[j];
				result[j] += ggl[l] * a[i];
			}
		}
	}

	// Добавить элемент в матрицу
	void Matrix::AddElement(int i, int j, double element)
	{
		int id;
		bool flag;

		if (i == j)	di[i] += element;
		else {
			if (i < j) {
				flag = false;
				for (id = ig[j]; !flag && id < ig[j + 1]; id++)
					if (jg[id] == i) flag = true;

				if (flag) ggu[id - 1] += element;
			}
			else {
				flag = false;
				for (id = ig[i]; !flag && id < ig[i + 1]; id++)
					if (jg[id] == j) flag = true;

				if (flag) ggl[id - 1] += element;
			}
		}
	}

	// Заменить элемент в матрице
	void Matrix::ChangeElement(int i, int j, double element)
	{
		int id;
		bool flag;

		if (i == j)	di[i] = element;
		else {
			if (i < j) {
				flag = false;
				for (id = ig[j]; !flag && id < ig[j + 1]; id++)
					if (jg[id] == i) flag = true;

				if (flag) ggu[id - 1] = element;
			}
			else {
				flag = false;
				for (id = ig[i]; !flag && id < ig[i + 1]; id++)
					if (jg[id] == j) flag = true;

				if (flag) ggl[id - 1] = element;
			}
		}
	}

	// Сложение элементов в строке (для проверки)
	void Matrix::Sum(vector <double> &result)
	{
		int i, j, l, ik, iend, k;

		for (i = 0; i < n; i++) {
			// начало i-ой строки(столбца)
			l = ig[i];
			// начало (i+1)-ой строки(столбца)
			iend = ig[i + 1];
			// количество элементов в i строке(столбце)
			ik = iend - l;

			result[i] = di[i];

			// проходим по всем элементам i строки (столбца)
			for (k = 0; k < ik; k++, l++) {
				j = jg[l];
				result[i] += ggl[l];
				result[j] += ggu[l];
			}
		}
	}
	void Matrix::LU()
	{
		int i, i0, j0, iend, num, ki, kj, jend;
		double suml, sumu, sumdg;

		L = ggl;
		U = ggu;
		D = di;

		for (i = 0; i < n; i++) {
			i0 = ig[i];
			iend = ig[i + 1];

			for (num = i0, sumdg = 0; num < iend; num++) {
				j0 = ig[jg[num]];
				jend = ig[jg[num] + 1];
				ki = i0;
				kj = j0;
				//для num учитываются все предыдущие элементы
				for (suml = 0, sumu = 0, ki = i0; ki < num; ki++) {
					for (int m = kj; m < jend; m++)
						//ищем соответствующие ненулевые элементы для умножения
						if (jg[ki] == jg[m]) {
							suml += L[ki] * U[m];
							sumu += L[m] * U[ki];
						}
				}
				L[num] -= suml;
				U[num] = (U[num] - sumu) / D[jg[num]];
				//умножаются симметричные элементы
				sumdg += L[num] * U[num];
			}
			D[i] -= sumdg;
		}
	}
	void Matrix::LYF(const vector<double>& C, vector<double>& yl)
	{
		int i, i0, iend; //i0-адрес начала строки, iend-адрес конца строки
		double sum;
		for (i = 0; i < n; i++)
		{
			i0 = ig[i]; iend = ig[i + 1];

			yl[i] = C[i];

			for (i0, sum = 0; i0 < iend; i0++)
				yl[i] -= yl[jg[i0]] * L[i0];
			yl[i] /= D[i];
		}
	}
	void Matrix::UXY(const vector<double>& C, vector<double>& yu)
	{
		int i, i0, iend;

		for (i = 0; i < n; i++)
			yu[i] = 0.0;

		for (i = n - 1; i >= 0; i--)//проход по столбцам с конца
		{
			yu[i] += C[i];

			i0 = ig[i]; iend = ig[i + 1]; iend--;
			for (; iend >= i0; iend--)//идём по столбцу с конца
				yu[jg[iend]] -= yu[i] * U[iend];
		}
	}
}