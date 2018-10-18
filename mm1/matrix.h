#pragma once
#pragma once
#include <algorithm>
#include "basis.h"

namespace matrix
{
	class Matrix
	{
		// Матрица 
		// Размерность матрицы
		int n;
		// Размер векторов ggl, ggu, jg
		int size;

		// Указатели начала строк(столбцов)
		vector <int> ig;
		// Номера столбцов внедиагональных элементов 
		vector <int> jg;
		// Диагональные элементы матрицы
		vector <double> di;
		// Внедиагональные элементы нижнего треугольника матрицы
		vector <double> ggl;
		// Внедиагональные элементы верхнего треугольника матрицы
		vector <double> ggu;

		//Компоненты матрицы с факторизацией
		vector <double> L;
		vector <double> D;
		vector <double> U;

	public:
		// Генерация портрета матрицы
		void CreatePortret();
		// Инициализация матрицы после генерации портрета
		void Initialize(int size1, int size2);
		// Умножение матрицы на вектор
		vector <double> operator*(vector<double> x);

		void Clear();

		// Умножение транспонированной матрицы на вектор
		void MultiplyATx(vector<double> a, vector<double>& result);
		// Добавить элемент в матрицу
		void AddElement(int i, int j, double element);
		// Заменить элемент в матрице
		void ChangeElement(int i, int j, double element);
		// Сложение элементов в строке (для проверки)
		void Sum(vector<double>& result);
		//LU-факторизация
		void LU();
		//Вспомогательные функции для решателя
		void LYF(const vector<double>& C, vector<double>& yl);
		void UXY(const vector<double>& C, vector<double>& yu);
	};
}