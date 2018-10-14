#include "integration.h"
using namespace matrix;
using namespace basis;
using namespace integration;

namespace slae
{
	struct SLAE
	{
		//Размерность задачи
		int n;
		//Глобальная матрица
		Matrix A;
		//Вектор приближенного решения
		vector <double> u;
		//Глобальный вектор правой части
		vector <double> F;
		//Локальные матрицы
		//Матрица жёсткости
		array<array<double, 4>, 4> G;
		//Матрица массы
		array<array<double, 4>, 4> M;
		//Локальный вектор правой части
		array <double, 4> locF;

		double dphix(int i, int elementNumber, double ksi, double etta);
		double dphiy(int i, int elementNumber, double ksi, double etta);
		double GetAbsJ(int elementNumber, double ksi, double etta);
		double GetAbsBQJ(int elementNumber, double ksi, double etta);

		//Сборка локальных матриц жёсткости
		void CalculateG(int elementNumber);
		//Сборка локальных матриц масс
		void CalculateM(int elementNumber);
		//Сборка локальных правых частей
		void CalculateLocalF(int elementNumber);
		//Сборка локальных матриц(векторов) и добавление в глобальные
		void CalculateLocals(int elementNumber);

		//Вектор праввой части для первого краевого 
		array<double, 2> g;
		//Нахождение правой части для 1ого краевого условия
		void Calculate_g(int formNumber, int orientation, int elNumber);
		//Вычисление 1ого краевого условия для одного узла
		void CalculateBoundaries1ForNode(int node, double gi, double weight);
		//Учёт первого краевого условия
		void CalculateBoundaries1(int number);

		
		//Вектор невязки
		vector <double> r;
		//Вектор спуска
		vector <double> z;

		//Вычисление нормы вектора
		double Norm(const vector<double>& x);
		//Скалярное произведение векторов
		double Scalar(const vector<double>& x, const vector<double>& y);

		//Генерация СЛАУ на i-ой итерации по времени
		void GenerateSLAE();
		
		double Rel_Discrepancy();
		//Решатель ЛОС с LU-факторизацией
		void LULOS();
		SLAE();
		void Solve();
		~SLAE() {};
	};
}