#include "solver.h"
namespace solver
{
	Solver::Solver()
	{
		qSpline.resize(slae.n);
		qSolution.resize(slae.n);
	}
	void Solver::CalculateLocals(int elementNumber)
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
	void Solver::CalculateA(int elementNumber, const Point * points, int kPoints)
	{
		Element element = grid.elements[elementNumber];

		vector <vector<double>> A;
		A.resize(nSplineFunc);
		for (int i = 0; i < nSplineFunc; i++)
			A[i].resize(nSplineFunc);

		double hx = grid.hx(elementNumber);
		double hy = grid.hy(elementNumber);

		for (int i = 0; i < nSplineFunc; i++) {
			for (int j = i; j < nSplineFunc; j++) {
				// Цикл по точкам
				for (int l = 0; l < kPoints; l++) {
					double p_ksi = points[l].x,
						p_etta = points[l].y;
					A[i][j] += splBasis.phi[i](p_ksi, p_etta, hx, hy) * splBasis.phi[j](p_ksi, p_etta, hx, hy);
				}
			}
		}

		for (int i = 1; i < nSplineFunc; i++)
			for (int j = 0; j < i; j++)
				A[i][j] = A[j][i];

		for (int i = 0; i < nSplineFunc; i++) {
			int id_i = element.dof[i];
			for (int j = 0; j < nSplineFunc; j++) {
				int id_j = element.dof[j];
				slae.A.AddElement(id_i, id_j, A[i][j]);
			}
		}
	}
	void Solver::CalculateAlphaMatrix(int elementNumber)
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
				slae.A.AddElement(id_i, id_j, AM[i][j]);
			}
		}
	}
	void Solver::CalculateBettaMatrix(int elementNumber)
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
				slae.A.AddElement(id_i, id_j, BM[i][j]);
			}
		}
	}
	void Solver::CalculateF(int elementNumber, const Point * points, int kPoints)
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
				double f = get_solution_in_point_bilinear(hx * p_ksi + x0, hy * p_etta + y0, elementNumber, qSolution);
				F[i] += splBasis.phi[i](p_ksi, p_etta, hx, hy) * f;
			}
		}

		for (int i = 0; i < nSplineFunc; i++) {
			int id_i = element.dof[i];
			slae.F[id_i] += F[i];
		}
	}
}
