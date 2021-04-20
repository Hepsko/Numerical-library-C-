#pragma once
#include <string>

namespace MathLibrary
{
	class Numeric
	{
	public:
		
		static double lagrange(double* X, double* Y, int n, double x);

		static bool unique(double X[], int length);

		static double newtonDifferenceQuotient(double* tabX, double* tabY, int row, int i);

		static double newton(double* tabX, double* tabY, double newPoint, int num);

		static double rectangleIntegration(double fun_form[], double start, double end, unsigned int prec,  size_t size);

		static double trapezeIntegration(double fun_form[], double start, double end, unsigned int prec,  size_t size);

		static double valueAtPoint(double fun[], double  point, size_t size);

		static double simpsonIntegration(double fun_form[], double start, double end, unsigned int steps, size_t size);

		static double monteCarloIntegration(double fun_form[], double start, double end, unsigned int prec, size_t size);

		static double* gaussianElimination(double **equation,  size_t size);

		static double** readTheEquation(std::string filename);

		static double* readPoints(std::string filename);

		static double readDataAndClculateLagrange(std::string filename, double point);
	};
}