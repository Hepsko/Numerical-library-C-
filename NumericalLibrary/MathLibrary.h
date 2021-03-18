#pragma once
#include <string>

namespace MathLibrary
{
	class Numeric
	{
	public:
		static double* readPoints(std::string filename);

		static double readDataAndClculateLagrange(std::string filename, double point);

		static double lagrange(double* X, double* Y, int n, double x);

		static bool unique(double X[], int length);

		static double newton(double* tabX, double* tabY, int row, int i);

		static double calculateNewton(double* tabX, double* tabY, int point, int num);
	};
}