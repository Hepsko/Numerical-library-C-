#include "MathLibrary.h"
#include <string>
#include <fstream>
#include <algorithm>
#include <iostream>
namespace MathLibrary
{



	double Numeric::valueAtPoint(double fun[], double  point,  size_t size) {

		double y = 0.;
		for (unsigned int i = 0; i < size; i++) {
			y += fun[i] * pow(point, size - 1 - i);
		}
		return y;
	}

	double Numeric::rectangleIntegration(double fun_form[], double  start, double end, int prec,  size_t size) {
		double dx = (end - start) / prec;
		double output = 0.;

		for (unsigned int i = 0; i < prec; i++) {

			start += dx;
			output += valueAtPoint(fun_form, start,  size);

		}

		return dx * output;
	}

	double	Numeric::trapezeIntegration(double fun_form[], double start, double end, int prec,  size_t size){
		double dx = (end - start) / prec;
		double output = 0.;
		
		

		for (unsigned int i = 0; i < prec; i++) {

			double numerator;

			numerator = valueAtPoint(fun_form, start, size);
			start += dx;
			numerator += valueAtPoint(fun_form, start, size);
			output += numerator / 2;

		}
		return dx * output;
}

	double Numeric::newton(double* tabX, double* tabY, double newPoint, int num) {
		double result = tabY[0];  
		double tmp = 0;			
		for (int i = 0; i < num; i++)
		{
			tmp = 0;				
			for (int j = 0; j < i; j++) {	
				if (j == 0){
					tmp = newPoint - tabX[0];
				}
				else
				{
					tmp *= newPoint - tabX[j];
				}	
			}
			result += tmp * Numeric::newtonDifferenceQuotient(tabX, tabY, i, 0);
		}
		return result;
	}


	double Numeric::newtonDifferenceQuotient(double* tabX, double* tabY, int row, int i) {
		if (row == 0) {
			return 1;
		}
		if (row == 1) {
			return  (tabY[i + 1] - tabY[i]) / (tabX[i + 1] - tabX[i]);
		}
		return	(newtonDifferenceQuotient(tabX, tabY, row - 1, i + 1) - newtonDifferenceQuotient(tabX, tabY, row - 1, i)) / (tabX[i + row] - tabX[i]);
	}


	
	
	double Numeric::lagrange(double* X, double* Y,  int n, double x) {
			double y = 0;
		
			for (int i = 0; i < n; i++)
			{
				double tmp = 1;
				for (int j = 0; j < n; j++)
				{
					if (i != j)
					{
						 tmp *= ((x - X[j]) / (X[i] - X[j]));
					}
				}
				y += tmp * Y[i];
			}
			return y;
		}

	bool Numeric::unique(double X[], int length) {
		for (int i = 0; i < length; i++) {
			for (int j = 0; j < length; j++) {
				if (X[i] == X[j] && i != j) {return false;}
			}
		}
		return true;
	}

	double Numeric::readDataAndClculateLagrange(std::string filename, double point) {
		/*
			File pattern:
			X Y
			X Y
		*/
		int data = 0;
		std::string line;
		std::fstream handle(filename, std::ios::in);
		if (handle.good() == true) {
			//counting input data
			while (getline(handle, line)) { data++; }

			if (data <= 1) {
				handle.close();
				throw "Error -> Not enough data";
			}

			//handle point reset
			handle.clear();
			handle.seekg(0);

			//data transfer
			double* X = new double[data];
			double* Y = new double[data];
			for (int i = 0; i < data; i++) {
				handle >> X[i];
				handle >> Y[i];
			}
			handle.close();

			//check unique of x
			if (!unique(X, data)) { throw " Error -> x point should be unique"; }

			double* max = std::max_element(X, X + data);
			double* min = std::min_element(X, X + data);
			if (*max < point || *min > point) { throw " Error -> Point is out of range"; }

			//execute lagrange form MathLibrary
			return lagrange(X, Y, data, point);
		}
		else {
			//Throw err when there is a problem with reading
			handle.close();
			throw "Error -> problem with opening data";
		}
	}
}