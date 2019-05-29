#pragma once
#include <vector>
#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <time.h>


class numeric {

public:

	class interpolation {

	public:

		class newtonCacheEntry {
		public:
			unsigned int size;
			unsigned int pos;

			double value;
		};

		static std::vector<newtonCacheEntry> newtonCache;
		static int newtonCacheLength;
		
		static double lagrange(const std::vector<double> & x, const std::vector<double> & y, double newPoint);

		static double newton(const std::vector<double> & x, const std::vector<double> & y, double newPoint);
		static double newtonDifferenceQuotient(unsigned int size, const std::vector<double> & x, const std::vector<double> & y, unsigned int pos);

	
	};

	class integration {

	public:
		static double rectangular(double (*function)(double), int start, int end, unsigned int precision);
		static double trapezoidal(double (*function)(double), int start, int end, unsigned int precision);

		static double simpson(double (*function)(double), double start, double end, unsigned int steps);

		static double monteCarlo(double (*function)(double), double start, double end, unsigned int precision);
		
		static double gaussianQuadrature(const std::vector<double> & xi, const std::vector<double> & yi, const std::vector<double> & weight, const std::vector<double> & point);

	};
	
	class systemsOfEquasions {
	public:

		static double matrixDeterminant(std::vector<std::vector<double>> * matrixPtr);
		static std::vector<double> polynomial(std::vector<std::vector<double>> * matrixPtr);

		static std::vector<double> gaussianElimination(std::vector<std::vector<double>> * matrixPtr);
	};

	class zeroOfAFunction{
	public:
		static double bisection(double (*f)(double), double a, double b, double precision);
		static double newtonRaphson(double (*f)(double), double (*fDerrivative)(double), double a, double b, double precision);
		
	};


};

