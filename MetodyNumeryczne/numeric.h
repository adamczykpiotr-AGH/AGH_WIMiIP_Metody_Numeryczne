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
		
		static double lagrange(std::vector<double> x, std::vector<double> fx, double newPoint);

		static double newton(std::vector<double> x, std::vector<double> fx, double newPoint);
		static double newtonDifferenceQuotient(unsigned int size, std::vector<double> x, std::vector<double> y, unsigned int pos);

	
	};

	/*DOT¥D OK*/


	class integration {

	public:
		static double rectangular(double (*function)(double), int start, int end, unsigned int precision);
		static double trapezoidal(double (*function)(double), int start, int end, unsigned int precision);

	/*+*/	static double simpson(double (*function)(double), double start, double end, unsigned int steps);

	/*+*/	static double monteCarlo(double (*function)(double), double start, double end, unsigned int precision);
		
		static inline double polynomialFunction(const std::vector<double> & coefficients, double argument);

	};
	
	class systemsOfEquasions {
	public:

		static double matrixDeterminant(std::vector < std::vector<double>> matrix);
		static std::vector<double> polynomial(std::vector<std::vector<double>> matrix);

		static std::vector<double> gaussianElimination(std::vector<std::vector<double>> matrix);
	};

	/*PONI¯EJ GIT*/
	class zeroOfAFunction{
	public:
		static double bisection(double (*f)(double), double a, double b, double precision);
		static double newtonRaphson(double (*f)(double), double (*fDerrivative)(double), double a, double b, double precision);
		
	};


};

