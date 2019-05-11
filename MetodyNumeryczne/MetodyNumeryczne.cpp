#include "pch.h"
#include <iostream>
#include "numeric.h"

std::vector<double> coefficients = { 1.,0.,0.,2. };

double f(double x) {
	unsigned int coeffSize = coefficients.size();
	double y = 0.;
	for (unsigned int i = 0; i < coeffSize; i++) {
		y += abs(coefficients[i] * pow(x, coeffSize - i - 1));
	}
	return y;
}

int main() {

	
	//std::vector<std::vector<double>> matrix = { {4,-2,4,-2,		8}, {3,1,4,2,	7}, {2,4,2,1,	10}, {2,-2,4,2,		2} };
	//std::vector<std::vector<double>> matrix = { {10,-7,0,	6}, {-3,2,6,	4}, {5,-1,5,	3} };
	std::vector<std::vector<double>> matrix = { {2,5,3,	5}, {4,2,5,	4}, {3,8,4,	9} };
	std::vector<double> ret = numeric::systemsOfEquasions::gaussianElimination(matrix);

	std::vector<double> ret2 = numeric::systemsOfEquasions::polynomial(matrix);


	for (unsigned int i = 0; i < ret.size(); i++) {
		std::cout << "x" << i + 1 << " = " << ret[i] << "\n";
	}

	for (unsigned int i = 0; i < ret2.size(); i++) {
		std::cout << "x" << i + 1 << " = " << ret2[i] << "\n";
	}

	

		
}