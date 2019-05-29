#include "pch.h"
#include <iostream>
#include "numeric.h"

int main() {

	std::vector<double> xi = { 0.,5.,5.,0. };
	std::vector<double> yi = { 0.,0.,5.,5. };
	std::vector<double> weight = { 1.,1. };
	std::vector<double> point = { -0.5773502692, 0.5773502692 };
	std::cout << "Area: " << numeric::integration::gaussianQuadrature(xi, yi, weight, point);
		
}