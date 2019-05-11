#include "pch.h"
#include "numeric.h"

int numeric::interpolation::newtonCacheLength = 0;
std::vector<numeric::interpolation::newtonCacheEntry> numeric::interpolation::newtonCache;

double numeric::interpolation::newton(std::vector<double> x, std::vector<double> y, double newPoint) {
	
	unsigned int size = x.size() - 1;

	double ret = y[0];

	std::vector<double> coefficients;

	for (unsigned int i = 0; i < size; i++) {
		coefficients.push_back(newtonDifferenceQuotient(i+1, x, y, 0));
	}
	
	for (unsigned int i = 0; i < size; i++) {
	
		double temp = 1;
		for (unsigned int j = 0; j <= i; j++) {
			temp *= (newPoint - x[j]);
		}

		ret += coefficients[i] * temp;
	}

	return ret;

}

double numeric::interpolation::newtonDifferenceQuotient(unsigned int size, std::vector<double> x, std::vector<double> y, unsigned int pos) {

	if (size == 1) {
		
		double ret = y[pos+1] - y[pos];
		ret /= (x[pos + 1] - x[pos]);

		return ret;
	}

	for (int i = 0; i < newtonCacheLength; i++) {
		if (newtonCache[i].pos == pos && newtonCache[i].size == size) {
			return newtonCache[i].value;
		}

	}

	double ret = newtonDifferenceQuotient(size - 1, x, y, pos + 1) - newtonDifferenceQuotient(size - 1, x, y, pos);
	ret /= (x[pos+size] - x[pos]);

	newtonCache.push_back({ pos,size,ret });
	newtonCacheLength++;

	return ret;

}


double numeric::systemsOfEquasions::matrixDeterminant(std::vector<std::vector<double>> matrix) {

	unsigned int matrixSize = matrix.size();

	switch (matrixSize) {		
		case 1:
			return matrix[0][0];
		case 2:
			return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
		case 3:
			return matrix[0][0] * matrix[1][1] * matrix[2][2] + matrix[0][1] * matrix[1][2] * matrix[2][0] + matrix[0][2] * matrix[1][0] * matrix[2][1] - (matrix[0][2] * matrix[1][1] * matrix[2][0] + matrix[0][0] * matrix[1][2] * matrix[2][1] + matrix[0][1] * matrix[1][0] * matrix[2][2]);
		case 4:
			return matrix[0][0] * (matrix[1][1] * matrix[2][2] * matrix[3][3] + matrix[1][2] * matrix[2][3] * matrix[3][1] + matrix[1][3] * matrix[2][1] * matrix[3][2] - matrix[1][3] * matrix[2][2] * matrix[3][1] - matrix[1][1] * matrix[2][3] * matrix[3][2] - matrix[1][2] * matrix[2][1] * matrix[3][3]) - matrix[0][1] * (matrix[1][0] * matrix[2][2] * matrix[3][3] + matrix[1][2] * matrix[2][3] * matrix[3][0] + matrix[1][3] * matrix[2][0] * matrix[3][2] - matrix[1][3] * matrix[2][2] * matrix[3][0]- matrix[1][0] * matrix[2][3] * matrix[3][2]-matrix[1][2] * matrix[2][0] * matrix[3][3]) + matrix[0][2] * (matrix[1][0] * matrix[2][1] * matrix[3][3] + matrix[1][1] * matrix[2][3] * matrix[3][0] + matrix[1][3] * matrix[2][0] * matrix[3][1]- matrix[1][3] * matrix[2][1] * matrix[3][0]-matrix[1][0] * matrix[2][3] * matrix[3][1] - matrix[1][1] * matrix[2][0] * matrix[3][3]) - matrix[0][3] * (matrix[1][0] * matrix[2][1] * matrix[3][2] + matrix[1][1] * matrix[2][2] * matrix[3][0] + matrix[1][2] * matrix[2][0] * matrix[3][1] - matrix[1][2] * matrix[2][1] * matrix[3][0]- matrix[1][0] * matrix[2][2] * matrix[3][1] - matrix[1][1] * matrix[2][0] * matrix[3][2]);
		default:
			return INFINITY;
	}

}

std::vector<double> numeric::systemsOfEquasions::polynomial(std::vector<std::vector<double>> matrix) {

	std::vector<std::vector<double>> wx;
	std::vector<std::vector<double>> wy;
	std::vector<std::vector<double>> wz;

	unsigned int matrixSize = matrix.size();

	for (unsigned int i = 0; i < matrixSize; i++) {

		std::vector<double> wxTemp;
		std::vector<double> wyTemp;
		std::vector<double> wzTemp;

		for (unsigned int j = 0; j < matrixSize; j++) {


			if (j == 0) {
				wxTemp.push_back(matrix[i][matrixSize]);
			}
			else {
				wxTemp.push_back(matrix[i][j]);
			}


			if (j == 1) {
				wyTemp.push_back(matrix[i][matrixSize]);
			}
			else {
				wyTemp.push_back(matrix[i][j]);
			}


			if (j == 2) {
				wzTemp.push_back(matrix[i][matrixSize]);
			}
			else {
				wzTemp.push_back(matrix[i][j]);
			}

		}

		wx.push_back(wxTemp);
		wy.push_back(wyTemp);
		wz.push_back(wzTemp);
	}

	double detW = matrixDeterminant(matrix);
	double detWx = matrixDeterminant(wx);
	double detWy = matrixDeterminant(wy);
	double detWz = matrixDeterminant(wz);


	if (detW == 0) {

		if (detWx == 0. || detWy == 0. || detWz == 0.) {
			std::cerr << "No solution!\n";
			return {};
		}

		if (detWx == 0. && detWy == 0. && detWz == 0.) {
			std::cerr << "Infinitely many solutions!\n";
			return {};
		}

	}

	std::vector<double> ret;
	ret.reserve(3);

	ret.push_back(detWx / detW);
	ret.push_back(detWy / detW);
	ret.push_back(detWz / detW);
	
	return ret;
}


double numeric::interpolation::lagrange(std::vector<double> x, std::vector<double> fx, double newPoint) {

	unsigned int size = x.size();

	std::vector<double> li;
	li.reserve(size);

	for (unsigned int i = 0; i < size; i++) {

		double l = 1;

		for (unsigned int j = 0; j < size; j++) {

			if (j == i) {
				continue;
			}

			l*= static_cast<double>((newPoint - x[j])) / (x[i] - x[j]);

		}

		li.push_back(l);
	}

	double ret = 0.;

	for (unsigned int i = 0; i < size; i++) {
	
		ret += fx[i] * li[i];
	
	}

	return ret;

}

double numeric::integration::rectangular(double (*function)(double), int start, int end, unsigned int precision) {

	double dx = end;
	dx -= start;
	dx /= precision;

	double outcome = 0.;	
	double x = start;

	for (unsigned int i = 0; i < precision; i++) {

		x += dx;
		outcome += function(x);
	
	}
	 
	return dx * outcome;

}

double numeric::integration::trapezoidal(double (*function)(double), int start, int end, unsigned int precision){

	double dx = end;
	dx -= start;
	dx /= precision;

	double outcome = 0.;
	double x = start;

	for (unsigned int i = 0; i < precision; i++) {

		double numerator;

		//X_i
		numerator = function(x);

		x += dx;

		//X_i+1
		numerator += function(x);


		// outcome += (Xi + Xi+1) / 2
		outcome += numerator / 2;

	}

	return dx * outcome;
		
}

double numeric::integration::simpson(double(*function)(double), double start, double end, unsigned int steps) {

	double avgLen = (end - start) / steps;
	std::vector<double> startPos;
	
	double temp = start;
	for (unsigned int i = 0; i <= steps; i++) {
		startPos.push_back(temp);
		temp += avgLen;
	}

	double outcome = 0.;
	for (unsigned int i = 0; i < steps; i++) {
		
		unsigned int startIndex = i;
		unsigned int middleIndex = ++i;
		unsigned int endIndex = ++i;

		double h = startPos[endIndex] - startPos[middleIndex];

		double temp = function(startPos[startIndex]) + 4 * function(startPos[middleIndex]) + function(startPos[endIndex]);
		temp *= (h / 3);
			
		outcome += temp;
			
		i--;
	}
	   
	return outcome;
}

double numeric::integration::monteCarlo(double(*function)(double), double start, double end, unsigned int precision) {

	srand(static_cast<unsigned int>(time(NULL)));

	double avg = 0.;

	for (unsigned int i = 0; i < precision; i++) {
	
		double f = static_cast<double>(rand()) / RAND_MAX;
		double x = start + f * (end - start);

		avg += function(x);
		
	}

	avg /= precision;


	return avg * abs(end - start);

}

double numeric::integration::polynomialFunction(const std::vector<double>& coefficients, double argument) {
	
	unsigned int coeffSize = coefficients.size();

	double y = 0.;
	for (unsigned int i = 0; i < coeffSize; i++) {

		y += abs(coefficients[i] * pow(argument, coeffSize - i - 1));
			
	}

	return y;
}

double numeric::zeroOfAFunction::bisection(double (*f)(double), double a, double b, double precision){

	
	if (f(a) * f(b) >= 0) {
		//No solution between said points
		return INFINITY;
	}

	double middle = a + b;
	middle /= 2;

	if (f(middle) == 0) {
		return middle;
	}

	while (abs(a - b) >= precision) {
	
		if (f(middle) * f(a) < 0) {
			b = middle;
		}
	
		if (f(middle) * f(b) < 0) {
			a = middle;
		}
	
		middle = a + b;
		middle /= 2;
	
	}

	return middle;
}

double numeric::zeroOfAFunction::newtonRaphson(double (*f)(double), double (*fDerrivative)(double), double a, double b, double precision) {

	if (f(a) * f(b) >= 0) {
		//No solution between said points
		return INFINITY;
	}

	double middle = a + b;
	middle /= 2;

	if (f(middle) == 0) {
		return middle;
	}

	double xn = a;
		
	while (abs(f(xn)) >= precision) {
		xn = xn - (f(xn) / fDerrivative(xn));
	}

	return xn;
}

std::vector<double> numeric::systemsOfEquasions::gaussianElimination(std::vector<std::vector<double>> matrix) {

	unsigned int matrixSize = matrix.size();
	unsigned int matrixCols = matrix[0].size();
	unsigned int startPos = 0;
	std::vector<double> coefficients;
	coefficients.reserve(matrixSize);

	//Triangulate matrix
	while (coefficients.size() < 1) {
	
		for (unsigned int i = startPos+1; i < matrixSize; i++) {
			double coeff = matrix[i][startPos] / matrix[startPos][startPos];
			for (unsigned int j = startPos; j < matrix[i].size(); j++) {
				matrix[i][j] -= matrix[startPos][j] * coeff;
			}
		}

		//If last row has only 1 coefficient, push value to vector thus end loop
		if (matrix[matrixSize-1][matrixCols-3] == 0) {
			double xn = matrix[matrixSize - 1][matrixCols - 1] / matrix[matrixSize - 1][matrixCols - 2];
			coefficients.push_back(xn);
		}
		startPos++;
	}

	//Start searching for other coefficients
	unsigned int coefficientsSize = coefficients.size();
	while (coefficientsSize < matrixSize) {
		double xn = matrix[matrixSize - coefficientsSize - 1][matrixCols - 1];
		for (unsigned int i = 0; i < coefficientsSize; i++) {
			xn -= matrix[matrixSize - coefficientsSize - 1][matrixCols - 2 - i] * coefficients[i];
		}
		xn /= matrix[matrixSize - coefficientsSize - 1][matrixCols - coefficientsSize - 2];
		coefficients.push_back(xn);
		coefficientsSize = coefficients.size();
	}

	//Reverse vector
	std::vector<double> ret;
	for (unsigned int i = 0; i < matrixSize; i++) {
		ret.push_back(coefficients[matrixSize - 1 - i]);
	}

	return ret;
}
