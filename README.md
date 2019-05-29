# AGH_WIMiIP_Metody_Numeryczne
Laboratoria z Metod Numerycznych na AGH WIMiIP [prowadzący: dr hab. inż. Marcin Hojny]

## Zaimplementowane metody:


## Interpolacja:
Przykładowe dane:

|  Xi |  Yi |
| --- | --- |
| -2  | -1  |
| -1  |  0  |
|  0  |  5  |
|  2  |  99 |
|  4  | -55 |

```cpp
double newPoint = 1; //punkt w którym chcemy wyznaczyć wartość
std::vector<double> xi = { -2.,-1.,0.,2.,4. };
std::vector<double> yi = { -1.,0.,5.,99.,-55. };

```
* Interpolacja Lagrange'a 
```cpp
double numeric::interpolation::lagrange(xi, yi, newPoint);
```
* Interpolacja Newtona
```cpp
double numeric::interpolation::newton(xi, yi, newPoint);
```

## Rozwiązywanie liniowych układów równań:
Przykładowy układ równań:\
2x + 5y + 3z = 5\
4x + 2y + 5z = 4\
3x + 8y + 4z = 9
```cpp
std::vector<std::vector<double>> matrix = { {2,5,3,	5}, {4,2,5,	4}, {3,8,4,	9} };	
```

* Metoda oparta o wzory Cramera
```cpp
std::vector<double> numeric::systemsOfEquasions::polynomial(matrix);
```
* Metoda eliminacji Gaussa
```cpp
std::vector<double> numeric::systemsOfEquasions::gaussianElimination(matrix);
```


## Wyznaczanie miejsca zerowego
Przykładowa funkcja:\
f(x) = x^2 - 1\
f'(x) = 2x\
Obszar poszukiwań miejsca zerowego: < -2 ; 0 >\
Prezyzja: znaleziony punkt różni się o ± 0.01 od rzeczywistego

```cpp
double f(double x){
	return (x*x - 1);
}

double fDerivative(double x){
	return 2 * x;
}
```
* Metoda bisekcji
```cpp
double numeric::zeroOfAFunction::bisection(f, -2., 0., 0.1);
```
* Metoda Newtona-Raphsona
```cpp
double numeric::zeroOfAFunction::newtonRaphson(f, fDerivative, -2., 0., 0.1);
```

## Całkowanie
Przykładowa funkcja:\
1*x^3 + 0*x^2 + 0*x^1 + 2*x^0\
Całkowanie na przedziale < 1 ; 4 >\
Precyzja: ilość iteracji pętli


```cpp
std::vector<double> coefficients = { 1.,0.,0.,2. };

//Funkcja obliczająca wartość wielomianu w punkcie
double f(double x) {
	unsigned int coeffSize = coefficients.size();
	double y = 0.;
	for (unsigned int i = 0; i < coeffSize; i++) {
		y += abs(coefficients[i] * pow(x, coeffSize - i - 1));
	}
	return y;
}

```
* Całkowanie metodą prostokątów
```cpp
double numeric::integration::rectangular(f, 1, 4, 300);
```
* Całkowanie metodą trapezową
```cpp
double numeric::integration::trapezoidal(f, 1, 4, 300);
```
* Całkowanie metodą Simpsona
```cpp
double numeric::integration::simpson(f, 1, 4, 300);
```
* Całkowanie metodą Monte Carlo
```cpp
double numeric::integration::monteCarlo(f, 1, 4, 300);
```
