#include <stdio.h>
#include <omp.h>
#include <math.h>

const double L = 10.0, //длина стороны квадратной области
	h = 0.01; //шаг сетки
const int kmax = 1000; //максимальное количество итераций

void sequentialSolution()
{
	//Размерность сетки
	int n = (int)(L / (double)h) + 1, 
		m = (int)(L / (double)h) + 1; 
	
	//динамический двумерный массив
	double** U = new double* [n]; //хранит значение функции в узле сетки.
	for (int i = 0; i < n; i++)
		U[i] = new double[m];

	//По умолчанию нули
	for (int i = 0; i < n; i++)
		for (int j = 0; j < m; j++)
			U[i][j] = 0.0;

	double x = 0.0, y = 0.0; //Граничные условия
	//AD и BC задаем значения функции на левой и правой сторонах квадрата
	for (int i = 0; i < n; i++)
	{
		U[i][0] = 30.0 * x * (100.0 - x);	// левая граница (AD) 30x(100 - x)
		U[i][m - 1] = 20.0 * sqrt(x);		// правая граница (BC) 20Sqrt(x)
		x += h;
	} 
	//AB и CD Верхняя и нижняя границы:
	for (int j = 0; j < m; j++)
	{
		U[0][j] = 30.0 * (1 - y);	// верхняя граница (AB) 30(1-y)
		U[n - 1][j] = 20.0 * y;		// нижняя граница (CD) 20y
		y += h;
	}

	//
	double tn = omp_get_wtime();
	double k = 0;
	while (k != kmax)
	{
		for (int i = 1; i < n - 1; i++)
		{
			for (int j = 1; j < m - 1; j++)
			{
				//вычисляется новое значение U[i][j] как среднее четырёх соседей:
				U[i][j] = 0.25 * (U[i - 1][j] + U[i + 1][j] +
					U[i][j - 1] + U[i][j + 1]); 
			}
		}

		k++;
	}
	double tk = omp_get_wtime();
	double deltat = tk - tn;

	//Вывод первых 100x100 значений
	for (int i = 0; i < 100; i++)
	{
		for (int j = 0; j < 100; j++)
		{
			printf("U[%d][%d] = %lf\n", i, j, U[i][j]);
		}
	}
	printf("Time = %lf\n", deltat);

	//Освобождение памяти
	for (int i = 0; i < n; i++)
	{
		delete[](U[i]);
	}
	delete[](U);

}
