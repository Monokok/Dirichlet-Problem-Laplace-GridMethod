#include <Windows.h>
#include <mpi.h>
#include <stdio.h>
#include <sstream>
#include <fstream>

#include <iostream>

using namespace std;

int parallelSolution(int argc, char** argv)
{
    //шапка
    SetConsoleOutputCP(1251);
    int size, rank;
    if (MPI_Init(&argc, &argv) != MPI_SUCCESS)
        return 1;
    if (MPI_Comm_size(MPI_COMM_WORLD, &size) != MPI_SUCCESS)
    {
        MPI_Finalize();
        return 2;
    }
    if (MPI_Comm_rank(MPI_COMM_WORLD, &rank) != MPI_SUCCESS)
    {
        MPI_Finalize();
        return 3;
    }
    MPI_Status status;
    int msgtag = 1;

    //кол-во итераций
    int Kmax = 1000;

    double L = 10.0, //длина стороны квадратной области
        h = 0.01; //шаг сетки
    int N, n, m;

    //–азмерность сетки
    N = (int)(L / h) + 1;
    m = (int)(L / h) + 1;

    //ƒеление сетки между процессами
    int l1 = N / size; //сколько строк (по оси X) получает каждый процесс минимум.
    int l2 = N % size; //остаток, который нужно распределить (если N не делитс€ на size).

    int* kol = new int[size]; //kol[i] хранит, сколько строк у каждого процесса.
    for (int i = 0; i < size; i++)
        kol[i] = l1;

    //раскидываем остаток
    if (l2)
    {
        if (rank == size - 1)
            l1++;   //последний +1
        l2--;       // уменьшаем остаток
        kol[size - 1]++;    // фиксируем, что последнему досталось больше
    }
    if (l2) // если остаток всЄ ещЄ осталс€
    {
        if (rank < l2) // раздаЄм оставшиес€ строки первым процессам
            l1++;
        for (int i = 0; i < rank; i++)
            kol[i]++;
    }
    
    //”чЄт "граничных строк" (обмен между сосед€ми)
    n = l1;
    if ((!rank) || (rank == size - 1)) //0-й и последний
        n++;
    else
        n += 2; //не крайние процессы


    //u[n][m] Ч локальный кусок сетки.
    double** u = new double* [n];
    for (int i = 0; i < n; i++)
    {
        u[i] = new double[m];
    }
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            u[i][j] = 0.0;


    //количество дл€ каждого
    int l = 0;
    double x = 0.0, y = 0.0;
    if (rank)
    {
        int sum = 0;
        for (int i = 0; i < rank; i++)
            sum += kol[i];
        l = sum;
        x = h * (sum - 1); //x с которого начинаем
    }

    //”становка граничных условий
    for (int i = 0; i < n; i++)
    {
        u[i][0] = 30.0 * x * (100.0 - x);//ad лева€ граница
        u[i][m - 1] = 20.0 * sqrt(x);//bc права€ граница
        x += h;
    }
    for (int j = 0; j < m; j++)
    {
        if (rank == 0)
            u[0][j] = 30.0 * (1 - y); //ab верхн€€ граница
        if (rank == size - 1)
            u[n - 1][j] = 20.0 * y;//cd нижн€€ граница
        y += h;
    }
    double tn = MPI_Wtime();

    for (int k = 0; k < Kmax; k++)
    {
        // аждый процесс обмениваетс€ строками с сосед€ми :
        //if (rank > 0)
        if (rank != 0) //все кто rank > 0 принимают от левого
            MPI_Recv(u[0], m, MPI_DOUBLE, rank - 1, msgtag, MPI_COMM_WORLD, &status);

        //if (k > 0 && rank < size - 1) //все кто > 0 и < последнего принимают от правого
        if (k != 0 && rank != size - 1) // со второй итерации мы уже принимаем от правых соседей (кроме последнего)
            MPI_Recv(u[n - 1], m, MPI_DOUBLE, rank + 1, msgtag, MPI_COMM_WORLD, &status);

        //просчет в зоне ответственности
        for (int i = 1; i < n - 1; i++)
        {
            for (int j = 1; j < m - 1; j++)
            {
                u[i][j] = 0.25 * (u[i - 1][j] + u[i + 1][j] + u[i][j - 1] + u[i][j + 1]);
            }
            if (i == 2  //означает, что мы только что посчитали строку u[1] (второе значение = первое в зоне ответственности процесса). верхн€€ внутренн€€ строка текущего процесса.
                    && 
                rank != 0) //rank > 0 - отдаем Ћ≈¬ќћ”. не надо делать отправку у самого первого процесса (у него нет соседа слева).
                MPI_Send(u[1], m, MPI_DOUBLE, rank - 1, msgtag, MPI_COMM_WORLD); //ѕосле того как процесс обновил вторую строку, он сразу отправл€ет еЄ вверх (в соседний процесс rank-1).
        }

        //if (rank < size - 1)
        if (rank != size - 1) //не последние отправл€ют правому rank < last
            MPI_Send(u[n - 2], m, MPI_DOUBLE, rank + 1, msgtag, MPI_COMM_WORLD);
    }

    //врем€
    double tk = MPI_Wtime();

    //вывод
    if (!rank)
    {
        char fn[256];
        sprintf_s(fn, "resultTemp.txt");
        std::ofstream result(fn);

        if (!result.is_open()) { perror("Error"); exit(1); }

        result << "Time = " << tk - tn << endl;

        for (int i = 0; i < 100; i++)
        {
            for (int j = 0; j < 100; j++)
            {
                result << "u[" << i << "][" << j << "] = " << u[i][j] << "\n";
            }
            result << endl;
        }
        result.close();
    }

    MPI_Finalize();
}