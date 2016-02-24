#include "Gauss.h"
#include <cstring>
#include <omp.h>

Gauss::Gauss()
{
    A = NULL;
}

Gauss::~Gauss()
{
    delete_matrix();
}

void Gauss::convert_matrix()
{
    n_gauss = n;

    delete_matrix();
    A = new complex<double> * [n_gauss];
    for(size_t i = 0; i < n_gauss; i++)
        A[i] = new complex<double> [n_gauss + 1];

    for(size_t i = 0; i < n_gauss; i++)
    {
        for(size_t k = gi[i]; k < gi[i + 1]; k++)
        {
            A[i][gj[k]] = gl[k];
            A[gj[k]][i] = gu[k];
        }
        A[i][i] = di[i];
        A[i][n_gauss] = rp[i];
    }
}

void Gauss::delete_matrix()
{
    if(A)
    {
        for(size_t i = 0; i < n_gauss; i++)
            delete [] A[i];
        delete [] A;
        A = NULL;
    }
}

// Инициализация несимметричная
void Gauss::init(size_t * gi_s, size_t * gj_s, complex<double> * di_s,
                 complex<double> * gl_s, complex<double> * gu_s, size_t n_s)
{
    n = n_s;
    gi = gi_s;
    gj = gj_s;
    di = di_s;
    gl = gl_s;
    gu = gu_s;
}

// Инициализация симметричная
void Gauss::init(size_t * gi_s, size_t * gj_s, complex<double> * di_s,
                 complex<double> * gg_s, size_t n_s)
{
    init(gi_s, gj_s, di_s, gg_s, gg_s, n_s);
}

// Получение решения
void Gauss::solve(complex<double> * solution, complex<double> * rp_s, double eps, size_t max_iter)
{
    rp = rp_s;
    (void)eps;
    (void)max_iter;
    convert_matrix();

//    long N = n_gauss;
//    //верхний треугольный вид
//    for(long i = 0; i < N; i++)
//    {
//        for(long j = N; j >= i; j--)
//            A[i][j] = A[i][j] / A[i][i];
//        for(long j = i + 1; j < N; j++)
//            for(long k = N; k >= i; k--)
//                A[j][k] -= A[i][k] * A[j][i];
//    }
//    //диагональный вид
//    for(long i = N - 1; i > 0; i--)
//        for(long j = i - 1; j >= 0; j--)
//            A[j][N] -= A[j][i] * A[i][N];

    int N = n_gauss;
    int M = N;

    // Заказываем количество параллельных процессов
    omp_set_num_threads(4);
    // Задаем параллельный блок из (4-х) процессов
#pragma omp parallel
    {
        int size = omp_get_num_threads();
        int Np = M / size;
        int rank = omp_get_thread_num();
        int uk = Np * rank;
        // Прямой ход
        for(int k = 0; k < M; k++)
        {
            // На каждом шаге цикла вычисляем полосы строк (dn-начало, dk-конец полосы),
            // обрабатываемых каждым процессом
            int dn = (uk * (k < uk)) + ((k + 1) * (k >= uk));
            int dk = uk + Np - 1;
            // Синронизация процессов
#pragma omp barrier
            // Деление на коэффициент элементов текущей строки процессом master
#pragma omp master
            {
                complex<double> MAD = 1.0 / A[k][k];
                for(int j = M; j >= k; j--)
                    A[k][j] *= MAD;
            }
            // Синронизация процессов для того чтобы master успел обработать сроку
#pragma omp barrier
            // Обработка строк параллельными процессами
            for(int i = dn; i <= dk; i++)
            {
                for(int j = M; j >= k; j--)
                    A[i][j] -= A[i][k] * A[k][j];
            }
        }
        // Обратный ход
        uk = (M - 1 - Np * (size - rank - 1));
        for(int k = M - 1; k >= 0; k--)
        {
            // На каждом шаге цикла вычисляем полосы строк (dn-начало, dk-конец полосы),
            // обрабатываемых каждым процессом
            int dn = ((uk * (k > uk)) + ((k-1) * (k <= uk)));
            int dk = uk - Np + 1;
            // Синронизация процессов, т.к. получение корней – процесс последовательный
#pragma omp barrier
            // Обработка строк параллельными процессами
            for(int i = dn; i >= dk; i--)
                A[i][M] -= A[k][M] * A[i][k];
        }
    } // Конец параллельного блока

    for(size_t i = 0; i < n_gauss; i++)
        solution[i] = A[i][n_gauss];

    delete_matrix();
}
