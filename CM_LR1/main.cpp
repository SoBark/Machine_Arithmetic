#include "Sparse_Matrix_CSLR.h"
#include "Dense_Matrix.h"
#include <iostream>
#include <cmath>
#include <chrono> // для функций из std::chrono
#define STEPS 330
class Timer
{
private:
	// Псевдонимы типов используются для удобного доступа к вложенным типам
	using clock_t = std::chrono::high_resolution_clock;
	using second_t = std::chrono::duration<double, std::ratio<1> >;

	std::chrono::time_point<clock_t> m_beg;

public:
	Timer() : m_beg(clock_t::now())
	{
	}

	void reset()
	{
		m_beg = clock_t::now();
	}

	double elapsed() const
	{
		return std::chrono::duration_cast<second_t>(clock_t::now() - m_beg).count();
	}
};

double EuclidianNorm(std::vector<double>& Y);

void Special_case_div(Sparse_Matrix_CSLR& Matrix);

void Special_case_mult(Sparse_Matrix_CSLR& Matrix);

void Special_case_div_improved(Sparse_Matrix_CSLR& Matrix);

void Special_case_mult_improved(Sparse_Matrix_CSLR& Matrix);

void txtMatrix2binMatrix(const std::string PATH);

int main() {
	const std::string PATH = "";
	txtMatrix2binMatrix(PATH);
	Sparse_Matrix_CSLR CSLR_Matrix(PATH);
	Dense_Matrix Dense_Matrix;
	CSLR_Matrix.transfer_to_dense(Dense_Matrix.Matrix);		//получение плотной матрицы из CSLR

	////вывод плотной матрицы
	//for (int i = 0; i < Dense_Matrix.Matrix.size(); i++) {
	//	for (int j = 0; j < Dense_Matrix.Matrix.size(); j++)
	//		std::cout << Dense_Matrix.Matrix[i][j] << "\t";
	//	std::cout << "\n";
	//}
	//std::cout << std::endl;

	//std::vector<double> res_CSLR(CSLR_Matrix.N),
	//	res_Dense(Dense_Matrix.Matrix.size()),
	//	X(CSLR_Matrix.N, 1.0);
	////Умножение матрицы на вектор. Время работы
	//Timer T;
	//T.reset();
	//CSLR_Matrix.matrix_mult_vector(X, res_CSLR);
	//std::cout << "CSLR time: " << T.elapsed() << std::endl;
	//T.reset();
	//Dense_Matrix.matrix_mult_vector(X, res_Dense);
	//std::cout << "Dense time: " << T.elapsed() << std::endl;
	////вывод ответов при умножении матриц на вектор
	//for (int i = 0; i < res_CSLR.size(); i++)
	//	std::cout << res_CSLR[i] << " ";
	//std::cout << std::endl;
	//for (int i = 0; i < res_Dense.size(); i++)
	//	std::cout << res_Dense[i] << " ";
	//std::cout << std::endl;

	Special_case_mult(CSLR_Matrix);
	Special_case_mult_improved(CSLR_Matrix);
	Special_case_div(CSLR_Matrix);
	Special_case_div_improved(CSLR_Matrix);
}

//евклидова норма
double EuclidianNorm(std::vector<double>& Y) 
{
	double Norm = 0.0;
	for (int i = 0; i < Y.size(); i++)
		Norm += Y[i] * Y[i];
	return sqrt(Norm);
}

void Special_case_div(Sparse_Matrix_CSLR& Matrix)
{
	int steps = STEPS; //количество итераций
	std::vector<double> X(Matrix.N, 0.0), Y(Matrix.N);
	X[0] = 1.0;
	std::cout << "division\n||Y||\tx1" << std::endl;
	for (int i = 0; i < steps; i++)
	{
		Matrix.matrix_mult_vector(X, Y);
		std::cout << EuclidianNorm(Y) << "\t" << X[0] << std::endl;
		X[0] /= 10; 
	}
	std::cout << std::endl;
}

void Special_case_mult(Sparse_Matrix_CSLR& Matrix)
{
	int steps = STEPS;
	std::vector<double> X(Matrix.N, 0.0), Y(Matrix.N);
	X[0] = 1.0;
	std::cout << "multiplication\n||Y||\tx1" << std::endl;
	for (int i = 0; i < steps; i++)
	{
		Matrix.matrix_mult_vector(X, Y);
		std::cout << EuclidianNorm(Y) << "\t" << X[0] << std::endl;
		X[0] *= 10;
	}
	std::cout << std::endl;
}

void Special_case_mult_improved(Sparse_Matrix_CSLR& Matrix)
{
	int steps = STEPS; //количество итераций
	std::vector<double> X(Matrix.N, 0.0), Y(Matrix.N);
	X[0] = 1.0;
	std::cout << "multiplication\n||Y||\tx1" << std::endl;
	for (int i = 0; i < steps; i++)
	{
		Matrix.matrix_mult_vector(X, Y);
		double Norm = 0.0, nextNorm = 0.0, factor = 1.0;
		for (int j = 0; j < Y.size(); j++) {
			nextNorm = Norm + std::pow(Y[j]*factor, 2.0);
			if (nextNorm >= std::numeric_limits<double>::max())
			{
				//подбор множителя
				while (nextNorm >= std::numeric_limits<double>::max() && factor >= std::numeric_limits<double>::min()) {
					factor /= 10; 
					Norm /= 100; //вынесение множителя
					nextNorm = Norm + std::pow(factor * Y[j], 2.0);
				}
			}
			Norm = nextNorm;
		}
		std::cout << sqrt(Norm)/factor << "\t" << X[0] << std::endl;
		X[0] *= 10;
	}
	std::cout << std::endl;
}

void Special_case_div_improved(Sparse_Matrix_CSLR& Matrix)
{
	int steps = STEPS; //количество итераций
	std::vector<double> X(Matrix.N, 0.0), Y(Matrix.N);
	X[0] = 1.0;
	std::cout << "multiplication\n||Y||\tx1" << std::endl;
	for (int i = 0; i < steps; i++)
	{
		Matrix.matrix_mult_vector(X, Y);
		double Norm = 0.0, nextNorm = 0.0, factor = 1.0;
		for (int j = 0; j < Y.size(); j++) {
			nextNorm = Norm + std::pow(Y[j] * factor, 2.0);
			if (nextNorm <= std::numeric_limits<double>::min())
			{
				//подбор множителя
				while (nextNorm <= std::numeric_limits<double>::min() && factor <= std::numeric_limits<double>::max()) {
					factor *= 10;
					Norm *= 100; //вынесение множителя
					nextNorm = Norm + std::pow(factor * Y[j], 2.0);
				}
			}
			Norm = nextNorm;
		}
		std::cout << sqrt(Norm) / factor << "\t" << X[0] << std::endl;
		X[0] /= 10;
	}
	std::cout << std::endl;
}

void txtMatrix2binMatrix(const std::string PATH)
{
	//размер матрицы
	int N = 0;
	std::ofstream Writer(PATH + "size.bin", std::ios::binary);
	std::ifstream Reader(PATH + "size.txt");
	if (!Reader.is_open())
		throw std::exception("File size.txt was not found...");
	Reader >> N;
	Writer.write((char*)&N, sizeof(N));
	Reader.close();
	Writer.close();
	//di
	Reader.open(PATH + "di.txt");
	Writer.open(PATH + "di.bin", std::ios::binary);
	if (!Reader.is_open())
		throw std::exception("File di.txt was not found...");
	double num_double;
	for (int i = 0; i < N; i++) {
		Reader >> num_double;
		Writer.write((char*)&num_double, sizeof(double));
	}
	Reader.close();
	Writer.close();
	//iptr
	int num_int;
	Reader.open(PATH + "iptr.txt");
	Writer.open(PATH + "iptr.bin", std::ios::binary);
	if (!Reader.is_open())
		throw std::exception("File iptr.txt was not found...");
	for (int i = 0; i < N + 1; i++) {
		Reader >> num_int;
		Writer.write((char*)&num_int, sizeof(int));
	}
	Reader.close();
	Writer.close();
	//jptr
	Reader.open(PATH + "jptr.txt");
	Writer.open(PATH + "jptr.bin", std::ios::binary);
	if (!Reader.is_open())
		throw std::exception("File jptr.txt was not found...");
	int JPTR_SIZE = num_int - 1;
	for (int i = 0; i < JPTR_SIZE; i++) {
		Reader >> num_int;
		Writer.write((char*)&num_int, sizeof(int));
	}
	Reader.close();
	Writer.close();
	//altr
	Reader.open(PATH + "altr.txt");
	Writer.open(PATH + "altr.bin", std::ios::binary);
	if (!Reader.is_open())
		throw std::exception("File altr.txt was not found...");
	for (int i = 0; i < JPTR_SIZE; i++) {
		Reader >> num_double;
		Writer.write((char*)&num_double, sizeof(double));
	}
	Reader.close();
	Writer.close();
	//autr
	Reader.open(PATH + "autr.txt");
	Writer.open(PATH + "autr.bin", std::ios::binary);
	if (!Reader.is_open())
		throw std::exception("File autr.txt was not found...");
	for (int i = 0; i < JPTR_SIZE; i++) {
		Reader >> num_double;
		Writer.write((char*)&num_double, sizeof(double));
	}
	Reader.close();
	Writer.close();
}

 