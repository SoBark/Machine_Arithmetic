#include "Sparse_Matrix_CSLR.h"

Sparse_Matrix_CSLR::Sparse_Matrix_CSLR(const std::string PATH)
{
	//размер матрицы
	std::ifstream Reader(PATH + "size.bin", std::ios::binary);
	if (!Reader.is_open())
		throw std::exception("File size.bin was not found...");
	Reader.read((char*)&N, sizeof(N));
	Reader.close();
	//di
	Reader.open(PATH + "di.bin", std::ios::binary);
	if (!Reader.is_open())
		throw std::exception("File di.bin was not found...");
	di.resize(N);
	for (int i = 0; i < N; i++)
		Reader.read((char*)&di[i], sizeof(double));
	Reader.close();
	//iptr
	Reader.open(PATH + "iptr.bin", std::ios::binary);
	if (!Reader.is_open())
		throw std::exception("File iptr.bin was not found...");
	iptr.resize(N + 1);
	for (int i = 0; i < N + 1; i++)
		Reader.read((char*)&iptr[i], sizeof(int));
	Reader.close();
	//jptr
	Reader.open(PATH + "jptr.bin", std::ios::binary);
	if (!Reader.is_open())
		throw std::exception("File jptr.bin was not found...");
	int JPTR_SIZE = iptr[N] - 1;
	jptr.resize(JPTR_SIZE);
	for (int i = 0; i < JPTR_SIZE; i++)
		Reader.read((char*)&jptr[i], sizeof(int));
	Reader.close();
	//altr
	Reader.open(PATH + "altr.bin", std::ios::binary);
	if (!Reader.is_open())
		throw std::exception("File altr.bin was not found...");
	altr.resize(JPTR_SIZE);
	for (int i = 0; i < JPTR_SIZE; i++)
		Reader.read((char*)&altr[i], sizeof(double));
	Reader.close();
	//autr
	Reader.open(PATH + "autr.bin", std::ios::binary);
	if (!Reader.is_open())
		throw std::exception("File autr.bin was not found...");
	autr.resize(JPTR_SIZE);
	for (int i = 0; i < JPTR_SIZE; i++)
		Reader.read((char*)&autr[i], sizeof(double));
	Reader.close();
}

void Sparse_Matrix_CSLR::transfer_to_dense(std::vector<std::vector<double>>& Matrix)
{
	Matrix.resize(N);
	for (int i = 0; i < N; i++)
		Matrix[i].resize(N, 0.0);
	//di
	for (int i = 0; i < N; i++)
		Matrix[i][i] = di[i];
	//Заполнение верхнего и нижнего треугольников
	for (int i = 0; i < N; i++)
		for (int j = iptr[i] - 1; j < iptr[i + 1] - 1; j++)
		{
			Matrix[i][jptr[j] - 1] = altr[j];
			Matrix[jptr[j] - 1][i] = autr[j];
		} 
}

void Sparse_Matrix_CSLR::matrix_mult_vector(std::vector<double> & X, std::vector<double>& Res)
{
	//инициализация результата через умножения вектора на диагональ
	for (int i = 0; i < N; i++) Res[i] = X[i] * di[i];
	//проход по всем строкам и столбцам с учётом формата
	for (int i = 0; i < N; i++)
		for (int j = iptr[i] - 1; j < iptr[i + 1] - 1; j++)
		{
			Res[i] += X[jptr[j] - 1] * altr[j];
			Res[jptr[j] - 1] += X[i] * autr[j];
		}
}
