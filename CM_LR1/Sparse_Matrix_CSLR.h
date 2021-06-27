#pragma once
#include <string>
#include <fstream>
#include <vector>
class Sparse_Matrix_CSLR
{
public:
	int N;						//размер матрицы
	std::vector <int>	iptr,	//индексный массив начала строк
						jptr;	//индексный массив начала столбцов

	std::vector <double>	altr,	//нижний треугольник
							autr,	//верхний треугольник
							di;		//диагональ

	Sparse_Matrix_CSLR(const std::string PATH);

	void transfer_to_dense(std::vector<std::vector<double>> &Matrix);

	void matrix_mult_vector(std::vector<double> &X, std::vector<double> & Res);


};

