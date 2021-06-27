#pragma once
#include <string>
#include <fstream>
#include <vector>
class Sparse_Matrix_CSLR
{
public:
	int N;						//������ �������
	std::vector <int>	iptr,	//��������� ������ ������ �����
						jptr;	//��������� ������ ������ ��������

	std::vector <double>	altr,	//������ �����������
							autr,	//������� �����������
							di;		//���������

	Sparse_Matrix_CSLR(const std::string PATH);

	void transfer_to_dense(std::vector<std::vector<double>> &Matrix);

	void matrix_mult_vector(std::vector<double> &X, std::vector<double> & Res);


};

