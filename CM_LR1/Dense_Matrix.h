#pragma once
#include <vector>
class Dense_Matrix
{
public:
	std::vector <std::vector<double>> Matrix;

	void matrix_mult_vector(std::vector<double>& X, std::vector<double> &Res);
};

