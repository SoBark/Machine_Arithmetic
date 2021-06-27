#include "Dense_Matrix.h"

void Dense_Matrix::matrix_mult_vector(std::vector<double>& X, std::vector<double>& Res)
{
	for (int i = 0; i < Res.size(); i++) Res[i] = 0;
	for (int i = 0; i < Matrix.size(); i++)
		for (int j = 0; j < Matrix.size(); j++)
			Res[i] += Matrix[i][j] * X[j];			
}
