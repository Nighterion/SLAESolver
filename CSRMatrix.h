#pragma once

extern int threads;

#include <math.h>
#include "IncompatibleDimException.h"
#include <iostream>
#include <vector>
#include <chrono>
#include <random>
#include <omp.h>

class CSRMatrix
{
private:

	std::vector<double> f_Mat;
	std::vector<int> f_indexes;
	std::vector<int> f_amounts;

	int f_size;

private:

public:

	CSRMatrix();

	CSRMatrix(unsigned int _size, const std::vector<double>& p_Mat, const std::vector<int>& p_idx, const std::vector<int>& p_amts);

	CSRMatrix(unsigned int p_size, int matsize, int idxsize, int amountssize);

	CSRMatrix(const CSRMatrix& _M);

	~CSRMatrix()
	{
		f_Mat.clear();
		f_indexes.clear();
		f_amounts.clear();
	}

	CSRMatrix operator* (double _d_to_mult) const;

	CSRMatrix& operator= (const CSRMatrix& _M);

	friend std::ostream& operator<<(std::ostream& p_out, const CSRMatrix& _M);

	friend CSRMatrix operator* (const double _d_to_mult, CSRMatrix& _M);

	std::vector<double> operator* (const std::vector<double>& p_vec);

	int Get_size() const;

	int Get_amount() const;

	CSRMatrix transpose() const;

	std::vector<double> trace();

	std::vector<double> getMat() const
	{
		return f_Mat;
	}

	std::vector<int> getIndexes() const
	{
		return f_indexes;
	}

	std::vector<int> getAmounts() const
	{
		return f_amounts;
	}

};