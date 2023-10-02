#include "BiCGSTAB.h"
#include <math.h>


BiCGSTAB::BiCGSTAB():
	f_size(0),
	f_counter(0)
{}

BiCGSTAB::BiCGSTAB(const CSRMatrix& p_A, const std::vector<double>& p_b):
	f_A(p_A),
	f_b(p_b),
	f_size(p_b.size()),
	f_counter(0)
{}

BiCGSTAB::BiCGSTAB(const std::string & p_pathA, const std::string & p_pathb):
	f_size(0),
	f_counter(0)
{
	setA(p_pathA);
	setb(p_pathb);
}

BiCGSTAB::BiCGSTAB(std::string&& p_pathA, std::string&& p_pathb):
	f_size(0),
	f_counter(0)
{
	setA(p_pathA);
	setb(p_pathb);
}


void BiCGSTAB::setA(std::string&& p_path)
{
	std::vector<double> mat;

	std::vector<int> indexes;
	std::vector<int> amounts;
	std::vector<std::tuple<int, int, double>> tupleVector;



	std::string tmp;

	std::ifstream file_stream(p_path);

	if (!file_stream.is_open())
	{
		throw(std::exception("Matrix A file stream was not opened."));
	}

	while (std::getline(file_stream, tmp))
	{
		std::regex numRegEx("[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+)");

		std::sregex_token_iterator cutIter(
			tmp.cbegin(),
			tmp.cend(),
			numRegEx
		);

		int row = std::stoi(*cutIter);
		cutIter++;

		int col = std::stoi(*cutIter);
		cutIter++;

		double num = std::stod(*cutIter);

		tupleVector.push_back(
			std::make_tuple(row, col, num)
		);
	}

	std::sort(
		tupleVector.begin(),
		tupleVector.end(),
		[](const auto& a, const auto& b)
		{
			if (std::get<0>(a) == std::get<0>(b))
				return std::get<1>(a) < std::get<1>(b);
			return std::get<0>(a) < std::get<0>(b);
		}

	);

	int prev_row = -1;



	for (int i = 0; i < tupleVector.size(); i++)
	{
		mat.push_back(std::get<2>(tupleVector[i]));
		indexes.push_back(std::get<1>(tupleVector[i]) - 1);

		if (std::get<0>(tupleVector[i]) != prev_row)
		{
			amounts.push_back(i);
			prev_row = std::get<0>(tupleVector[i]);
		}
	}

	amounts.push_back(mat.size());

	f_A = CSRMatrix(amounts.size() - 1, mat, indexes, amounts);
	f_size = amounts.size() - 1;

}

void BiCGSTAB::setA(const std::string& p_path)
{
	std::string path = p_path;

	setA(std::move(path));
}

int BiCGSTAB::getSize() const
{
	return f_size;
}

CSRMatrix BiCGSTAB::getA() const
{
	return f_A;
}

void BiCGSTAB::setb(std::string&& p_path)
{
	if (f_size == 0)
		throw(std::exception("Tried to initialize vector b before matrix A."));

	std::vector<double> numbers;

	std::ifstream file_stream(p_path);

	if (!file_stream.is_open())
	{
		throw(std::exception("Vector b file stream was not opened."));
	}

	std::string tmp;

	std::getline(file_stream, tmp);

	std::regex numRegEx("[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+)");

	std::transform(std::sregex_token_iterator{ tmp.cbegin(), tmp.cend(), numRegEx },
				{},
				std::back_inserter(numbers),
				[](const auto& val)
				{
					return std::stod(val.str());
				}
	); 

	if (numbers.size() != f_size)
	{
		throw(std::exception("Matrix A and vector b sizes were incompatible."));
	}

	f_b = numbers;
}

void BiCGSTAB::setb(const std::string& p_path)
{
	std::string path = p_path;

	setb(std::move(path));
}

std::vector<double> BiCGSTAB::getb() const
{
	return f_b;
}

void BiCGSTAB::solve(double p_precision)
{
	
	f_counter = 0;

	f_x = f_A.trace();

	std::vector<double> tech;

	double techNum = 0;

	std::vector<double> r = vectorSubstract(f_b, f_A * f_x);

	std::vector<double> r_0 = r;

	std::vector<double> h;

	std::vector<double> s;

	std::vector<double> t;

	long double alpha = 1;

	long double omega = 1;

	long double rho = 1;

	long double rho_prev = 1;

	long double beta;

	std::vector<double> v(f_size, 0);

	std::vector<double> p(f_size, 0);

	for (int i = 0; i < 999999; i++)
	{
		rho_prev = rho;

		rho = vectorScalarMultpl(r_0, r);

		if (rho == 0)
		{
			throw std::exception("Variable rho turned out to be zero. Program stopped due to possible division by zero.");
		}

		beta = (rho / rho_prev) * (alpha / omega);

		p = specOp2(r, p, v, omega, beta);

		v = f_A * p;

		techNum = vectorScalarMultpl(r_0, v);

		if (techNum == 0)
		{
			throw std::exception("Variable vectorScalarMultpl(r_0, v) turned out to be zero. Program stopped due to possible division by zero.");
		}

		alpha = rho / techNum;

		h = specOp3(f_x, p, alpha, 1);

		s = specOp3(r, v, alpha, -1);

		t = f_A * s;

		techNum = vectorScalarMultpl(t, t);

		if (techNum == 0)
		{
			throw std::exception("Variable vectorScalarMultpl(t, t) turned out to be zero. Program stopped due to possible division by zero.");
		}

		omega = vectorScalarMultpl(t, s) / techNum;

		f_x = specOp3(h, s, omega, 1);

		tech = vectorSubstract(f_b, f_A * f_x);

		if (sqrt(
			vectorScalarMultpl(
				tech,
				tech)) < p_precision)
			break;

		r = specOp3(s, t, omega, -1);


		f_counter++;
	}
	
}

std::vector<double> BiCGSTAB::specOp2(const std::vector<double>& p_r, const std::vector<double>& p_p, const std::vector<double>& p_v, double omega, double beta) const
{
	std::vector<double> res(f_size, 0);
	if (f_size > 1000)
	{
#pragma omp parallel for
		for (int at = 0; at < f_size; at++)
		{
			res[at] = p_r[at] + (p_p[at] - p_v[at] * omega) * beta;
		}
	}
	else
	{
		for (int at = 0; at < f_size; at++)
		{
			res[at] = p_r[at] + (p_p[at] - p_v[at] * omega) * beta;
		}
	}

	return res;
}

std::vector<double> BiCGSTAB::specOp3(const std::vector<double>& p_s, const std::vector<double>& p_t, double num, int sign) const
{
	std::vector<double> res(f_size, 0);
	if (f_size > 1000)
	{
#pragma omp parallel for
		for (int at = 0; at < f_size; at++)
		{
			res[at] = p_s[at] + sign * (p_t[at] * num);
		}
	}
	else
	{
		for (int at = 0; at < f_size; at++)
		{
			res[at] = p_s[at] + sign * (p_t[at] * num);
		}
	}

	return res;
}

std::vector<double> BiCGSTAB::getSolution()
{
	return f_x;
}

unsigned long BiCGSTAB::getIterCount()
{
	return f_counter;
}

double BiCGSTAB::vectorScalarMultpl(const std::vector<double>& p_a, const std::vector<double>& p_b)
{
	double res = 0;
	int size = p_a.size();
	if (f_size > 1000)
	{
//#pragma omp parallel for reduction(+:res)
		for (int i = 0; i < size; i++)
		{
			res += p_a[i] * p_b[i];
		}
	}
	else
	{
		for (int i = 0; i < size; i++)
		{
			res += p_a[i] * p_b[i];
		}
	}
	return res;
}

std::vector<double> BiCGSTAB::vectorSubstract(const std::vector<double>& p_a, const std::vector<double>& p_b)
{
	int size = p_a.size();
	std::vector<double> res(size); 
	
	if (f_size > 1000)
	{
#pragma omp parallel for
		for (int w = 0; w < size; w++)
		{
			res[w] = p_a[w] - p_b[w];
		}
	}
	else
	{
		for (int w = 0; w < size; w++)
		{
			res[w] = p_a[w] - p_b[w];
		}
	}

	return res;
}

std::vector<double> BiCGSTAB::vectorAdd(const std::vector<double>& p_a, const std::vector<double>& p_b)
{

	int size = p_a.size();
	std::vector<double> res(size);
	

	if (f_size > 1000)
	{
#pragma omp parallel for
		for (int w = 0; w < size; w++)
		{
			res[w] = p_a[w] + p_b[w];
		}
	}
	else
	{
		for (int w = 0; w < size; w++)
		{
			res[w] = p_a[w] + p_b[w];
		}
	}

	return res;
}

std::vector<double> BiCGSTAB::vectorMultiply(const std::vector<double>& p_a, double p_omega)
{
	int size = p_a.size();
	std::vector<double> res(size);
	if (f_size > 1000)
	{
#pragma omp parallel for
		for (int w = 0; w < size; w++)
		{
			res[w] = p_a[w] * p_omega;
		}
	}
	else
	{
		for (int w = 0; w < size; w++)
		{
			res[w] = p_a[w] * p_omega;
		}
	}


	return res;

}

double BiCGSTAB::evaluateSolution()
{
	return sqrt(
		vectorScalarMultpl(
			vectorSubstract(f_b, f_A * f_x),
			vectorSubstract(f_b, f_A * f_x)));
}



