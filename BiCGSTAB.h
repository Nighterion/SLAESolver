#pragma once
#include "CSRMatrix.h"
#include <regex>
#include <fstream>
#include <tuple>


class BiCGSTAB
{
private:

	CSRMatrix f_A;
	std::vector<double> f_b;
	std::vector<double> f_x;
	int f_size;
	unsigned long f_counter;

public:

	BiCGSTAB();

	BiCGSTAB(const CSRMatrix& p_A, const std::vector<double>& p_b);

	BiCGSTAB(const std::string& p_pathA, const std::string& p_pathb);

	BiCGSTAB(std::string&& p_pathA, std::string&& p_pathb);

	void setA(std::string&& p_path);

	void setA(const std::string& p_path);

	int getSize() const;

	CSRMatrix getA() const;

	void setb(std::string&& p_path);

	void setb(const std::string& p_path);

	std::vector<double> getb() const;

	void solve(double p_precision);

	std::vector<double> specOp2(const std::vector<double>& p_r, const std::vector<double>& p_p, const std::vector<double>& p_v, double omega, double beta) const;

	std::vector<double> specOp3(const std::vector<double>& p_s, const std::vector<double>& p_t, double num, int sign) const;
	
	std::vector<double> getSolution();

	unsigned long getIterCount();

	double vectorScalarMultpl(const std::vector<double>&, const std::vector<double>&);

	std::vector<double> vectorSubstract(const std::vector<double>&, const std::vector<double>&);

	std::vector<double> vectorAdd(const std::vector<double>&, const std::vector<double>&);

	std::vector<double> vectorMultiply(const std::vector<double>&, double p_omega);

	double evaluateSolution();
};