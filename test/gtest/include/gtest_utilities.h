/*-----------------------------------------------------------------------*\
|    ___                   ____  __  __  ___  _  _______                  |
|   / _ \ _ __   ___ _ __ / ___||  \/  |/ _ \| |/ / ____| _     _         |
|  | | | | '_ \ / _ \ '_ \\___ \| |\/| | | | | ' /|  _| _| |_ _| |_       |
|  | |_| | |_) |  __/ | | |___) | |  | | |_| | . \| |__|_   _|_   _|      |
|   \___/| .__/ \___|_| |_|____/|_|  |_|\___/|_|\_\_____||_|   |_|        |
|        |_|                                                              |
|                                                                         |
|   Author: Alberto Cuoci <alberto.cuoci@polimi.it>                       |
|   CRECK Modeling Group <http://creckmodeling.chem.polimi.it>            |
|   Department of Chemistry, Materials and Chemical Engineering           |
|   Politecnico di Milano                                                 |
|   P.zza Leonardo da Vinci 32, 20133 Milano                              |
|                                                                         |
|-------------------------------------------------------------------------|
|                                                                         |
|   This file is part of OpenSMOKE++.                                     |
|                                                                         |
|   Copyright(C) 2018  Alberto Cuoci                                      |
|   Source-code or binary products cannot be resold or distributed        |
|   Non-commercial use only                                               |
|   Cannot modify source-code for any purpose (cannot create              |
|   derivative works)                                                     |
|                                                                         |
\*-----------------------------------------------------------------------*/

#ifndef _GTEST_UTILITIES_H_

#define _GTEST_UTILITIES_H_

bool CheckErrors(const Eigen::VectorXd& solution, const Eigen::VectorXd& exact_solution,
	const double mean_threshold = 1.e-15, const double max_threshold = 1.e-15)
{
	const unsigned int n = static_cast<unsigned int>(solution.size());

	Eigen::VectorXd e(n);
	for (unsigned int i = 0; i < n; i++)
		e(i) = solution(i) - exact_solution(i);

	const double mean_error = e.norm() / (n*n);
	const double max_error = e.maxCoeff();

	{
		std::ostringstream s_mean_error;  s_mean_error << mean_error;
		std::ostringstream s_max_error;   s_max_error << (max_error);
		const std::string message1 = "Mean error = " + s_mean_error.str();
		const std::string message2 = "Max error  = " + s_max_error.str();

		TEST_COUT << message1.c_str() << std::endl;
		TEST_COUT << message2.c_str() << std::endl;
	}

	if (mean_error < mean_threshold && max_error < max_threshold)
		return true;
	else
		return false;
}

bool CheckErrors(const double solution, const double exact_solution,
	const double threshold = 1.e-15)
{
	if ( std::fabs(solution-exact_solution)/std::fabs(exact_solution) < threshold)
		return true;
	else
		return false;
}

unsigned int MatrixFromFile(const std::string file_name, Eigen::MatrixXd& M)
{
	std::ifstream fInput(file_name.c_str(), std::ios::in);

	unsigned int n = 0;
	unsigned int non_zero = 0;

	fInput >> n;
	fInput >> n;
	fInput >> non_zero;

	// Linear system matrix
	M.resize(n, n);
	M.setConstant(0.);
	for (unsigned int k = 0; k < non_zero; k++)
	{
		unsigned int i, j;
		fInput >> i;
		fInput >> j;
		fInput >> M(i - 1, j - 1);
	}

	fInput.close();

	return n;
}

template<typename Matrix>
unsigned int VectorsFromFile(const std::string file_name, Matrix& M)
{
	std::ifstream fInput(file_name.c_str(), std::ios::in);

	unsigned int n = 0;
	unsigned int m = 0;

	fInput >> n;
	fInput >> m;

	// Linear system matrix
	M.resize(n, m);
	for (unsigned int k = 0; k < n; k++)
		for (unsigned int j = 0; j < m; j++)
			fInput >> M(k, j);

	fInput.close();

	return n*m;
}

template<typename T>
bool compare_vectors(const T& v1, const T& v2)
{
	if (v1.size() != v2.size())
		return false;

	for (unsigned int i = 0; i < v1.size(); i++)
		if (v1[i] != v2[i])
			return false;

	return true;
}

#endif	// _G_TEST_UTILITIES_H
