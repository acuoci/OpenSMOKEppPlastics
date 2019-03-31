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
|   Copyright(C) 2019  Alberto Cuoci                                      |
|   Source-code or binary products cannot be resold or distributed        |
|   Non-commercial use only                                               |
|   Cannot modify source-code for any purpose (cannot create              |
|   derivative works)                                                     |
|                                                                         |
\*-----------------------------------------------------------------------*/

#include "Utilities.h"

namespace opensmokepp::plastics
{

	int SearchForLC(std::vector<double>& list_boiling_temperature, const double T)
	{
		for (int i = 0; i < list_boiling_temperature.size(); i++)
			if (list_boiling_temperature[i] > T)
				return i + 1;

		std::cout << "Fatal error: Maximum allowed temperature is " << list_boiling_temperature[list_boiling_temperature.size() - 1] << std::endl;
		std::cout << "Press enter to exit...";
		getchar();
		exit(-1);
	}

	double ConversionActivationEnergy(const double E, const std::string units)
	{
		if (units == "cal/mol")
			return E;
		else if (units == "J/mol")
			return E * 4.186;
		else
		{
			std::cout << "Fatal error: wrong units for activation energy. Available units: cal/mol | J/mol" << std::endl;
			std::cout << "Press enter to exit..." << std::endl;
			getchar();
			exit(-1);
		}

		return 0;
	}

	double ConversionFrequencyFactor(const double A, const std::string units)
	{
		return A;
	}

	double SchultzFloryDistribution(const double gamma, const double p, const int i)
	{
		return gamma * std::pow(p, i - 1)*std::pow(1. - p, 2.);
	}

	double SchultzDistribution(const double gamma, const double p, const int i)
	{
		return gamma / p * (4.*i) / std::pow(p - 1., 2.)*std::pow((p - 1.) / (p + 1.), i + 1);
	}


	void InitialDistribution(	const PolymerDistribution distribution,
								const double MW_monomer, const double MWE, 
								const int lumping_start, const int lumping_step,
								const double epsilon, const int stretching_coefficient,
								Eigen::VectorXd& y, int& N,
								const std::string folder_name)
	{

		// Maximum number of monomeric units to be tracked
		const int nm = stretching_coefficient * static_cast<int>(MWE / MW_monomer);

		// Normalization factor
		const double mass_initial = 100.;

		// Allocate memory for PDF
		y.resize(nm);
		y.setZero();

		int nsom = 0;
		if (distribution == SCHULTZ || distribution == SCHULTZ_FLORY)
		{
			const double chi = MWE / MW_monomer;	// averaged degree of polymerization
			const double p = 1. - 1. / chi;
			const double a = 1. - p;				// fraction of unreacted monomer

			double sum = 0.;
			for (int i = 1; i <= nm; i++)
			{
				if (distribution == SCHULTZ_FLORY)
					y(i - 1) = SchultzFloryDistribution(mass_initial / MW_monomer, p, i);
				else if (distribution == SCHULTZ)
					y(i - 1) = SchultzDistribution(mass_initial / MW_monomer, chi, i);

				sum += y(i - 1)*MW_monomer*i;

				if (std::fabs(100. - sum) < epsilon)
				{
					nsom = lumping_start + lumping_step * (static_cast<int>((i - lumping_start) / lumping_step) + 1);

					for (int j = i - 1; j <= nsom; j++)
					{
						if (distribution == SCHULTZ_FLORY)
							y(j - 1) = SchultzFloryDistribution(mass_initial / MW_monomer, p, j);
						else if (distribution == SCHULTZ)
							y(j - 1) = SchultzDistribution(mass_initial / MW_monomer, chi, j);
					}

					break;
				}
			}

			if (distribution == SCHULTZ_FLORY)
			{
				std::cout << "Residual fraction:        " << a << std::endl;
				std::cout << "Degree of polymerization: " << chi << std::endl;
				std::cout << "Mean:                     " << 2. / a - 1. << std::endl;
				std::cout << "Variance:                 " << (2. - 2.*a) / a / a << std::endl;
				std::cout << "Std deviation:            " << std::sqrt((2. - 2.*a) / a / a) << std::endl;
				std::cout << "Total number of chains:   " << nsom << std::endl;
			}
			else if (distribution == SCHULTZ_FLORY)
			{
			}
		}
		else if (distribution == DIRAC_DELTA)
		{
			const int lmax = static_cast<int>(MWE / MW_monomer);
			y(lmax - 1) = mass_initial / (MW_monomer*lmax + 2);
			nsom = lumping_start + lumping_step * (static_cast<int>((lmax - lumping_start) / lumping_step) + 1);
		}

		// Normalization
		{
			double sum = 0.;
			for (int i = 1; i <= nsom; i++)
				sum += y(i - 1)*(MW_monomer*i + 2.);

			for (int i = 1; i <= nsom; i++)
				y(i - 1) *= mass_initial / sum;
		}

		// Total number of equations (i.e. chains to be tracked)
		N = nsom + 1;

		// Write on file
		bool is_write_on_file = true;
		if (is_write_on_file == true)
		{
			const std::string file_name = folder_name + "/" + "Distribution.out";

			std::ofstream fOut(file_name.c_str(), std::ios::out);
			fOut.setf(std::ios::scientific);

			fOut << std::setw(8)  << std::left << "Units";
			fOut << std::setw(16) << std::left << "PDF(mass)";
			fOut << std::setw(16) << std::left << "CDF(mass)";
			fOut << std::endl;

			double total_mass = 0.;
			for (int i = 1; i <= nsom; i++)
			{
				total_mass += y(i - 1)*(MW_monomer*i + 2.);
				fOut << std::setw(8)  << std::left << i;
				fOut << std::setw(16) << std::left << y(i - 1);
				fOut << std::setw(16) << std::left << total_mass;
				fOut << std::endl;
			}

			fOut.close();
		}

	}
	

	void LumpingSetup(	const bool is_lumping, const int lumping_start, const int lumping_step,
						const double MW_monomer, Eigen::VectorXd& y, int& N, const std::string folder_name)
	{
		if (is_lumping == true)
		{
			if (lumping_start >= N)
			{
				std::cout << "Lumping is impossible since the maximum number of units to be tracked is too small" << std::endl;
				std::cout << "Press enter to exit..." << std::endl;
				getchar();
				exit(-1);
			}

			int j = lumping_start;
			int j1 = j;
			int j2 = j + 1;
			int yk = lumping_step;
			int flag = 0;
			double yl = 0.;
			double ym = 0.;

			for (;;)
			{
				yl = (MW_monomer*j + 2)*y(j - 1)*yk / lumping_step + yl;
				ym = (MW_monomer*j + 2)*y(j - 1)*(lumping_step - yk) / lumping_step + ym;

				if (j == N)
					break;

				j++;

				if (flag == 0)
					yk--;
				else if (flag == 1)
					yk++;

				if (yk == 0)
				{
					const int py = lumping_start + (j1 - lumping_start)*lumping_step;
					y(j1 - 1) = yl / (py*MW_monomer + 2.);
					j1 += 2;
					yl = 0.;
					flag = 1;
				}
				else if (yk == lumping_step)
				{
					const int py = lumping_start + (j2 - lumping_start)*lumping_step;
					y(j2 - 1) = ym / (py*MW_monomer + 2.);
					j2 += 2;
					ym = 0.;
					flag = 0;
				}
			}

			int nl = 0;

			if (flag == 1)
			{
				const int py = lumping_start + (j2 - lumping_start)*lumping_step;
				y(j2 - 1) = ym / (py * MW_monomer + 2.);
				nl = j2;
			}
			else if (flag == 0)
			{
				const int py = lumping_start + (j1 - lumping_start)*lumping_step;
				y(j1 - 1) = yl / (py * MW_monomer + 2.);
				nl = j1;
			}

			for (int kl = nl + 1; kl <= N; kl++)
				y(kl - 1) = 0.;

			// Update the maximum number of units to be tracked
			N = nl + 1;

			// Write on file
			{
				const std::string file_name = folder_name + "/" + "Lumping.out";

				std::ofstream fOut(file_name.c_str(), std::ios::out);
				fOut.setf(std::ios::scientific);

				fOut << std::setw(8) << std::left  << "Index";
				fOut << std::setw(8) << std::left  << "Units";
				fOut << std::setw(16) << std::left << "PDF(mass)";
				fOut << std::setw(16) << std::left << "CDF(mass)";
				fOut << std::endl;

				double cdf = 0;
				for (int jj = 1; jj <= N - 1; jj++)
				{
					int jf = jj;
					if (jj > lumping_start)
						jf = lumping_start + lumping_step * (jj - lumping_start);

					cdf += y(jj - 1)*(jf * MW_monomer + 2.);

					fOut << std::setw(8)  << std::left << jj;
					fOut << std::setw(8)  << std::left << jf;
					fOut << std::setw(16) << std::left << y(jj - 1);
					fOut << std::setw(16) << std::left << cdf;
					fOut << std::endl;
				}

				fOut.close();
			}
		}

		
		// Reshaping the distribution
		{
			// Total number of equations (including a dummy final equation)
			const int neq = 5 * N + 1;

			Eigen::VectorXd tmp = y;
			y.resize(neq);
			y.setZero();
			for (int i = 0; i < std::min(neq, static_cast<int>(tmp.size())); i++)
				y(i) = tmp(i);
		}
	}


	void InitialDistribution(Eigen::VectorXd& y, const double epsi, const double MWm, const double MWp,
		int& N, double& wg, const std::string folder_name)
	{
		y.setZero();

		const int ls = static_cast<int>((y.size() - 1) / 3.);

		const double p = 1. - MWm / MWp;

		const int nn = static_cast<int>(1. / (1. - p));


		double sum = 0.;
		for (int i = 1; i <= 1000000; i++)
		{
			y(i - 1) = wg / MWm * std::pow(p, i - 1)*std::pow(1. - p, 2.);

			sum += y(i - 1)*(MWm*i + 2.);

			if (std::fabs(wg - sum) / wg < epsi)
			{
				N = i;
				break;

			}

			if (i == ls)
			{
				std::cout << "Initial distribution: Max dimension error" << std::endl;
				std::cout << "Press enter to exit..." << std::endl;
				getchar();
				exit(-1);
			}

			N = i;
		}

		const int ne = N * 3 + 1;

		Eigen::VectorXd ytmp = y;
		y.resize(ne);
		for (unsigned int i = 0; i < std::min(ytmp.size(), y.size()); i++)
			y(i) = ytmp(i);

		ytmp = y;
		y.resize(ne + 33);
		y.setZero();
		for (unsigned int i = 0; i < ytmp.size(); i++)
			y(i) = ytmp(i);


		// ndim è la dimensione del vettore inizializzato a seconda del polimero
		// Fornisce la prima distribuzione asimmetrica del polimero nel file.ris

		const std::string file_name = folder_name + "/" + "Distribution.out";

		std::ofstream fDistribution(file_name.c_str(), std::ios::out);
		fDistribution.setf(std::ios::scientific);

		double t = 0.;
		for (int i = 1; i <= N; i++)
		{
			const double s = y(i - 1)*(i*MWm + 2.);
			t += s;

			fDistribution << i * MWm << " " << s << " " << t << " " << y(i) << std::endl;
		}

		wg = t;	// TODO

		fDistribution.close();

		int NPA = N;
		int NOL = N;
		int NDO = N;

		y(NPA + NOL - 1) = 0.;
		y(ne - 1 - 1) = 0.;
		y(ne - 2 - 1) = 0.;
	}
}