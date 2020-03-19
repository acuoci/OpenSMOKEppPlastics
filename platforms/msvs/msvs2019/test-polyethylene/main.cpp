// Eigen C++ library
#include <Eigen/Dense>
#include <fstream>
#include <iostream>
#include <iomanip>

#include "PolyethyleneKinetics.h"

// Folder containing data
std::string test_folder = "../../../../data/";

const double Rgas = 1.987;

void FatalErrorMessage(const std::string message)
{
	std::cout << message << std::endl;
	std::cout << "Press eneter to exit..." << std::endl;
	getchar();
	exit(-1);
}

void InitialDistribution(Eigen::VectorXd& y, const double epsi, const double MWm, const double MWp,
	int& N, double& wg);

void Corrections(const Eigen::VectorXd& Y, Eigen::VectorXd& CORRO, Eigen::VectorXd& CORRP,
	const int NPA, const int LC, const double CA, const double CD);

void Aggregation(const Eigen::VectorXd& Y, const int J,
	double& SOM1, double& SOM2, double& SOM3, double& SOM4, double& SOM5, double& SOM6, double& SOM7,
	Eigen::VectorXd& SY, Eigen::VectorXd& SY1, Eigen::VectorXd& SY12, Eigen::VectorXd& SY2,
	const int NPA, const int LC);

int main()
{
	//bool seq_code = false;

	//if (seq_code == false)
	{
		// Initial distribution from input file
		Eigen::VectorXd y(600001);
		const double epsi = 0.01;
		const double MW_monomer = 14.;
		const double MW_polymer = 5000.;

		int N = 0;
		double wg = 1.;
		InitialDistribution(y, epsi, MW_monomer, MW_polymer, N, wg);

		// TODO
		/*
		for (unsigned int i = 1; i <= 7030; i++)
			y(i - 1) = 1.e-5*(1. + 0.75*std::sin(2.*3.14*i / 10.));

		{
			std::ifstream fInput("InputDistribution.inp", std::ios::in);

			for (;;)
			{
				int index;
				double value;
				fInput >> index;
				if (index == 0)
					break;
				fInput >> value;

				y(index - 1) = value;

				fInput >> value;
			}

			fInput.close();
		}
		*/
		// Define polystyrene kinetics
		opensmokepp::plastics::PolyethyleneKinetics kinetics;
		//ptkinetics = &kinetics;
		kinetics.SetVerbose(false);
		kinetics.SetMaxNumberOfUnits(N);
		kinetics.SetMolecularWeightPolymer(MW_polymer);
		kinetics.SetMolecularWeightMonomer(MW_monomer);
		kinetics.SetInitialAccelerationFactor(true);


		const double T = 600. + 273.15;
		const double P = 1.;

		const int LC = kinetics.MinNumberOfUnits(T, P);
		const double Teb = kinetics.BoilingTemperature(LC,P);

		Eigen::VectorXd x = y / y.sum();
		const double Ctot = kinetics.LiquidDensity(T)/kinetics.MolecularWeightMonomer();	// (kmol/m3 or mol/l)
		Eigen::VectorXd c = Ctot * x;	// (kmol/m3 or mol/l)

		kinetics.SetMinNumberOfUnits(LC);
		kinetics.SetStatus(T, P, c);

		

		Eigen::VectorXd n = y;

		const double nG = kinetics.SumGas(n);							// mol
		const double nL = kinetics.SumLiquid(n);						// mol
		const double mG = kinetics.SumGasMW(n);							// g
		const double mL = kinetics.SumLiquidMW(n);						// g
		const double VL = (mL / 1000.) / kinetics.LiquidDensity(T);		// m3

		kinetics.UpdateInitialAccelerationCoefficient(VL*1000., wg);
		kinetics.KineticConstants();
		kinetics.FormationRates();

		std::cout << wg << " " << mL << " " << mG + mL << std::endl;

		std::cout << 0 << " " << 0 << " " << T << " " << mG << " " << mL << " " << mG + mL << std::endl;

		const double mL0 = kinetics.SumLiquidMW(n);						// g
		const double VL0 = (mL0 / 1000.) / kinetics.LiquidDensity(T);		// m3

		bool integrate = true;
		if (integrate == true)
		{
			std::ofstream fOut("History.out", std::ios::out);
			fOut.setf(std::ios::scientific);

			
			double t = 0.;
			const double dt = 1e-4;
			const double tEnd = 1.;
			const int nsteps = static_cast<int>(tEnd / dt);
			for (int k = 0; k <= nsteps; k++)
			{
				t += dt;

				const double nG = kinetics.SumGas(n);							// mol
				const double nL = kinetics.SumLiquid(n);						// mol
				const double mG = kinetics.SumGasMW(n);							// g
				const double mL = kinetics.SumLiquidMW(n);						// g
				const double VL = (mL / 1000.) / kinetics.LiquidDensity(T);		// m3

				x = n / n.sum();		// mol fractions
				c = n / (VL*1000.);		// concentrations (in kmol/m3 or mol/l)

				//T = temperature_profile.Interpolate(t);
				kinetics.SetStatus(T, P, c);
				kinetics.UpdateInitialAccelerationCoefficient(VL*1000., wg);
				kinetics.KineticConstants();
				kinetics.FormationRates();
				
				n += dt * kinetics.R()*(VL*1000.);			// moles
				for (int i = 0; i < n.size(); i++)
					if (n(i) < 0.) n(i) = 0.;

				if (k % 10 == 1)
					std::cout << k << " " << t << " " << T << " " << mL/wg << " " << mG << " " << mL << " " << mG + mL << std::endl;

				if (k % 10 == 1)
					fOut << k << " " << t << " " << T << " " << mL / wg << " " << mG << " " << mL << " " << mG + mL << std::endl;
			}

			fOut.close();
		}
	}


	//if (seq_code == true)
	{
	Eigen::VectorXd a(10);
	Eigen::VectorXd aa(10);
	Eigen::VectorXd an(10);
	Eigen::VectorXd e(10);

	a(0) = 7.94328;	aa(0) = 14.; an(0) = 0.0; e(0) = 78000.0;	// initiation 1
	a(1) = 4.00;	aa(1) = 13.; an(1) = 0.0; e(1) = 68000.0;	// initiation 2
	a(2) = 1.00;	aa(2) = 10.2; an(2) = 1.0; e(2) = 14500.0;	// termination
	a(3) = 1.50;	aa(3) = 8.5; an(3) = 0.0; e(3) = 12000.0;	// abstraction 1
	a(4) = 1.50;	aa(4) = 8.5; an(4) = 0.0; e(4) = 13500.0;	// abstraction 2
	a(5) = 1.00;	aa(5) = 14.1; an(5) = 0.0; e(5) = 30500.0;	// beta-decomposition
	a(6) = 1.00;	aa(6) = 11.0; an(6) = 0.0; e(6) = 20600.0;	// bck 1-4
	a(7) = 1.60;	aa(7) = 10.0; an(7) = 0.0; e(7) = 14500.0;	// bck 1-5
	a(8) = 5.00;	aa(8) = 9.0; an(8) = 0.0; e(8) = 14500.0;	// bck 1-6
	a(9) = 1.00;	aa(9) = 0.0; an(9) = 0.0; e(9) = 180000.0;	// ??? TODO

	int ne = 600001;
	const int nre = 6;

	const double epsi = 0.01;
	const double MWm = 14.;
	double wg = 1.;
	const double rho = 760;
	const double MWp = 5000.0;

	Eigen::VectorXd y(ne);
	int N;
	InitialDistribution(y, epsi, MWm, MWp, N, wg);
	
	int NPA = N;
	int NOL = N;
	int NDO = N;

	// TODO
	for (unsigned int i = 1; i <= 7030; i++)
		y(i - 1) = 1.e-5*(1. + 0.75*std::sin(2.*3.14*i / 10.));

	{
		std::ifstream fInput("InputDistribution.inp", std::ios::in);

		for (;;)
		{
			int index;
			double value;
			fInput >> index;
			if (index == 0)
				break;
			fInput >> value;

			y(index - 1) = value;

			fInput >> value;
		}

		fInput.close();
	}

	const double T = 400. + 273.15;
	const double P = 1.;

	const int LC = static_cast<int>(
			5.46e-5*std::pow(T, 2.)*std::pow(1.-std::log(P)/10.5, 2.));

	const double Teb = std::sqrt((LC+1.)/5.46e-5) *1./(1.-std::log(P)/10.5);


	const double alfa = 1.7;
	const double beta = 100000.0;
	const double gamma = 0.91;
	const double fi = 0.016;
	double PP = 0.;

	//System
	Eigen::VectorXd DY(y.size());
	DY.setZero();
	{
		if (PP <= 783)
			PP = static_cast<int>(24. + 12.*(783. - T) / 110.);
		else
			PP = 24;

	//	const double FN = 950. * std::sqrt(PP);

		Eigen::VectorXd k_(10);
		k_(0) = a(0) * std::pow(10., aa(0))*std::exp(-e(0) / Rgas / T);
		k_(1) = a(1) * std::pow(10., aa(1))*std::exp(-e(1) / Rgas / T);
		k_(2) = std::pow(10., 12.8)*T/400. * 14.*PP / rho * std::exp(-1258.*std::sqrt(PP) / Rgas / T)*fi;
		k_(3) = a(3) * std::pow(10., aa(3))*std::exp(-e(3) / Rgas / T);
		k_(4) = a(4) * std::pow(10., aa(4))*std::exp(-e(4) / Rgas / T);
		k_(5) = a(5) * std::pow(10., aa(5))*std::exp(-e(5) / Rgas / T);
		k_(6) = a(6) * std::pow(10., aa(6))*std::exp(-e(6) / Rgas / T);
		k_(7) = a(7) * std::pow(10., aa(7))*std::exp(-e(7) / Rgas / T);
		k_(8) = a(8) * std::pow(10., aa(8))*std::exp(-e(8) / Rgas / T);
		k_(9) = a(9) * std::pow(10., aa(9))*std::exp(-e(9) / Rgas / T);

		const double CD = 2.;
		const double CA = 2.;

		const double AKI = k_(5)*k_(3) / (k_(4)*rho/14.+k_(5));
		const double AKT = 1.0 * std::pow(10, 10.)*std::exp(-6000. / Rgas / T)*T / 400.;

		double SOMMAR1 = 0.;
		double SOMMAR2 = 0.;
		double SOMMAR3 = 0.;


		for (int i = NPA; i >= LC; i--)
		{
			SOMMAR1 += (i - 3)*y(i - 1) + (i - 5)*y(i + NPA - 1) + (i - 7)*y(i + NPA + NOL - 1);
			SOMMAR3 += i*(y(i-1) + y(i + NPA-1) + y(i + NPA + NOL-1));
		}

		for (int i = NPA + NOL; i >= NPA + LC; i--)
			SOMMAR2 += y(i - 1) + 2.*y(i + NOL - 1);

		const double CSI = SOMMAR3 * 14. / wg;

		const double ART = (0.9 - CSI) / 0.05;
		
		double COEFF = 1.;
		if (ART < -2.)
			COEFF = std::sqrt(40.);
		else if (ART >= -2. && ART < 2.)
			COEFF = std::sqrt(1. + 19.5*(1. - std::tanh(ART)));

		k_(0) *= COEFF;

		const double AKS = (k_(0)*SOMMAR1 + k_(1)*SOMMAR2) / SOMMAR3;

		const double RAD = std::sqrt(2.*rho*AKS / 14. / AKT);

		const double AKP = RAD * AKI;


		double SOM1 = 0.;
		double SOM2 = 0.;
		double SOM3 = 0.;
		double SOM4 = 0.;
		double SOM5 = 0.;
		double SOM6 = 0.;
		double SOM7 = 0.;
		double SJ1 = 0.;
		double SJ2 = 0.;

		Eigen::VectorXd SY(22);		SY.setZero();
		Eigen::VectorXd SY1(33);	SY1.setZero();
		Eigen::VectorXd SY12(11);	SY12.setZero();
		Eigen::VectorXd SY2(22);	SY2.setZero();

		double FR = 0.5;

		const double AKD = 0.5*AKP*FR*CA;
		const double AKA = 0.5*AKP*CA;

		

		Eigen::VectorXd CORRO(20);	CORRO.setZero();
		Eigen::VectorXd CORRP(20);	CORRP.setZero();

		Corrections(y, CORRO, CORRP, NPA, LC, CA, CD);
		//WRITE(108, 500) Y(NE), (CORRO(I), I = 1, 20);
		const double sum1 = CORRO.sum();
		const double sum2 = CORRP.sum();

		Eigen::VectorXd AKES(2);
		for (int i = 1; i <= 2; i++)
			AKES(i-1) = a(i + 3 - 1)*std::pow(10., aa(i + 3 - 1))*std::exp(-e(i + 3 - 1) / (Rgas*T));

		Eigen::VectorXd AKIS(3);
		for (int i = 1; i <= 3; i++)
			AKIS(i - 1) = a(i + 6 - 1)*std::pow(10., aa(i + 6 - 1)) * std::exp(-e(i + 6 - 1) / (Rgas*T));

		Eigen::VectorXd alpha(3);
		Eigen::VectorXd beta(3);
		for (int i = 1; i <= 3; i++)
		{
			alpha(i - 1) = AKIS(i-1) / (AKIS(i-1) + AKES(1-1)*(rho / 14.));
			beta(i - 1) = AKIS(i-1) / (AKIS(i-1) + AKES(2-1)*(rho / 14.));
		}

		Eigen::VectorXd ABB(11);
		ABB(0) = alpha(0);
		ABB(1) = alpha(1);
		ABB(2) = alpha(2);
		ABB(3) = alpha(0)*beta(0);
		ABB(4) = alpha(0)*beta(1) + alpha(1)*beta(0);
		ABB(5) = alpha(0)*beta(2) + alpha(2)*beta(0) + alpha(1)*beta(1);
		ABB(6) = alpha(1)*beta(2) + alpha(2)*beta(1) + 2.*alpha(0)*std::pow(beta(0),2.);
		ABB(7) = alpha(2)*beta(2) + alpha(0)*beta(0)*beta(1) +
				 alpha(0)*beta(1)*beta(0) + 2.*alpha(1)*std::pow(beta(0), 2.);
		ABB(8) = alpha(0)*beta(0)*beta(2) + alpha(0)*beta(2)*beta(0) + 2.*alpha(2) *
				 std::pow(beta(0),2.) + alpha(1)*beta(1)*beta(0) + 
				 2.*alpha(0)*std::pow(beta(1), 2.) + alpha(1)*beta(0)*beta(1);
		ABB(9) = 6.*alpha(0)*std::pow(beta(0), 3.) + alpha(0)*beta(1)*beta(2) +
			      alpha(0)*beta(2)*beta(1) + alpha(1)*beta(0)*beta(2) +
			      alpha(1)*beta(2)*beta(0) + alpha(2)*beta(1)*beta(0) +
			      alpha(2)*beta(0)*beta(1) + 2.*alpha(1)*std::pow(beta(1), 2.);
		ABB(10) = 2.*alpha(0)*std::pow(beta(0),2.) * beta(1) +
			      2.*alpha(0)*beta(0)*beta(1)*beta(0) +
			      2.*alpha(0)*beta(1)*std::pow(alpha(0), 2.) +
			      6.*alpha(1)*std::pow(beta(0),3.) + 2.*alpha(2)*std::pow(beta(1),2.) +
			      alpha(1)*beta(2)*beta(1) + alpha(1)*beta(1)*beta(2) +
			      alpha(2)*beta(2)*beta(0) + alpha(2)*beta(0)*beta(2) +
			      2.*alpha(0)*std::pow(beta(2),2.);

		double ATOT = 0.0;
		for (int i = 1; i <= 11; i++)
			ATOT += ABB(i - 1);


		Eigen::VectorXi M(17);	M.setZero();
		Eigen::VectorXi LW(16); LW.setZero();
		Eigen::VectorXi N(19);	N.setZero();

		M(1-1) = 1;
		N(1-1) = 1;
		LW(1-1) = 1;

		int IV = 1;
		const int III = 1;
		const int II = 1;

		double TOT = 0.;
		for (int i = LC; i <= NPA; i++)
			TOT += y(NPA + i - 1);


		double TOT0 = 0.;
		double TOT1 = 0.;
		double TOT2 = 0.;

		// Liquid phase
		Eigen::VectorXd SC(7);
		SC.setZero();
		for (int j = NPA; j >= LC; j--)
		{

			Aggregation(y, j,
				SOM1, SOM2, SOM3, SOM4, SOM5, SOM6, SOM7, SY, SY1, SY12, SY2,
				NPA, LC);

			const int J1 = j + NPA;
			const int J2 = j + NPA * 2;

			if (j >= LC)
			{
				if (j >= NPA - 1)
				{
					M(IV-1) = 1;
					N(IV-1) = 1;
					LW(IV-1) = 1;
				}
				else 
				{
					IV++;

					if (IV > 19)
					{
						// do nothing
					}
					else
					{
						N(IV - 1) = 1;
					}

					if (IV > 16)
					{
						// do nothing
					}
					else
					{
						LW(IV - 1) = 1;
					}

					if (IV >= 5)
					{
						if (IV > 19)
						{
							// leave the if region
						}
						else
						{
							M(IV - 2 - 1) = 1;
						}
					}
					else
					{
						M(IV-1) = 1;
					}
				}

				for (int i = 1; i <= 7; i++)
					SC(i - 1) = 0.;

				for (int ISC = 0; ISC < 11; ISC++)
				{
					SC(0) += ABB(ISC)*SY(ISC);
					SC(1) += (ABB(ISC)*(M(ISC + 5)*CA*y(J1 + ISC + 8) + SY(ISC + 11)));
					SC(2) += ABB(ISC)*SY1(ISC);
					SC(3) += (ABB(ISC)*(N(ISC + 4)*CA*y(J1 + ISC + 5) + SY1(ISC + 11) + SY12(ISC)));
					SC(4) += (ABB(ISC)*(N(ISC + 7)*CA*y(J2 + ISC + 8) + SY1(ISC + 22)));
					SC(5) += (ABB(ISC)*(LW(ISC + 4)*CA*y(J2 + ISC + 5) + SY2(ISC + 11)));
					SC(6) += ABB(ISC)*SY2(ISC);
				}

				DY(j - 1) = AKP * (-(j - 4)*y(j-1) + SOM1 * (1. - ATOT)
					+ (1. - FR)*(SOM2 + M(3-1)*CA*y(J1 + 4 - 1))*(1. - ATOT)
					+ 0.5*SC(1-1) + 0.5*FR*SC(2-1));


				DY(J1 - 1) = AKP * (-((j - 6) + FR * CD + (1. - FR)*CA)*y(J1-1)
					+ SOM3 + FR * (N(2-1)*CD*y(J1 + 3 - 1) + SOM4)
					+ (N(3-1)*CA*y(J2 + 4 - 1) + SOM5)*(1. - ATOT)
					+ FR * (N(2-1)*y(J1 + 3 - 1) + SOM4)*(1. - ATOT)
					+ 0.5*SC(3-1) + 0.5*FR*SC(4-1) + 0.5*SC(5-1));


				DY(J2 - 1) = AKP * (-((j - 8) + CA + CD)*y(J2 - 1) + (1. - FR)*SOM6
					+ SOM7 + LW(2-1)*CD*y(J2 + 3 - 1)
					+ 0.5*SC(6-1) + 0.5*FR*SC(7-1));

				const double aaa = 0;
			}

		}	// end cycle for (int j = NPA; j >= LC; j--) (Liquid phase)

		// Gas phase
		{
			const double ABB1 = ABB(0);
			const double ABB2 = ABB(1);
			const double ABB3 = ABB(2);
			const double ABB4 = ABB(3);
			const double ABB5 = ABB(4);
			const double ABB6 = ABB(5);
			const double ABB7 = ABB(6);
			const double ABB8 = ABB(7);
			const double ABB9 = ABB(8);
			const double ABB10 = ABB(9);
			const double ABB11 = ABB(10);

			for (int J = LC - 1; J >= 2; J--)
			{
				Aggregation(y, J,
					SOM1, SOM2, SOM3, SOM4, SOM5, SOM6, SOM7, SY, SY1, SY12, SY2,
					NPA, LC);

				const double SY_1 = SY(0);
				const double SY_2 = SY(1);
				const double SY_3 = SY(2);
				const double SY_4 = SY(3);
				const double SY_5 = SY(4);
				const double SY_6 = SY(5);
				const double SY_7 = SY(6);
				const double SY_8 = SY(7);
				const double SY_9 = SY(8);
				const double SY_10 = SY(9);
				const double SY_11 = SY(10);
				const double SY_12 = SY(11);
				const double SY_13 = SY(12);
				const double SY_14 = SY(13);
				const double SY_15 = SY(14);
				const double SY_16 = SY(15);
				const double SY_17 = SY(16);
				const double SY_18 = SY(17);
				const double SY_19 = SY(18);
				const double SY_20 = SY(19);
				const double SY_21 = SY(20);
				const double SY_22 = SY(21);

				const double SY1_1 = SY1(0);
				const double SY1_2 = SY1(1);
				const double SY1_3 = SY1(2);
				const double SY1_4 = SY1(3);
				const double SY1_5 = SY1(4);
				const double SY1_6 = SY1(5);
				const double SY1_7 = SY1(6);
				const double SY1_8 = SY1(7);
				const double SY1_9 = SY1(8);
				const double SY1_10 = SY1(9);
				const double SY1_11 = SY1(10);
				const double SY1_12 = SY1(11);
				const double SY1_13 = SY1(12);
				const double SY1_14 = SY1(13);
				const double SY1_15 = SY1(14);
				const double SY1_16 = SY1(15);
				const double SY1_17 = SY1(16);
				const double SY1_18 = SY1(17);
				const double SY1_19 = SY1(18);
				const double SY1_20 = SY1(19);
				const double SY1_21 = SY1(20);
				const double SY1_22 = SY1(21);
				const double SY1_23 = SY1(22);
				const double SY1_24 = SY1(23);
				const double SY1_25 = SY1(24);
				const double SY1_26 = SY1(25);
				const double SY1_27 = SY1(26);
				const double SY1_28 = SY1(27);
				const double SY1_29 = SY1(28);
				const double SY1_30 = SY1(29);
				const double SY1_31 = SY1(30);
				const double SY1_32 = SY1(31);
				const double SY1_33 = SY1(32);

				const double SY12_1 = SY12(0);
				const double SY12_2 = SY12(1);
				const double SY12_3 = SY12(2);
				const double SY12_4 = SY12(3);
				const double SY12_5 = SY12(4);
				const double SY12_6 = SY12(5);
				const double SY12_7 = SY12(6);
				const double SY12_8 = SY12(7);
				const double SY12_9 = SY12(8);
				const double SY12_10 = SY12(9);
				const double SY12_11 = SY12(10);

				const double SY2_1 = SY2(0);
				const double SY2_2 = SY2(1);
				const double SY2_3 = SY2(2);
				const double SY2_4 = SY2(3);
				const double SY2_5 = SY2(4);
				const double SY2_6 = SY2(5);
				const double SY2_7 = SY2(6);
				const double SY2_8 = SY2(7);
				const double SY2_9 = SY2(8);
				const double SY2_10 = SY2(9);
				const double SY2_11 = SY2(10);
				const double SY2_12 = SY2(11);
				const double SY2_13 = SY2(12);
				const double SY2_14 = SY2(13);
				const double SY2_15 = SY2(14);
				const double SY2_16 = SY2(15);
				const double SY2_17 = SY2(16);
				const double SY2_18 = SY2(17);
				const double SY2_19 = SY2(18);
				const double SY2_20 = SY2(19);
				const double SY2_21 = SY2(20);
				const double SY2_22 = SY2(21);




				const int J1 = J + NPA - 1;
				const int J2 = J + NPA * 2 - 1;

				if (J >= 16)
				{
					DY(J - 1) = AKP * (
						SOM1*(1. - (ATOT))
						+ (1. - FR)*(SOM2)*(1. - (ATOT))
						+ 0.5*(ABB1*SY_1 + ABB2 * SY_2 + ABB3 * SY_3
							+ ABB4 * SY_4 + ABB5 * SY_5 + ABB6 * SY_6
							+ ABB7 * SY_7 + ABB8 * SY_8 + ABB9 * SY_9
							+ ABB10 * SY_10 + ABB11 * SY_11)
						+ 0.5*FR*(ABB1*(SY_12)
							+ABB2 * (SY_13)
							+ABB3 * (SY_14)
							+ABB4 * (SY_15)
							+ABB5 * (SY_16)
							+ABB6 * (SY_17)
							+ABB7 * (SY_18)
							+ABB8 * (SY_19)
							+ABB9 * (SY_20)
							+ABB10 * (SY_21)
							+ABB11 * (SY_22)));

					if (J >= LC) DY(J - 1) += -AKP * (J - 4)*y(J - 1);
					if (J + 4 >= LC) DY(J - 1) += AKP * (1. - FR)*CA*y(J1 + 4) * (1. - ATOT);
					if (J + 9 >= LC) DY(J - 1) += AKD * ABB1*y(J1 + 9);
					if (J + 10 >= LC) DY(J - 1) += AKD * ABB2*y(J1 + 10);
					if (J + 11 >= LC) DY(J - 1) += AKD * ABB3*y(J1 + 11);
					if (J + 12 >= LC) DY(J - 1) += AKD * ABB4*y(J1 + 12);
					if (J + 13 >= LC) DY(J - 1) += AKD * ABB5*y(J1 + 13);
					if (J + 14 >= LC) DY(J - 1) += AKD * ABB6*y(J1 + 14);
					if (J + 15 >= LC) DY(J - 1) += AKD * ABB7*y(J1 + 15);
					if (J + 16 >= LC) DY(J - 1) += AKD * ABB8*y(J1 + 16);
					if (J + 17 >= LC) DY(J - 1) += AKD * ABB9*y(J1 + 17);
					if (J + 18 >= LC) DY(J - 1) += AKD * ABB10*y(J1 + 18);
					if (J + 19 >= LC) DY(J - 1) += AKD * ABB11*y(J1 + 19);


					DY(J1) = AKP * (SOM3 +
						FR * (SOM4)
						+(SOM5)*(1. - (ATOT))
						+ FR * (SOM4)*(1. - (ATOT))
						+ 0.5*(ABB1*SY1_1 + ABB2 * SY1_2 + ABB3 * SY1_3
							+ ABB4 * SY1_4
							+ ABB5 * SY1_5 + ABB6 * SY1_6 + ABB7 * SY1_7
							+ ABB8 * SY1_8 + ABB9 * SY1_9 + ABB10 * SY1_10
							+ ABB11 * SY1_11)
						+ 0.5*FR*(ABB1*(SY1_12 + SY12_1)
							+ ABB2 * (SY1_13 + SY12_2)
							+ ABB3 * (SY1_14 + SY12_3)
							+ ABB4 * (SY1_15 + SY12_4)
							+ ABB5 * (SY1_16 + SY12_5)
							+ ABB6 * (SY1_17 + SY12_6)
							+ ABB7 * (SY1_18 + SY12_7)
							+ ABB8 * (SY1_19 + SY12_8)
							+ ABB9 * (SY1_20 + SY12_9)
							+ ABB10 * (SY1_21 + SY12_10)
							+ ABB11 * (SY1_22 + SY12_11))
						+ 0.5*(ABB1*(SY1_23)
							+ABB2 * (SY1_24)
							+ABB3 * (SY1_25)
							+ABB4 * (SY1_26)
							+ABB5 * (SY1_27)
							+ABB6 * (SY1_28)
							+ABB7 * (SY1_29)
							+ABB8 * (SY1_30)
							+ABB9 * (SY1_31)
							+ABB10 * (SY1_32)
							+ABB11 * (SY1_33)));

					if (J >= LC) DY(J1) += -AKP * ((J - 6) + FR * CD + (1. - FR) * CA)*y(J1);
					if (J + 3 >= LC) DY(J1) += AKP * FR*(CD*y(J1 + 3) + y(J1 + 3)*(1. - ATOT));
					if (J + 4 >= LC) DY(J1) += AKP * CA*y(J2 + 4)*(1. - ATOT);
					if (J + 6 >= LC) DY(J1) += AKD * ABB1*y(J1 + 6);
					if (J + 7 >= LC) DY(J1) += AKD * ABB2*y(J1 + 7);
					if (J + 8 >= LC) DY(J1) += AKD * ABB3*y(J1 + 8);
					if (J + 9 >= LC) DY(J1) += AKD * ABB4*y(J1 + 9);
					if (J + 10 >= LC) DY(J1) += AKD * ABB5*y(J1 + 10);
					if (J + 11 >= LC) DY(J1) += AKD * ABB6*y(J1 + 11);
					if (J + 12 >= LC) DY(J1) += AKD * ABB7*y(J1 + 12);
					if (J + 13 >= LC) DY(J1) += AKD * ABB8*y(J1 + 13);
					if (J + 14 >= LC) DY(J1) += AKD * ABB9*y(J1 + 14);
					if (J + 15 >= LC) DY(J1) += AKD * ABB10*y(J1 + 15);
					if (J + 16 >= LC) DY(J1) += AKD * ABB11*y(J1 + 16);
					if (J + 9 >= LC)  DY(J1) += AKA * ABB1*y(J2 + 9);
					if (J + 10 >= LC) DY(J1) += AKA * ABB2*y(J2 + 10);
					if (J + 11 >= LC) DY(J1) += AKA * ABB3*y(J2 + 11);
					if (J + 12 >= LC) DY(J1) += AKA * ABB4*y(J2 + 12);
					if (J + 13 >= LC) DY(J1) += AKA * ABB5*y(J2 + 13);
					if (J + 14 >= LC) DY(J1) += AKA * ABB6*y(J2 + 14);
					if (J + 15 >= LC) DY(J1) += AKA * ABB7*y(J2 + 15);
					if (J + 16 >= LC) DY(J1) += AKA * ABB8*y(J2 + 16);
					if (J + 17 >= LC) DY(J1) += AKA * ABB9*y(J2 + 17);
					if (J + 18 >= LC) DY(J1) += AKA * ABB10*y(J2 + 18);
					if (J + 19 >= LC) DY(J1) += AKA * ABB11*y(J2 + 19);

					DY(J2) = AKP * ((1. - FR)*SOM6 + SOM7
						+ 0.5*(ABB1*(SY2_12)
							+ABB2 * (SY2_13)
							+ABB3 * (SY2_14)
							+ABB4 * (SY2_15)
							+ABB5 * (SY2_16)
							+ABB6 * (SY2_17)
							+ABB7 * (SY2_18)
							+ABB8 * (SY2_19)
							+ABB9 * (SY2_20)
							+ABB10 * (SY2_21)
							+ABB11 * (SY2_22))
						+ 0.5*FR*(ABB1*SY2_1 + ABB2 * SY2_2 + ABB3 * SY2_3
							+ ABB4 * SY2_4 + ABB5 * SY2_5 + ABB6 * SY2_6
							+ ABB7 * SY2_7 + ABB8 * SY2_8 + ABB9 * SY2_9
							+ ABB10 * SY2_10 + ABB11 * SY2_11));

					if (J + 3 >= LC) DY(J2) += AKP * CD*y(J2 + 3);
					if (J + 6 >= LC) DY(J2) += AKA * ABB1*y(J2 + 6);
					if (J + 7 >= LC) DY(J2) += AKA * ABB2*y(J2 + 7);
					if (J + 8 >= LC) DY(J2) += AKA * ABB3*y(J2 + 8);
					if (J + 9 >= LC) DY(J2) += AKA * ABB4*y(J2 + 9);
					if (J + 10 >= LC) DY(J2) += AKA * ABB5*y(J2 + 10);
					if (J + 11 >= LC) DY(J2) += AKA * ABB6*y(J2 + 11);
					if (J + 12 >= LC) DY(J2) += AKA * ABB7*y(J2 + 12);
					if (J + 13 >= LC) DY(J2) += AKA * ABB8*y(J2 + 13);
					if (J + 14 >= LC) DY(J2) += AKA * ABB9*y(J2 + 14);
					if (J + 15 >= LC) DY(J2) += AKA * ABB10*y(J2 + 15);
					if (J + 16 >= LC) DY(J2) += AKA * ABB11*y(J2 + 16);

					const int aaa = 1;

				} // end if (J >= 16)

				else if (J == 15)
				{
					DY(J - 1) = AKP * (
						SOM1*(1. - (ATOT))
						+ (1. - FR)*(SOM2)*(1. - (ATOT))
						+ 0.5*(ABB1*SY_1 + ABB2 * SY_2 + ABB3 * SY_3
							+ ABB4 * SY_4 + ABB5 * SY_5 + ABB6 * SY_6
							+ ABB7 * SY_7 + ABB8 * SY_8 + ABB9 * SY_9
							+ ABB10 * SY_10 + ABB11 * SY_11)
						+ 0.5*FR*(ABB1*(SY_12)
							+ABB2 * (SY_13)
							+ABB3 * (SY_14)
							+ABB4 * (SY_15)
							+ABB5 * (SY_16)
							+ABB6 * (SY_17)
							+ABB7 * (SY_18)
							+ABB8 * (SY_19)
							+ABB9 * (SY_20)
							+ABB10 * (SY_21)
							+ABB11 * (SY_22)));

					if (J >= LC) DY(J - 1) += -AKP * (J - 4)*y(J - 1);
					if (J + 4 >= LC) DY(J - 1) += AKP * (1. - FR)*CA*y(J1 + 4) * (1 - ATOT);
					if (J + 9 >= LC) DY(J - 1) += AKD * ABB1*y(J1 + 9);
					if (J + 10 >= LC) DY(J - 1) += AKD * ABB2*y(J1 + 10);
					if (J + 11 >= LC) DY(J - 1) += AKD * ABB3*y(J1 + 11);
					if (J + 12 >= LC) DY(J - 1) += AKD * ABB4*y(J1 + 12);
					if (J + 13 >= LC) DY(J - 1) += AKD * ABB5*y(J1 + 13);
					if (J + 14 >= LC) DY(J - 1) += AKD * ABB6*y(J1 + 14);
					if (J + 15 >= LC) DY(J - 1) += AKD * ABB7*y(J1 + 15);
					if (J + 16 >= LC) DY(J - 1) += AKD * ABB8*y(J1 + 16);
					if (J + 17 >= LC) DY(J - 1) += AKD * ABB9*y(J1 + 17);
					if (J + 18 >= LC) DY(J - 1) += AKD * ABB10*y(J1 + 18);
					if (J + 19 >= LC) DY(J - 1) += AKD * ABB11*y(J1 + 19);


					DY(J1) = AKP * (SOM3 +
						FR * (SOM4)
						+(SOM5)*(1. - (ABB1 + ABB2 + ABB3 + ABB4
							+ ABB5 + ABB6 + ABB7 + ABB8 + ABB9 + ABB10))
						+ FR * (SOM4)*(1. - (ABB1 + ABB2 + ABB3 + ABB4
							+ ABB5 + ABB6 + ABB7 + ABB8 + ABB9 + ABB10))
						+ 0.5*(ABB1*SY1_1 + ABB2 * SY1_2 + ABB3 * SY1_3
							+ ABB4 * SY1_4
							+ ABB5 * SY1_5 + ABB6 * SY1_6 + ABB7 * SY1_7
							+ ABB8 * SY1_8 + ABB9 * SY1_9 + ABB10 * SY1_10)
						+ 0.5*FR*(ABB1*(SY1_12 + SY12_1)
							+ ABB2 * (SY1_13 + SY12_2)
							+ ABB3 * (SY1_14 + SY12_3)
							+ ABB4 * (SY1_15 + SY12_4)
							+ ABB5 * (SY1_16 + SY12_5)
							+ ABB6 * (SY1_17 + SY12_6)
							+ ABB7 * (SY1_18 + SY12_7)
							+ ABB8 * (SY1_19 + SY12_8)
							+ ABB9 * (SY1_20 + SY12_9)
							+ ABB10 * (SY1_21 + SY12_10))
						+ 0.5*(ABB1*(SY1_23)
							+ABB2 * (SY1_24)
							+ABB3 * (SY1_25)
							+ABB4 * (SY1_26)
							+ABB5 * (SY1_27)
							+ABB6 * (SY1_28)
							+ABB7 * (SY1_29)
							+ABB8 * (SY1_30)
							+ABB9 * (SY1_31)
							+ABB10 * (SY1_32))
						+ ABB11 * CORRO(J - 1));

					if (J >= LC) DY(J1) = DY(J1) - AKP * ((J - 6) + FR * CD + (1. - FR) * CA)*y(J1);
					if (J + 3 >= LC) DY(J1) += AKP * FR*(CD*y(J1 + 3) + y(J1 + 3)*(1. - (ATOT - ABB11)));
					if (J + 4 >= LC) DY(J1) += AKP * CA*y(J2 + 4) * (1. - (ATOT - ABB11));
					if (J + 6 >= LC) DY(J1) += AKD * ABB1*y(J1 + 6);
					if (J + 7 >= LC) DY(J1) += AKD * ABB2*y(J1 + 7);
					if (J + 8 >= LC) DY(J1) += AKD * ABB3*y(J1 + 8);
					if (J + 9 >= LC) DY(J1) += AKD * ABB4*y(J1 + 9);
					if (J + 10 >= LC) DY(J1) += AKD * ABB5*y(J1 + 10);
					if (J + 11 >= LC) DY(J1) += AKD * ABB6*y(J1 + 11);
					if (J + 12 >= LC) DY(J1) += AKD * ABB7*y(J1 + 12);
					if (J + 13 >= LC) DY(J1) += AKD * ABB8*y(J1 + 13);
					if (J + 14 >= LC) DY(J1) += AKD * ABB9*y(J1 + 14);
					if (J + 15 >= LC) DY(J1) += AKD * ABB10*y(J1 + 15);
					if (J + 9 >= LC) DY(J1) += 0.5*AKP*ABB1*CA*y(J2 + 9);
					if (J + 10 >= LC) DY(J1) += 0.5*AKP*ABB2*CA*y(J2 + 10);
					if (J + 11 >= LC) DY(J1) += 0.5*AKP*ABB3*CA*y(J2 + 11);
					if (J + 12 >= LC) DY(J1) += 0.5*AKP*ABB4*CA*y(J2 + 12);
					if (J + 13 >= LC) DY(J1) += 0.5*AKP*ABB5*CA*y(J2 + 13);
					if (J + 14 >= LC) DY(J1) += 0.5*AKP*ABB6*CA*y(J2 + 14);
					if (J + 15 >= LC) DY(J1) += 0.5*AKP*ABB7*CA*y(J2 + 15);
					if (J + 16 >= LC) DY(J1) += 0.5*AKP*ABB8*CA*y(J2 + 16);
					if (J + 17 >= LC) DY(J1) += 0.5*AKP*ABB9*CA*y(J2 + 17);
					if (J + 18 >= LC) DY(J1) += 0.5*AKP*ABB10*CA*y(J2 + 18);

					DY(J2) = AKP * ((1. - FR)*SOM6 + SOM7
						+ 0.5*(ABB1*(SY2_12)
							+ABB2 * (SY2_13)
							+ABB3 * (SY2_14)
							+ABB4 * (SY2_15)
							+ABB5 * (SY2_16)
							+ABB6 * (SY2_17)
							+ABB7 * (SY2_18)
							+ABB8 * (SY2_19)
							+ABB9 * (SY2_20)
							+ABB10 * (SY2_21)
							+ABB11 * (SY2_22))
						+ 0.5*FR*(ABB1*SY2_1 + ABB2 * SY2_2 + ABB3 * SY2_3
							+ ABB4 * SY2_4 + ABB5 * SY2_5 + ABB6 * SY2_6
							+ ABB7 * SY2_7 + ABB8 * SY2_8 + ABB9 * SY2_9
							+ ABB10 * SY2_10 + ABB11 * SY2_11));

					if (J >= LC) DY(J2) += -AKP * ((J - 8) + CA + CD)*y(J2);
					if (J + 3 >= LC) DY(J2) += AKP * CD*y(J2 + 3);
					if (J + 6 >= LC) DY(J2) += 0.5*AKP*ABB1*CA*y(J2 + 6);
					if (J + 7 >= LC) DY(J2) += 0.5*AKP*ABB2*CA*y(J2 + 7);
					if (J + 8 >= LC) DY(J2) += 0.5*AKP*ABB3*CA*y(J2 + 8);
					if (J + 9 >= LC) DY(J2) += 0.5*AKP*ABB4*CA*y(J2 + 9);
					if (J + 10 >= LC) DY(J2) += 0.5*AKP*ABB5*CA*y(J2 + 10);
					if (J + 11 >= LC) DY(J2) += 0.5*AKP*ABB6*CA*y(J2 + 11);
					if (J + 12 >= LC) DY(J2) += 0.5*AKP*ABB7*CA*y(J2 + 12);
					if (J + 13 >= LC) DY(J2) += 0.5*AKP*ABB8*CA*y(J2 + 13);
					if (J + 14 >= LC) DY(J2) += 0.5*AKP*ABB9*CA*y(J2 + 14);
					if (J + 15 >= LC) DY(J2) += 0.5*AKP*ABB10*CA*y(J2 + 15);
					if (J + 16 >= LC) DY(J2) += 0.5*AKP*ABB11*CA*y(J2 + 16);

					const int aaa = 1;	// TODO
				}

				else if (J == 14)
				{
					DY(J - 1) = AKP * (
						SOM1*(1. - (ABB1 + ABB2 + ABB3 + ABB4 + ABB5 + ABB6
							+ ABB7 + ABB8 + ABB9 + ABB10))
						+ (1. - FR)*(SOM2)*(1. - (ABB1 + ABB2 + ABB3 + ABB4
							+ ABB5 + ABB6 + ABB7 + ABB8 + ABB9 + ABB10))
						+ 0.5*(ABB1*SY_1 + ABB2 * SY_2 + ABB3 * SY_3
							+ ABB4 * SY_4 + ABB5 * SY_5 + ABB6 * SY_6
							+ ABB7 * SY_7 + ABB8 * SY_8 + ABB9 * SY_9
							+ ABB10 * SY_10 + ABB11 * SY_11)
						+ 0.5*FR*(ABB1*(SY_12)
							+ABB2 * (SY_13)
							+ABB3 * (SY_14)
							+ABB4 * (SY_15)
							+ABB5 * (SY_16)
							+ABB6 * (SY_17)
							+ABB7 * (SY_18)
							+ABB8 * (SY_19)
							+ABB9 * (SY_20)
							+ABB10 * (SY_21)
							+ABB11 * (SY_22)));


					if (J >= LC) DY(J - 1) += -AKP * (J - 4)*y(J - 1);
					if (J + 4 >= LC) DY(J - 1) += AKP * (1. - FR)*CA*y(J1 + 4)* (1 - (ATOT - ABB11));
					if (J + 9 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB1*CA*y(J1 + 9);
					if (J + 10 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB2*CA*y(J1 + 10);
					if (J + 11 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB3*CA*y(J1 + 11);
					if (J + 12 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB4*CA*y(J1 + 12);
					if (J + 13 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB5*CA*y(J1 + 13);
					if (J + 14 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB6*CA*y(J1 + 14);
					if (J + 15 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB7*CA*y(J1 + 15);
					if (J + 16 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB8*CA*y(J1 + 16);
					if (J + 17 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB9*CA*y(J1 + 17);
					if (J + 18 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB10*CA*y(J1 + 18);
					if (J + 19 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB11*CA*y(J1 + 19);


					DY(J1) = AKP * (SOM3 +
						FR * (SOM4)
						+(SOM5)*(1. - (ABB1 + ABB2 + ABB3 + ABB4
							+ ABB5 + ABB6 + ABB7 + ABB8 + ABB9))
						+ FR * (SOM4)*(1. - (ABB1 + ABB2 + ABB3 + ABB4
							+ ABB5 + ABB6 + ABB7 + ABB8 + ABB9))
						+ 0.5*(ABB1*SY1_1 + ABB2 * SY1_2 + ABB3 * SY1_3
							+ ABB4 * SY1_4
							+ ABB5 * SY1_5 + ABB6 * SY1_6 + ABB7 * SY1_7
							+ ABB8 * SY1_8 + ABB9 * SY1_9 + ABB11 * SY1_11)
						+ 0.5*FR*(ABB1*(SY1_12 + SY12_1)
							+ ABB2 * (SY1_13 + SY12_2)
							+ ABB3 * (SY1_14 + SY12_3)
							+ ABB4 * (SY1_15 + SY12_4)
							+ ABB5 * (SY1_16 + SY12_5)
							+ ABB6 * (SY1_17 + SY12_6)
							+ ABB7 * (SY1_18 + SY12_7)
							+ ABB8 * (SY1_19 + SY12_8)
							+ ABB9 * (SY1_20 + SY12_9)
							+ ABB11 * (SY1_22 + SY12_11))
						+ 0.5*(ABB1*(SY1_23)
							+ABB2 * (SY1_24)
							+ABB3 * (SY1_25)
							+ABB4 * (SY1_26)
							+ABB5 * (SY1_27)
							+ABB6 * (SY1_28)
							+ABB7 * (SY1_29)
							+ABB8 * (SY1_30)
							+ABB9 * (SY1_31)
							+ABB11 * (SY1_33))
						+ ABB10 * CORRO(J - 1));

					if (J >= LC) DY(J1) = DY(J1) - AKP * ((J - 6) + FR * CD + (1. - FR) *CA)*y(J1);
					if (J + 3 >= LC) DY(J1) += AKP * FR*(CD*y(J1 + 3) + y(J1 + 3)*(1. - (ATOT - ABB11 - ABB10)));
					if (J + 4 >= LC) DY(J1) += AKP * CA*y(J2 + 4) * (1. - (ATOT - ABB11 - ABB10));
					if (J + 6 >= LC) DY(J1) += 0.5*AKP*FR*ABB1*CA*y(J1 + 6);
					if (J + 7 >= LC) DY(J1) += 0.5*AKP*FR*ABB2*CA*y(J1 + 7);
					if (J + 8 >= LC) DY(J1) += 0.5*AKP*FR*ABB3*CA*y(J1 + 8);
					if (J + 9 >= LC) DY(J1) += 0.5*AKP*FR*ABB4*CA*y(J1 + 9);
					if (J + 10 >= LC) DY(J1) += 0.5*AKP*FR*ABB5*CA*y(J1 + 10);
					if (J + 11 >= LC) DY(J1) += 0.5*AKP*FR*ABB6*CA * y(J1 + 11);
					if (J + 12 >= LC) DY(J1) += 0.5*AKP*FR*ABB7*CA * y(J1 + 12);
					if (J + 13 >= LC) DY(J1) += 0.5*AKP*FR*ABB8*CA * y(J1 + 13);
					if (J + 14 >= LC) DY(J1) += 0.5*AKP*FR*ABB9*CA *y(J1 + 14);
					if (J + 16 >= LC) DY(J1) += 0.5*AKP*FR*ABB11*CA *y(J1 + 16);
					if (J + 9 >= LC) DY(J1) += 0.5*AKP*ABB1*CA*y(J2 + 9);
					if (J + 10 >= LC) DY(J1) += 0.5*AKP*ABB2*CA*y(J2 + 10);
					if (J + 11 >= LC) DY(J1) += 0.5*AKP*ABB3*CA*y(J2 + 11);
					if (J + 12 >= LC) DY(J1) += 0.5*AKP*ABB4*CA*y(J2 + 12);
					if (J + 13 >= LC) DY(J1) += 0.5*AKP*ABB5*CA*y(J2 + 13);
					if (J + 14 >= LC) DY(J1) += 0.5*AKP*ABB6*CA*y(J2 + 14);
					if (J + 15 >= LC) DY(J1) += 0.5*AKP*ABB7*CA*y(J2 + 15);
					if (J + 16 >= LC) DY(J1) += 0.5*AKP*ABB8*CA*y(J2 + 16);
					if (J + 17 >= LC) DY(J1) += 0.5*AKP*ABB9*CA*y(J2 + 17);
					if (J + 19 >= LC) DY(J1) += 0.5*AKP*ABB11*CA*y(J2 + 19);


					DY(J2) = AKP * ((1. - FR)*SOM6 + SOM7
						+ 0.5*(ABB1*(SY2_12)
							+ABB2 * (SY2_13)
							+ABB3 * (SY2_14)
							+ABB4 * (SY2_15)
							+ABB5 * (SY2_16)
							+ABB6 * (SY2_17)
							+ABB7 * (SY2_18)
							+ABB8 * (SY2_19)
							+ABB9 * (SY2_20)
							+ABB10 * (SY2_21)
							+ABB11 * (SY2_22))
						+ 0.5*FR*(ABB1*SY2_1 + ABB2 * SY2_2 + ABB3 * SY2_3
							+ ABB4 * SY2_4 + ABB5 * SY2_5 + ABB6 * SY2_6
							+ ABB7 * SY2_7 + ABB8 * SY2_8 + ABB9 * SY2_9
							+ ABB10 * SY2_10 + ABB11 * SY2_11));

					if (J >= LC) DY(J2) = DY(J2) - AKP * ((J - 8) + CA + CD)*y(J2);
					if (J + 3 >= LC) DY(J2) += AKP * CD*y(J2 + 3);
					if (J + 6 >= LC) DY(J2) += 0.5*AKP*ABB1*CA*y(J2 + 6);
					if (J + 7 >= LC) DY(J2) += 0.5*AKP*ABB2*CA*y(J2 + 7);
					if (J + 8 >= LC) DY(J2) += 0.5*AKP*ABB3*CA*y(J2 + 8);
					if (J + 9 >= LC) DY(J2) += 0.5*AKP*ABB4*CA*y(J2 + 9);
					if (J + 10 >= LC) DY(J2) += 0.5*AKP*ABB5*CA*y(J2 + 10);
					if (J + 11 >= LC) DY(J2) += 0.5*AKP*ABB6*CA*y(J2 + 11);
					if (J + 12 >= LC) DY(J2) += 0.5*AKP*ABB7*CA*y(J2 + 12);
					if (J + 13 >= LC) DY(J2) += 0.5*AKP*ABB8*CA*y(J2 + 13);
					if (J + 14 >= LC) DY(J2) += 0.5*AKP*ABB9*CA*y(J2 + 14);
					if (J + 15 >= LC) DY(J2) += 0.5*AKP*ABB10*CA*y(J2 + 15);
					if (J + 16 >= LC) DY(J2) += 0.5*AKP*ABB11*CA*y(J2 + 16);

					const int aaa = 1;	// TODO
				}

				else if (J == 13)
				{
					DY(J - 1) = AKP * (
						SOM1*(1. - (ABB1 + ABB2 + ABB3 + ABB4 + ABB5 + ABB6
							+ ABB7 + ABB8 + ABB9))
						+ (1. - FR)*(SOM2)*(1. - (ABB1 + ABB2 + ABB3 + ABB4
							+ ABB5 + ABB6 + ABB7 + ABB8 + ABB9))
						+ 0.5*(ABB1*SY_1 + ABB2 * SY_2 + ABB3 * SY_3
							+ ABB4 * SY_4 + ABB5 * SY_5 + ABB6 * SY_6
							+ ABB7 * SY_7 + ABB8 * SY_8 + ABB9 * SY_9
							+ ABB10 * SY_10 + ABB11 * SY_11)
						+ 0.5*FR*(ABB1*(SY_12)
							+ABB2 * (SY_13)
							+ABB3 * (SY_14)
							+ABB4 * (SY_15)
							+ABB5 * (SY_16)
							+ABB6 * (SY_17)
							+ABB7 * (SY_18)
							+ABB8 * (SY_19)
							+ABB9 * (SY_20)
							+ABB10 * (SY_21)
							+ABB11 * (SY_22)));

					if (J >= LC) DY(J - 1) += -AKP * (J - 4)*y(J - 1);
					if (J + 4 >= LC) DY(J - 1) += AKP * (1. - FR)*CA*y(J1 + 4)* (1 - (ATOT - ABB11 - ABB10));
					if (J + 9 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB1*CA*y(J1 + 9);
					if (J + 10 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB2*CA*y(J1 + 10);
					if (J + 11 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB3*CA*y(J1 + 11);
					if (J + 12 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB4*CA*y(J1 + 12);
					if (J + 13 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB5*CA*y(J1 + 13);
					if (J + 14 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB6*CA*y(J1 + 14);
					if (J + 15 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB7*CA*y(J1 + 15);
					if (J + 16 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB8*CA*y(J1 + 16);
					if (J + 17 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB9*CA*y(J1 + 17);
					if (J + 18 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB10*CA*y(J1 + 18);
					if (J + 19 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB11*CA*y(J1 + 19);


					DY(J1) = AKP * (SOM3 +
						FR * (SOM4)
						+(SOM5)*(1. - (ABB1 + ABB2 + ABB3 + ABB4
							+ ABB5 + ABB6 + ABB7 + ABB8))
						+ FR * (SOM4)*(1. - (ABB1 + ABB2 + ABB3 + ABB4
							+ ABB5 + ABB6 + ABB7 + ABB8))
						+ 0.5*(ABB1*SY1_1 + ABB2 * SY1_2 + ABB3 * SY1_3
							+ ABB4 * SY1_4
							+ ABB5 * SY1_5 + ABB6 * SY1_6 + ABB7 * SY1_7
							+ ABB8 * SY1_8 + ABB10 * SY1_10 + ABB11 * SY1_11)
						+ 0.5*FR*(ABB1*(SY1_12 + SY12_1)
							+ ABB2 * (SY1_13 + SY12_2)
							+ ABB3 * (SY1_14 + SY12_3)
							+ ABB4 * (SY1_15 + SY12_4)
							+ ABB5 * (SY1_16 + SY12_5)
							+ ABB6 * (SY1_17 + SY12_6)
							+ ABB7 * (SY1_18 + SY12_7)
							+ ABB8 * (SY1_19 + SY12_8)
							+ ABB10 * (SY1_21 + SY12_10)
							+ ABB11 * (SY1_22 + SY12_11))
						+ 0.5*(ABB1*(SY1_23)
							+ABB2 * (SY1_24)
							+ABB3 * (SY1_25)
							+ABB4 * (SY1_26)
							+ABB5 * (SY1_27)
							+ABB6 * (SY1_28)
							+ABB7 * (SY1_29)
							+ABB8 * (SY1_30)
							+ABB10 * (SY1_32)
							+ABB11 * (SY1_33))
						+ ABB9 * CORRO(J - 1));

					if (J >= LC) DY(J1) = DY(J1) - AKP * ((J - 6) + FR * CD + (1. - FR) * CA)*y(J1);
					if (J + 3 >= LC) DY(J1) += AKP * FR*(CD*y(J1 + 3) + y(J1 + 3)*(1. - (ATOT - ABB11 - ABB10 - ABB9)));
					if (J + 4 >= LC) DY(J1) += AKP * CA*y(J2 + 4) *(1. - (ATOT - ABB11 - ABB10 - ABB9));
					if (J + 6 >= LC) DY(J1) += 0.5*AKP*FR*ABB1*CA*y(J1 + 6);
					if (J + 7 >= LC) DY(J1) += 0.5*AKP*FR*ABB2*CA*y(J1 + 7);
					if (J + 8 >= LC) DY(J1) += 0.5*AKP*FR*ABB3*CA*y(J1 + 8);
					if (J + 9 >= LC) DY(J1) += 0.5*AKP*FR*ABB4*CA*y(J1 + 9);
					if (J + 10 >= LC) DY(J1) += 0.5*AKP*FR*ABB5*CA*y(J1 + 10);
					if (J + 11 >= LC) DY(J1) += 0.5*AKP*FR*ABB6*CA * y(J1 + 11);
					if (J + 12 >= LC) DY(J1) += 0.5*AKP*FR*ABB7*CA * y(J1 + 12);
					if (J + 13 >= LC) DY(J1) += 0.5*AKP*FR*ABB8*CA *y(J1 + 13);
					if (J + 15 >= LC) DY(J1) += 0.5*AKP*FR*ABB10*CA *y(J1 + 15);
					if (J + 16 >= LC) DY(J1) += 0.5*AKP*FR*ABB11*CA * y(J1 + 16);
					if (J + 9 >= LC) DY(J1) += 0.5*AKP*ABB1*CA*y(J2 + 9);
					if (J + 10 >= LC) DY(J1) += 0.5*AKP*ABB2*CA*y(J2 + 10);
					if (J + 11 >= LC) DY(J1) += 0.5*AKP*ABB3*CA*y(J2 + 11);
					if (J + 12 >= LC) DY(J1) += 0.5*AKP*ABB4*CA*y(J2 + 12);
					if (J + 13 >= LC) DY(J1) += 0.5*AKP*ABB5*CA*y(J2 + 13);
					if (J + 14 >= LC) DY(J1) += 0.5*AKP*ABB6*CA*y(J2 + 14);
					if (J + 15 >= LC) DY(J1) += 0.5*AKP*ABB7*CA*y(J2 + 15);
					if (J + 16 >= LC) DY(J1) += 0.5*AKP*ABB8*CA*y(J2 + 16);
					if (J + 18 >= LC) DY(J1) += 0.5*AKP*ABB10*CA*y(J2 + 18);
					if (J + 19 >= LC) DY(J1) += 0.5*AKP*ABB11*CA*y(J2 + 19);


					DY(J2) = AKP * ((1. - FR)*SOM6 + SOM7
						+ 0.5*(ABB1*(SY2_12)
							+ABB2 * (SY2_13)
							+ABB3 * (SY2_14)
							+ABB4 * (SY2_15)
							+ABB5 * (SY2_16)
							+ABB6 * (SY2_17)
							+ABB7 * (SY2_18)
							+ABB8 * (SY2_19)
							+ABB9 * (SY2_20)
							+ABB10 * (SY2_21)
							+ABB11 * (SY2_22))
						+ 0.5*FR*(ABB1*SY2_1 + ABB2 * SY2_2 + ABB3 * SY2_3
							+ ABB4 * SY2_4 + ABB5 * SY2_5 + ABB6 * SY2_6
							+ ABB7 * SY2_7 + ABB8 * SY2_8 + ABB9 * SY2_9
							+ ABB10 * SY2_10 + ABB11 * SY2_11));

					if (J >= LC) DY(J2) = DY(J2) - AKP * ((J - 8) + CA + CD)*y(J2);
					if (J + 3 >= LC) DY(J2) += AKP * CD*y(J2 + 3);
					if (J + 6 >= LC) DY(J2) += 0.5*AKP*ABB1*CA*y(J2 + 6);
					if (J + 7 >= LC) DY(J2) += 0.5*AKP*ABB2*CA*y(J2 + 7);
					if (J + 8 >= LC) DY(J2) += 0.5*AKP*ABB3*CA*y(J2 + 8);
					if (J + 9 >= LC) DY(J2) += 0.5*AKP*ABB4*CA*y(J2 + 9);
					if (J + 10 >= LC) DY(J2) += 0.5*AKP*ABB5*CA*y(J2 + 10);
					if (J + 11 >= LC) DY(J2) += 0.5*AKP*ABB6*CA*y(J2 + 11);
					if (J + 12 >= LC) DY(J2) += 0.5*AKP*ABB7*CA*y(J2 + 12);
					if (J + 13 >= LC) DY(J2) += 0.5*AKP*ABB8*CA*y(J2 + 13);
					if (J + 14 >= LC) DY(J2) += 0.5*AKP*ABB9*CA*y(J2 + 14);
					if (J + 15 >= LC) DY(J2) += 0.5*AKP*ABB10*CA*y(J2 + 15);
					if (J + 16 >= LC) DY(J2) += 0.5*AKP*ABB11*CA*y(J2 + 16);

					const int aaa = 1;	// TODO
				}

				else if (J == 12)
				{
					DY(J - 1) = AKP * (
						SOM1*(1. - (ABB1 + ABB2 + ABB3 + ABB4 + ABB5 + ABB6
							+ ABB7 + ABB8))
						+ (1. - FR)*(SOM2)*(1. - (ABB1 + ABB2 + ABB3 + ABB4
							+ ABB5 + ABB6 + ABB7 + ABB8))
						+ 0.5*(ABB1*SY_1 + ABB2 * SY_2 + ABB3 * SY_3
							+ ABB4 * SY_4 + ABB5 * SY_5 + ABB6 * SY_6
							+ ABB7 * SY_7 + ABB8 * SY_8 + ABB9 * SY_9
							+ ABB10 * SY_10)
						+ 0.5*FR*(ABB1*(SY_12)
							+ABB2 * (SY_13)
							+ABB3 * (SY_14)
							+ABB4 * (SY_15)
							+ABB5 * (SY_16)
							+ABB6 * (SY_17)
							+ABB7 * (SY_18)
							+ABB8 * (SY_19)
							+ABB9 * (SY_20)
							+ABB10 * (SY_21))
						+ ABB11 * CORRP(3 + J - 1));


					if (J >= LC) DY(J - 1) += -AKP * (J - 4)*y(J - 1);
					if (J + 4 >= LC) DY(J - 1) += AKP * (1. - FR)*CA*y(J1 + 4)* (1 - (ATOT - ABB11 - ABB10 - ABB9));
					if (J + 9 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB1*CA*y(J1 + 9);
					if (J + 10 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB2*CA*y(J1 + 10);
					if (J + 11 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB3*CA*y(J1 + 11);
					if (J + 12 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB4*CA*y(J1 + 12);
					if (J + 13 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB5*CA*y(J1 + 13);
					if (J + 14 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB6*CA*y(J1 + 14);
					if (J + 15 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB7*CA*y(J1 + 15);
					if (J + 16 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB8*CA*y(J1 + 16);
					if (J + 17 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB9*CA*y(J1 + 17);
					if (J + 18 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB10*CA*y(J1 + 18);


					DY(J1) = AKP * (SOM3 +
						FR * (SOM4)
						+(SOM5)*(1. - (ABB1 + ABB2 + ABB3 + ABB4
							+ ABB5 + ABB6 + ABB7))
						+ FR * (SOM4)*(1. - (ABB1 + ABB2 + ABB3 + ABB4
							+ ABB5 + ABB6 + ABB7))
						+ 0.5*(ABB1*SY1_1 + ABB2 * SY1_2 + ABB3 * SY1_3
							+ ABB4 * SY1_4
							+ ABB5 * SY1_5 + ABB6 * SY1_6 + ABB7 * SY1_7
							+ ABB9 * SY1_9 + ABB10 * SY1_10 + ABB11 * SY1_11)
						+ 0.5*FR*(ABB1*(SY1_12 + SY12_1)
							+ ABB2 * (SY1_13 + SY12_2)
							+ ABB3 * (SY1_14 + SY12_3)
							+ ABB4 * (SY1_15 + SY12_4)
							+ ABB5 * (SY1_16 + SY12_5)
							+ ABB6 * (SY1_17 + SY12_6)
							+ ABB7 * (SY1_18 + SY12_7)
							+ ABB9 * (SY1_20 + SY12_9)
							+ ABB10 * (SY1_21 + SY12_10)
							+ ABB11 * (SY1_22 + SY12_11))
						+ 0.5*(ABB1*(SY1_23)
							+ABB2 * (SY1_24)
							+ABB3 * (SY1_25)
							+ABB4 * (SY1_26)
							+ABB5 * (SY1_27)
							+ABB6 * (SY1_28)
							+ABB7 * (SY1_29)
							+ABB9 * (SY1_31)
							+ABB10 * (SY1_32)
							+ABB11 * (SY1_33))
						+ ABB8 * CORRO(J - 1));

					if (J >= LC) DY(J1) = DY(J1) - AKP * ((J - 6) + FR * CD + (1. - FR) *CA)*y(J1);
					if (J + 3 >= LC) DY(J1) += AKP * FR*(CD*y(J1 + 3) + y(J1 + 3)*(1. - (ATOT - ABB11 - ABB10 - ABB9 - ABB8)));
					if (J + 4 >= LC) DY(J1) += AKP * CA*y(J2 + 4) * (1. - (ATOT - ABB11 - ABB10 - ABB9 - ABB8));
					if (J + 6 >= LC) DY(J1) += 0.5*AKP*FR*ABB1*CA*y(J1 + 6);
					if (J + 7 >= LC) DY(J1) += 0.5*AKP*FR*ABB2*CA*y(J1 + 7);
					if (J + 8 >= LC) DY(J1) += 0.5*AKP*FR*ABB3*CA*y(J1 + 8);
					if (J + 9 >= LC) DY(J1) += 0.5*AKP*FR*ABB4*CA*y(J1 + 9);
					if (J + 10 >= LC) DY(J1) += 0.5*AKP*FR*ABB5*CA*y(J1 + 10);
					if (J + 11 >= LC) DY(J1) += 0.5*AKP*FR*ABB6*CA * y(J1 + 11);
					if (J + 12 >= LC) DY(J1) += 0.5*AKP*FR*ABB7*CA *y(J1 + 12);
					if (J + 14 >= LC) DY(J1) += 0.5*AKP*FR*ABB9*CA * y(J1 + 14);
					if (J + 15 >= LC) DY(J1) += 0.5*AKP*FR*ABB10*CA *y(J1 + 15);
					if (J + 16 >= LC) DY(J1) += 0.5*AKP*FR*ABB11*CA *  y(J1 + 16);
					if (J + 9 >= LC) DY(J1) += 0.5*AKP*ABB1*CA*y(J2 + 9);
					if (J + 10 >= LC) DY(J1) += 0.5*AKP*ABB2*CA*y(J2 + 10);
					if (J + 11 >= LC) DY(J1) += 0.5*AKP*ABB3*CA*y(J2 + 11);
					if (J + 12 >= LC) DY(J1) += 0.5*AKP*ABB4*CA*y(J2 + 12);
					if (J + 13 >= LC) DY(J1) += 0.5*AKP*ABB5*CA*y(J2 + 13);
					if (J + 14 >= LC) DY(J1) += 0.5*AKP*ABB6*CA*y(J2 + 14);
					if (J + 15 >= LC) DY(J1) += 0.5*AKP*ABB7*CA*y(J2 + 15);
					if (J + 17 >= LC) DY(J1) += 0.5*AKP*ABB9*CA*y(J2 + 17);
					if (J + 18 >= LC) DY(J1) += 0.5*AKP*ABB10*CA*y(J2 + 18);
					if (J + 19 >= LC) DY(J1) += 0.5*AKP*ABB11*CA*y(J2 + 19);

					DY(J2) = AKP * ((1. - FR)*SOM6 + SOM7
						+ 0.5*(ABB1*(SY2_12)
							+ABB2 * (SY2_13)
							+ABB3 * (SY2_14)
							+ABB4 * (SY2_15)
							+ABB5 * (SY2_16)
							+ABB6 * (SY2_17)
							+ABB7 * (SY2_18)
							+ABB8 * (SY2_19)
							+ABB9 * (SY2_20)
							+ABB10 * (SY2_21)
							+ABB11 * (SY2_22))
						+ 0.5*FR*(ABB1*SY2_1 + ABB2 * SY2_2 + ABB3 * SY2_3
							+ ABB4 * SY2_4 + ABB5 * SY2_5 + ABB6 * SY2_6
							+ ABB7 * SY2_7 + ABB8 * SY2_8 + ABB9 * SY2_9
							+ ABB10 * SY2_10 + ABB11 * SY2_11));

					if (J >= LC) DY(J2) = DY(J2) - AKP * ((J - 8) + CA + CD)*y(J2);
					if (J + 3 >= LC) DY(J2) += AKP * CD*y(J2 + 3);
					if (J + 6 >= LC) DY(J2) += 0.5*AKP*ABB1*CA*y(J2 + 6);
					if (J + 7 >= LC) DY(J2) += 0.5*AKP*ABB2*CA*y(J2 + 7);
					if (J + 8 >= LC) DY(J2) += 0.5*AKP*ABB3*CA*y(J2 + 8);
					if (J + 9 >= LC) DY(J2) += 0.5*AKP*ABB4*CA*y(J2 + 9);
					if (J + 10 >= LC) DY(J2) += 0.5*AKP*ABB5*CA*y(J2 + 10);
					if (J + 11 >= LC) DY(J2) += 0.5*AKP*ABB6*CA*y(J2 + 11);
					if (J + 12 >= LC) DY(J2) += 0.5*AKP*ABB7*CA*y(J2 + 12);
					if (J + 13 >= LC) DY(J2) += 0.5*AKP*ABB8*CA*y(J2 + 13);
					if (J + 14 >= LC) DY(J2) += 0.5*AKP*ABB9*CA*y(J2 + 14);
					if (J + 15 >= LC) DY(J2) += 0.5*AKP*ABB10*CA*y(J2 + 15);
					if (J + 16 >= LC) DY(J2) += 0.5*AKP*ABB11*CA*y(J2 + 16);

					const int aaa = 1;	// TODO
				}

				else if (J == 11)
				{
					DY(J - 1) = AKP * (
						SOM1*(1. - (ABB1 + ABB2 + ABB3 + ABB4 + ABB5 + ABB6 + ABB7))
						+ (1. - FR)*(SOM2)*(1. - (ABB1 + ABB2 + ABB3 + ABB4
							+ ABB5 + ABB6 + ABB7))
						+ 0.5*(ABB1*SY_1 + ABB2 * SY_2 + ABB3 * SY_3
							+ ABB4 * SY_4 + ABB5 * SY_5 + ABB6 * SY_6
							+ ABB7 * SY_7 + ABB8 * SY_8 + ABB9 * SY_9
							+ ABB11 * SY_11)
						+ 0.5*FR*(ABB1*(SY_12)
							+ABB2 * (SY_13)
							+ABB3 * (SY_14)
							+ABB4 * (SY_15)
							+ABB5 * (SY_16)
							+ABB6 * (SY_17)
							+ABB7 * (SY_18)
							+ABB8 * (SY_19)
							+ABB9 * (SY_20)
							+ABB11 * (SY_22))
						+ ABB10 * CORRP(3 + J - 1));


					if (J >= LC) DY(J - 1) += -AKP * (J - 4)*y(J - 1);
					if (J + 4 >= LC) DY(J - 1) += AKP * (1. - FR)*CA*y(J1 + 4)* (1 - (ATOT - ABB11 - ABB10 - ABB9 - ABB8));
					if (J + 9 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB1*CA*y(J1 + 9);
					if (J + 10 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB2*CA*y(J1 + 10);
					if (J + 11 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB3*CA*y(J1 + 11);
					if (J + 12 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB4*CA*y(J1 + 12);
					if (J + 13 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB5*CA*y(J1 + 13);
					if (J + 14 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB6*CA*y(J1 + 14);
					if (J + 15 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB7*CA*y(J1 + 15);
					if (J + 16 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB8*CA*y(J1 + 16);
					if (J + 17 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB9*CA*y(J1 + 17);
					if (J + 19 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB11*CA*y(J1 + 19);


					DY(J1) = AKP * (SOM3 +
						FR * (SOM4)
						+(SOM5)*(1. - (ABB1 + ABB2 + ABB3 + ABB4 + ABB5
							+ ABB6))
						+ FR * (SOM4)*(1. - (ABB1 + ABB2 + ABB3 + ABB4 + ABB5
							+ ABB6))
						+ 0.5*(ABB1*SY1_1 + ABB2 * SY1_2 + ABB3 * SY1_3
							+ ABB4 * SY1_4
							+ ABB5 * SY1_5 + ABB6 * SY1_6 + ABB8 * SY1_8
							+ ABB9 * SY1_9 + ABB10 * SY1_10 + ABB11 * SY1_11)
						+ 0.5*FR*(ABB1*(SY1_12 + SY12_1)
							+ ABB2 * (SY1_13 + SY12_2)
							+ ABB3 * (SY1_14 + SY12_3)
							+ ABB4 * (SY1_15 + SY12_4)
							+ ABB5 * (SY1_16 + SY12_5)
							+ ABB6 * (SY1_17 + SY12_6)
							+ ABB8 * (SY1_19 + SY12_8)
							+ ABB9 * (SY1_20 + SY12_9)
							+ ABB10 * (SY1_21 + SY12_10)
							+ ABB11 * (SY1_22 + SY12_11))
						+ 0.5*(ABB1*(SY1_23)
							+ABB2 * (SY1_24)
							+ABB3 * (SY1_25)
							+ABB4 * (SY1_26)
							+ABB5 * (SY1_27)
							+ABB6 * (SY1_28)
							+ABB8 * (SY1_30)
							+ABB9 * (SY1_31)
							+ABB10 * (SY1_32)
							+ABB11 * (SY1_33))
						+ ABB7 * CORRO(J - 1));


					if (J >= LC) DY(J1) = DY(J1) - AKP * ((J - 6) + FR * CD + (1. - FR) *CA)*y(J1);
					if (J + 3 >= LC) DY(J1) += AKP * FR*(CD*y(J1 + 3) + y(J1 + 3)*(1. - (ATOT - ABB11 - ABB10 - ABB9 - ABB8 - ABB7)));
					if (J + 4 >= LC) DY(J1) += AKP * CA*y(J2 + 4) *(1. - (ATOT - ABB11 - ABB10 - ABB9 - ABB8 - ABB7));
					if (J + 6 >= LC) DY(J1) += 0.5*AKP*FR*ABB1*CA*y(J1 + 6);
					if (J + 7 >= LC) DY(J1) += 0.5*AKP*FR*ABB2*CA*y(J1 + 7);
					if (J + 8 >= LC) DY(J1) += 0.5*AKP*FR*ABB3*CA*y(J1 + 8);
					if (J + 9 >= LC) DY(J1) += 0.5*AKP*FR*ABB4*CA*y(J1 + 9);
					if (J + 10 >= LC) DY(J1) += 0.5*AKP*FR*ABB5*CA*y(J1 + 10);
					if (J + 11 >= LC) DY(J1) += 0.5*AKP*FR*ABB6*CA * y(J1 + 11);
					if (J + 13 >= LC) DY(J1) += 0.5*AKP*FR*ABB8*CA *y(J1 + 13);
					if (J + 14 >= LC) DY(J1) += 0.5*AKP*FR*ABB9*CA *y(J1 + 14);
					if (J + 15 >= LC) DY(J1) += 0.5*AKP*FR*ABB10*CA * y(J1 + 15);
					if (J + 16 >= LC) DY(J1) += 0.5*AKP*FR*ABB11*CA * y(J1 + 16);
					if (J + 9 >= LC)  DY(J1) += 0.5*AKP*ABB1*CA*y(J2 + 9);
					if (J + 10 >= LC) DY(J1) += 0.5*AKP*ABB2*CA*y(J2 + 10);
					if (J + 11 >= LC) DY(J1) += 0.5*AKP*ABB3*CA*y(J2 + 11);
					if (J + 12 >= LC) DY(J1) += 0.5*AKP*ABB4*CA*y(J2 + 12);
					if (J + 13 >= LC) DY(J1) += 0.5*AKP*ABB5*CA*y(J2 + 13);
					if (J + 14 >= LC) DY(J1) += 0.5*AKP*ABB6*CA*y(J2 + 14);
					if (J + 16 >= LC) DY(J1) += 0.5*AKP*ABB8*CA*y(J2 + 16);
					if (J + 17 >= LC) DY(J1) += 0.5*AKP*ABB9*CA*y(J2 + 17);
					if (J + 18 >= LC) DY(J1) += 0.5*AKP*ABB10*CA*y(J2 + 18);
					if (J + 19 >= LC) DY(J1) += 0.5*AKP*ABB11*CA*y(J2 + 19);


					DY(J2) = AKP * ((1. - FR)*SOM6 + SOM7
						+ 0.5*(ABB1*(SY2_12)
							+ABB2 * (SY2_13)
							+ABB3 * (SY2_14)
							+ABB4 * (SY2_15)
							+ABB5 * (SY2_16)
							+ABB6 * (SY2_17)
							+ABB7 * (SY2_18)
							+ABB8 * (SY2_19)
							+ABB9 * (SY2_20)
							+ABB10 * (SY2_21)
							+ABB11 * (SY2_22))
						+ 0.5*FR*(ABB1*SY2_1 + ABB2 * SY2_2 + ABB3 * SY2_3
							+ ABB4 * SY2_4 + ABB5 * SY2_5 + ABB6 * SY2_6
							+ ABB7 * SY2_7 + ABB8 * SY2_8 + ABB9 * SY2_9
							+ ABB10 * SY2_10 + ABB11 * SY2_11));


					if (J >= LC) DY(J2) = DY(J2) - AKP * ((J - 8) + CA + CD)*y(J2);
					if (J + 3 >= LC) DY(J2) += AKP * CD*y(J2 + 3);
					if (J + 6 >= LC) DY(J2) += 0.5*AKP*ABB1*CA*y(J2 + 6);
					if (J + 7 >= LC) DY(J2) += 0.5*AKP*ABB2*CA*y(J2 + 7);
					if (J + 8 >= LC) DY(J2) += 0.5*AKP*ABB3*CA*y(J2 + 8);
					if (J + 9 >= LC) DY(J2) += 0.5*AKP*ABB4*CA*y(J2 + 9);
					if (J + 10 >= LC) DY(J2) += 0.5*AKP*ABB5*CA*y(J2 + 10);
					if (J + 11 >= LC) DY(J2) += 0.5*AKP*ABB6*CA*y(J2 + 11);
					if (J + 12 >= LC) DY(J2) += 0.5*AKP*ABB7*CA*y(J2 + 12);
					if (J + 13 >= LC) DY(J2) += 0.5*AKP*ABB8*CA*y(J2 + 13);
					if (J + 14 >= LC) DY(J2) += 0.5*AKP*ABB9*CA*y(J2 + 14);
					if (J + 15 >= LC) DY(J2) += 0.5*AKP*ABB10*CA*y(J2 + 15);
					if (J + 16 >= LC) DY(J2) += 0.5*AKP*ABB11*CA*y(J2 + 16);

					const int aaa = 1;	// TODO
				}

				else if (J == 10)
				{
					DY(J - 1) = AKP * (
						SOM1*(1. - (ABB1 + ABB2 + ABB3 + ABB4 + ABB5 + ABB6))
						+ (1. - FR)*(SOM2)*(1. - (ABB1 + ABB2 + ABB3 + ABB4
							+ ABB5 + ABB6))
						+ 0.5*(ABB1*SY_1 + ABB2 * SY_2 + ABB3 * SY_3
							+ ABB4 * SY_4 + ABB5 * SY_5 + ABB6 * SY_6
							+ ABB7 * SY_7 + ABB8 * SY_8
							+ ABB10 * SY_10 + ABB11 * SY_11)
						+ 0.5*FR*(ABB1*(SY_12)
							+ABB2 * (SY_13)
							+ABB3 * (SY_14)
							+ABB4 * (SY_15)
							+ABB5 * (SY_16)
							+ABB6 * (SY_17)
							+ABB7 * (SY_18)
							+ABB8 * (SY_19)
							+ABB10 * (SY_21)
							+ABB11 * (SY_22))
						+ ABB9 * CORRP(3 + J - 1));


					if (J >= LC) DY(J - 1) += -AKP * (J - 4)*y(J - 1);
					if (J + 4 >= LC) DY(J - 1) += AKP * (1. - FR)*CA*y(J1 + 4)* (1 - (ATOT - ABB11 - ABB10 - ABB9 - ABB8 - ABB7));
					if (J + 9 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB1*CA*y(J1 + 9);
					if (J + 10 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB2*CA*y(J1 + 10);
					if (J + 11 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB3*CA*y(J1 + 11);
					if (J + 12 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB4*CA*y(J1 + 12);
					if (J + 13 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB5*CA*y(J1 + 13);
					if (J + 14 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB6*CA*y(J1 + 14);
					if (J + 15 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB7*CA*y(J1 + 15);
					if (J + 16 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB8*CA*y(J1 + 16);
					if (J + 18 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB10*CA*y(J1 + 18);
					if (J + 19 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB11*CA*y(J1 + 19);


					DY(J1) = AKP * (SOM3 +
						FR * (SOM4)
						+(SOM5)*(1. - (ABB1 + ABB2 + ABB3 + ABB4 + ABB5))
						+ FR * (SOM4)*(1. - (ABB1 + ABB2 + ABB3 + ABB4 + ABB5))
						+ 0.5*(ABB1*SY1_1 + ABB2 * SY1_2 + ABB3 * SY1_3
							+ ABB4 * SY1_4
							+ ABB5 * SY1_5 + ABB7 * SY1_7 + ABB8 * SY1_8
							+ ABB9 * SY1_9 + ABB10 * SY1_10 + ABB11 * SY1_11)
						+ 0.5*FR*(ABB1*(SY1_12 + SY12_1)
							+ ABB2 * (SY1_13 + SY12_2)
							+ ABB3 * (SY1_14 + SY12_3)
							+ ABB4 * (SY1_15 + SY12_4)
							+ ABB5 * (SY1_16 + SY12_5)
							+ ABB7 * (SY1_18 + SY12_7)
							+ ABB8 * (SY1_19 + SY12_8)
							+ ABB9 * (SY1_20 + SY12_9)
							+ ABB10 * (SY1_21 + SY12_10)
							+ ABB11 * (SY1_22 + SY12_11))
						+ 0.5*(ABB1*(SY1_23)
							+ABB2 * (SY1_24)
							+ABB3 * (SY1_25)
							+ABB4 * (SY1_26)
							+ABB5 * (SY1_27)
							+ABB7 * (SY1_29)
							+ABB8 * (SY1_30)
							+ABB9 * (SY1_31)
							+ABB10 * (SY1_32)
							+ABB11 * (SY1_33))
						+ ABB6 * CORRO(J - 1));


					if (J >= LC) DY(J1) = DY(J1) - AKP * ((J - 6) + FR * CD + (1. - FR) *CA)*y(J1);
					if (J + 3 >= LC) DY(J1) += AKP * FR*(CD*y(J1 + 3) + y(J1 + 3)*(1. - (ATOT - ABB11 - ABB10 - ABB9 - ABB8 - ABB7 - ABB6)));
					if (J + 4 >= LC) DY(J1) += AKP * CA*y(J2 + 4) *(1. - (ATOT - ABB11 - ABB10 - ABB9 - ABB8 - ABB7 - ABB6));
					if (J + 6 >= LC) DY(J1) += 0.5*AKP*FR*ABB1*CA*y(J1 + 6);
					if (J + 7 >= LC) DY(J1) += 0.5*AKP*FR*ABB2*CA*y(J1 + 7);
					if (J + 8 >= LC) DY(J1) += 0.5*AKP*FR*ABB3*CA*y(J1 + 8);
					if (J + 9 >= LC) DY(J1) += 0.5*AKP*FR*ABB4*CA*y(J1 + 9);
					if (J + 10 >= LC) DY(J1) += 0.5*AKP*FR*ABB5*CA*y(J1 + 10);
					if (J + 12 >= LC) DY(J1) += 0.5*AKP*FR*ABB7*CA * y(J1 + 12);
					if (J + 13 >= LC) DY(J1) += 0.5*AKP*FR*ABB8*CA *y(J1 + 13);
					if (J + 14 >= LC) DY(J1) += 0.5*AKP*FR*ABB9*CA *y(J1 + 14);
					if (J + 15 >= LC) DY(J1) += 0.5*AKP*FR*ABB10*CA *y(J1 + 15);
					if (J + 16 >= LC) DY(J1) += 0.5*AKP*FR*ABB11*CA *y(J1 + 16);
					if (J + 9 >= LC)  DY(J1) += 0.5*AKP*ABB1*CA*y(J2 + 9);
					if (J + 10 >= LC) DY(J1) += 0.5*AKP*ABB2*CA*y(J2 + 10);
					if (J + 11 >= LC) DY(J1) += 0.5*AKP*ABB3*CA*y(J2 + 11);
					if (J + 12 >= LC) DY(J1) += 0.5*AKP*ABB4*CA*y(J2 + 12);
					if (J + 13 >= LC) DY(J1) += 0.5*AKP*ABB5*CA*y(J2 + 13);
					if (J + 15 >= LC) DY(J1) += 0.5*AKP*ABB7*CA*y(J2 + 15);
					if (J + 16 >= LC) DY(J1) += 0.5*AKP*ABB8*CA*y(J2 + 16);
					if (J + 17 >= LC) DY(J1) += 0.5*AKP*ABB9*CA*y(J2 + 17);
					if (J + 18 >= LC) DY(J1) += 0.5*AKP*ABB10*CA*y(J2 + 18);
					if (J + 19 >= LC) DY(J1) += 0.5*AKP*ABB11*CA*y(J2 + 19);


					DY(J2) = AKP * ((1. - FR)*SOM6 + SOM7
						+ 0.5*(ABB1*(SY2_12)
							+ABB2 * (SY2_13)
							+ABB3 * (SY2_14)
							+ABB4 * (SY2_15)
							+ABB5 * (SY2_16)
							+ABB6 * (SY2_17)
							+ABB7 * (SY2_18)
							+ABB8 * (SY2_19)
							+ABB9 * (SY2_20)
							+ABB10 * (SY2_21)
							+ABB11 * (SY2_22))
						+ 0.5*FR*(ABB1*SY2_1 + ABB2 * SY2_2 + ABB3 * SY2_3
							+ ABB4 * SY2_4 + ABB5 * SY2_5 + ABB6 * SY2_6
							+ ABB7 * SY2_7 + ABB8 * SY2_8 + ABB9 * SY2_9
							+ ABB10 * SY2_10 + ABB11 * SY2_11));


					if (J >= LC) DY(J2) += -AKP * ((J - 8) + CA + CD)*y(J2);
					if (J + 3 >= LC) DY(J2) += AKP * CD*y(J2 + 3);
					if (J + 6 >= LC) DY(J2) += 0.5*AKP*ABB1*CA*y(J2 + 6);
					if (J + 7 >= LC) DY(J2) += 0.5*AKP*ABB2*CA*y(J2 + 7);
					if (J + 8 >= LC) DY(J2) += 0.5*AKP*ABB3*CA*y(J2 + 8);
					if (J + 9 >= LC) DY(J2) += 0.5*AKP*ABB4*CA*y(J2 + 9);
					if (J + 10 >= LC) DY(J2) += 0.5*AKP*ABB5*CA*y(J2 + 10);
					if (J + 11 >= LC) DY(J2) += 0.5*AKP*ABB6*CA*y(J2 + 11);
					if (J + 12 >= LC) DY(J2) += 0.5*AKP*ABB7*CA*y(J2 + 12);
					if (J + 13 >= LC) DY(J2) += 0.5*AKP*ABB8*CA*y(J2 + 13);
					if (J + 14 >= LC) DY(J2) += 0.5*AKP*ABB9*CA*y(J2 + 14);
					if (J + 15 >= LC) DY(J2) += 0.5*AKP*ABB10*CA*y(J2 + 15);
					if (J + 16 >= LC) DY(J2) += 0.5*AKP*ABB11*CA*y(J2 + 16);

					const int aaa = 1;	// TODO
				}

				else if (J == 9)
				{
					DY(J - 1) = AKP * (
						SOM1*(1. - (ABB1 + ABB2 + ABB3 + ABB4 + ABB5))
						+ (1. - FR)*(SOM2)*(1. - (ABB1 + ABB2 + ABB3 + ABB4
							+ ABB5))
						+ 0.5*(ABB1*SY_1 + ABB2 * SY_2 + ABB3 * SY_3
							+ ABB4 * SY_4 + ABB5 * SY_5 + ABB6 * SY_6
							+ ABB7 * SY_7 + ABB9 * SY_9
							+ ABB10 * SY_10 + ABB11 * SY_11)
						+ 0.5*FR*(ABB1*(SY_12)
							+ABB2 * (SY_13)
							+ABB3 * (SY_14)
							+ABB4 * (SY_15)
							+ABB5 * (SY_16)
							+ABB6 * (SY_17)
							+ABB7 * (SY_18)
							+ABB9 * (SY_20)
							+ABB10 * (SY_21)
							+ABB11 * (SY_22))
						+ ABB8 * CORRP(3 + J - 1));


					if (J >= LC) DY(J - 1) = DY(J - 1) - AKP * (J - 4)*y(J - 1);
					if (J + 4 >= LC) DY(J - 1) += AKP * (1. - FR)*CA*y(J1 + 4)* (1 - (ATOT - ABB11 - ABB10 - ABB9 - ABB8 - ABB7 - ABB6));
					if (J + 9 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB1*CA*y(J1 + 9);
					if (J + 10 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB2*CA*y(J1 + 10);
					if (J + 11 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB3*CA*y(J1 + 11);
					if (J + 12 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB4*CA*y(J1 + 12);
					if (J + 13 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB5*CA*y(J1 + 13);
					if (J + 14 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB6*CA*y(J1 + 14);
					if (J + 15 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB7*CA*y(J1 + 15);
					if (J + 17 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB9*CA*y(J1 + 17);
					if (J + 18 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB10*CA*y(J1 + 18);
					if (J + 19 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB11*CA*y(J1 + 19);


					DY(J1) = AKP * (SOM3 +
						FR * (SOM4)
						+(SOM5)*(1. - (ABB1 + ABB2 + ABB3 + ABB4))
						+ FR * (SOM4)*(1. - (ABB1 + ABB2 + ABB3 + ABB4))
						+ 0.5*(ABB1*SY1_1 + ABB2 * SY1_2 + ABB3 * SY1_3
							+ ABB4 * SY1_4
							+ ABB6 * SY1_6 + ABB7 * SY1_7 + ABB8 * SY1_8
							+ ABB9 * SY1_9 + ABB10 * SY1_10 + ABB11 * SY1_11)
						+ 0.5*FR*(ABB1*(SY1_12 + SY12_1)
							+ ABB2 * (SY1_13 + SY12_2)
							+ ABB3 * (SY1_14 + SY12_3)
							+ ABB4 * (SY1_15 + SY12_4)
							+ ABB6 * (SY1_17 + SY12_6)
							+ ABB7 * (SY1_18 + SY12_7)
							+ ABB8 * (SY1_19 + SY12_8)
							+ ABB9 * (SY1_20 + SY12_9)
							+ ABB10 * (SY1_21 + SY12_10)
							+ ABB11 * (SY1_22 + SY12_11))
						+ 0.5*(ABB1*(SY1_23)
							+ABB2 * (SY1_24)
							+ABB3 * (SY1_25)
							+ABB4 * (SY1_26)
							+ABB6 * (SY1_28)
							+ABB7 * (SY1_29)
							+ABB8 * (SY1_30)
							+ABB9 * (SY1_31)
							+ABB10 * (SY1_32)
							+ABB11 * (SY1_33))
						+ ABB5 * CORRO(J - 1));

					if (J >= LC) DY(J1) = DY(J1) - AKP * ((J - 6) + FR * CD + (1. - FR) *CA)*y(J1);
					if (J + 3 >= LC) DY(J1) += AKP * FR*(CD*y(J1 + 3) + y(J1 + 3)*(1. - (ABB1 + ABB2 + ABB3 + ABB4)));
					if (J + 4 >= LC) DY(J1) += AKP * CA*y(J2 + 4) *(1. - (ABB1 + ABB2 + ABB3 + ABB4));
					if (J + 6 >= LC) DY(J1) += 0.5*AKP*FR*ABB1*CA*y(J1 + 6);
					if (J + 7 >= LC) DY(J1) += 0.5*AKP*FR*ABB2*CA*y(J1 + 7);
					if (J + 8 >= LC) DY(J1) += 0.5*AKP*FR*ABB3*CA*y(J1 + 8);
					if (J + 9 >= LC) DY(J1) += 0.5*AKP*FR*ABB4*CA*y(J1 + 9);
					if (J + 11 >= LC) DY(J1) += 0.5*AKP*FR*ABB6*CA *y(J1 + 11);
					if (J + 12 >= LC) DY(J1) += 0.5*AKP*FR*ABB7*CA *y(J1 + 12);
					if (J + 13 >= LC) DY(J1) += 0.5*AKP*FR*ABB8*CA * y(J1 + 13);
					if (J + 14 >= LC) DY(J1) += 0.5*AKP*FR*ABB9*CA *y(J1 + 14);
					if (J + 15 >= LC) DY(J1) += 0.5*AKP*FR*ABB10*CA * y(J1 + 15);
					if (J + 16 >= LC) DY(J1) += 0.5*AKP*FR*ABB11*CA * y(J1 + 16);
					if (J + 9 >= LC)  DY(J1) += 0.5*AKP*ABB1*CA*y(J2 + 9);
					if (J + 10 >= LC) DY(J1) += 0.5*AKP*ABB2*CA*y(J2 + 10);
					if (J + 11 >= LC) DY(J1) += 0.5*AKP*ABB3*CA*y(J2 + 11);
					if (J + 12 >= LC) DY(J1) += 0.5*AKP*ABB4*CA*y(J2 + 12);
					if (J + 14 >= LC) DY(J1) += 0.5*AKP*ABB6*CA*y(J2 + 14);
					if (J + 15 >= LC) DY(J1) += 0.5*AKP*ABB7*CA*y(J2 + 15);
					if (J + 16 >= LC) DY(J1) += 0.5*AKP*ABB8*CA*y(J2 + 16);
					if (J + 17 >= LC) DY(J1) += 0.5*AKP*ABB9*CA*y(J2 + 17);
					if (J + 18 >= LC) DY(J1) += 0.5*AKP*ABB10*CA*y(J2 + 18);
					if (J + 19 >= LC) DY(J1) += 0.5*AKP*ABB11*CA*y(J2 + 19);


					DY(J2) = AKP * ((1. - FR)*SOM6 + SOM7
						+ 0.5*(ABB1*(SY2_12)
							+ABB2 * (SY2_13)
							+ABB3 * (SY2_14)
							+ABB4 * (SY2_15)
							+ABB5 * (SY2_16)
							+ABB6 * (SY2_17)
							+ABB7 * (SY2_18)
							+ABB8 * (SY2_19)
							+ABB9 * (SY2_20)
							+ABB10 * (SY2_21)
							+ABB11 * (SY2_22))
						+ 0.5*FR*(ABB1*SY2_1 + ABB2 * SY2_2 + ABB3 * SY2_3
							+ ABB4 * SY2_4 + ABB5 * SY2_5 + ABB6 * SY2_6
							+ ABB7 * SY2_7 + ABB8 * SY2_8 + ABB9 * SY2_9
							+ ABB10 * SY2_10 + ABB11 * SY2_11));


					if (J >= LC) DY(J2) = DY(J2) - AKP * ((J - 8) + CA + CD)*y(J2);
					if (J + 3 >= LC) DY(J2) += AKP * CD*y(J2 + 3);
					if (J + 6 >= LC) DY(J2) += 0.5*AKP*ABB1*CA*y(J2 + 6);
					if (J + 7 >= LC) DY(J2) += 0.5*AKP*ABB2*CA*y(J2 + 7);
					if (J + 8 >= LC) DY(J2) += 0.5*AKP*ABB3*CA*y(J2 + 8);
					if (J + 9 >= LC) DY(J2) += 0.5*AKP*ABB4*CA*y(J2 + 9);
					if (J + 10 >= LC) DY(J2) += 0.5*AKP*ABB5*CA*y(J2 + 10);
					if (J + 11 >= LC) DY(J2) += 0.5*AKP*ABB6*CA*y(J2 + 11);
					if (J + 12 >= LC) DY(J2) += 0.5*AKP*ABB7*CA*y(J2 + 12);
					if (J + 13 >= LC) DY(J2) += 0.5*AKP*ABB8*CA*y(J2 + 13);
					if (J + 14 >= LC) DY(J2) += 0.5*AKP*ABB9*CA*y(J2 + 14);
					if (J + 15 >= LC) DY(J2) += 0.5*AKP*ABB10*CA*y(J2 + 15);
					if (J + 16 >= LC) DY(J2) += 0.5*AKP*ABB11*CA*y(J2 + 16);

				}

				else if (J == 8)
				{
					DY(J - 1) = AKP * (
						SOM1*(1. - (ABB1 + ABB2 + ABB3 + ABB4))
						+ (1. - FR)*(SOM2)*(1. - (ABB1 + ABB2 + ABB3 + ABB4))
						+ 0.5*(ABB1*SY_1 + ABB2 * SY_2 + ABB3 * SY_3
							+ ABB4 * SY_4 + ABB5 * SY_5 + ABB6 * SY_6
							+ ABB8 * SY_8 + ABB9 * SY_9
							+ ABB10 * SY_10 + ABB11 * SY_11)
						+ 0.5*FR*(ABB1*(SY_12)
							+ABB2 * (SY_13)
							+ABB3 * (SY_14)
							+ABB4 * (SY_15)
							+ABB5 * (SY_16)
							+ABB6 * (SY_17)
							+ABB8 * (SY_19)
							+ABB9 * (SY_20)
							+ABB10 * (SY_21)
							+ABB11 * (SY_22))
						+ ABB7 * CORRP(3 + J - 1));


					if (J >= LC) DY(J - 1) = DY(J - 1) - AKP * (J - 4)*y(J - 1);
					if (J + 4 >= LC) DY(J - 1) += AKP * (1. - FR)*CA*y(J1 + 4)* (1 - (ABB1 + ABB2 + ABB3 + ABB4));
					if (J + 9 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB1*CA*y(J1 + 9);
					if (J + 10 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB2*CA*y(J1 + 10);
					if (J + 11 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB3*CA*y(J1 + 11);
					if (J + 12 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB4*CA*y(J1 + 12);
					if (J + 13 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB5*CA*y(J1 + 13);
					if (J + 14 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB6*CA*y(J1 + 14);
					if (J + 16 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB8*CA*y(J1 + 16);
					if (J + 17 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB9*CA*y(J1 + 17);
					if (J + 18 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB10*CA*y(J1 + 18);
					if (J + 19 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB11*CA*y(J1 + 19);



					DY(J1) = AKP * (SOM3 +
						FR * (SOM4)
						+(SOM5)*(1. - (ABB1 + ABB2 + ABB3))
						+ FR * (SOM4)*(1. - (ABB1 + ABB2 + ABB3))
						+ 0.5*(ABB1*SY1_1 + ABB2 * SY1_2 + ABB3 * SY1_3
							+ ABB5 * SY1_5
							+ ABB6 * SY1_6 + ABB7 * SY1_7 + ABB8 * SY1_8
							+ ABB9 * SY1_9 + ABB10 * SY1_10 + ABB11 * SY1_11)
						+ 0.5*FR*(ABB1*(SY1_12 + SY12_1)
							+ ABB2 * (SY1_13 + SY12_2)
							+ ABB3 * (SY1_14 + SY12_3)
							+ ABB5 * (SY1_16 + SY12_5)
							+ ABB6 * (SY1_17 + SY12_6)
							+ ABB7 * (SY1_18 + SY12_7)
							+ ABB8 * (SY1_19 + SY12_8)
							+ ABB9 * (SY1_20 + SY12_9)
							+ ABB10 * (SY1_21 + SY12_10)
							+ ABB11 * (SY1_22 + SY12_11))
						+ 0.5*(ABB1*(SY1_23)
							+ABB2 * (SY1_24)
							+ABB3 * (SY1_25)
							+ABB5 * (SY1_27)
							+ABB6 * (SY1_28)
							+ABB7 * (SY1_29)
							+ABB8 * (SY1_30)
							+ABB9 * (SY1_31)
							+ABB10 * (SY1_32)
							+ABB11 * (SY1_33))
						+ ABB4 * CORRO(J - 1));


					if (J >= LC) DY(J1) = DY(J1) - AKP * ((J - 6) + FR * CD + (1. - FR) *CA)*y(J1);
					if (J + 3 >= LC) DY(J1) += AKP * FR*(CD*y(J1 + 3) + y(J1 + 3)*(1. - (ABB1 + ABB2 + ABB3)));
					if (J + 4 >= LC) DY(J1) += AKP * CA*y(J2 + 4) * (1. - (ABB1 + ABB2 + ABB3));
					if (J + 6 >= LC) DY(J1) += 0.5*AKP*FR*ABB1*CA*y(J1 + 6);
					if (J + 7 >= LC) DY(J1) += 0.5*AKP*FR*ABB2*CA*y(J1 + 7);
					if (J + 8 >= LC) DY(J1) += 0.5*AKP*FR*ABB3*CA*y(J1 + 8);
					if (J + 10 >= LC) DY(J1) += 0.5*AKP*FR*ABB5*CA*y(J1 + 10);
					if (J + 11 >= LC) DY(J1) += 0.5*AKP*FR*ABB6*CA *y(J1 + 11);
					if (J + 12 >= LC) DY(J1) += 0.5*AKP*FR*ABB7*CA *y(J1 + 12);
					if (J + 13 >= LC) DY(J1) += 0.5*AKP*FR*ABB8*CA *y(J1 + 13);
					if (J + 14 >= LC) DY(J1) += 0.5*AKP*FR*ABB9*CA * y(J1 + 14);
					if (J + 15 >= LC) DY(J1) += 0.5*AKP*FR*ABB10*CA * y(J1 + 15);
					if (J + 16 >= LC) DY(J1) += 0.5*AKP*FR*ABB11*CA * y(J1 + 16);
					if (J + 9 >= LC)  DY(J1) += 0.5*AKP*ABB1*CA*y(J2 + 9);
					if (J + 10 >= LC) DY(J1) += 0.5*AKP*ABB2*CA*y(J2 + 10);
					if (J + 11 >= LC) DY(J1) += 0.5*AKP*ABB3*CA*y(J2 + 11);
					if (J + 13 >= LC) DY(J1) += 0.5*AKP*ABB5*CA*y(J2 + 13);
					if (J + 14 >= LC) DY(J1) += 0.5*AKP*ABB6*CA*y(J2 + 14);
					if (J + 15 >= LC) DY(J1) += 0.5*AKP*ABB7*CA*y(J2 + 15);
					if (J + 16 >= LC) DY(J1) += 0.5*AKP*ABB8*CA*y(J2 + 16);
					if (J + 17 >= LC) DY(J1) += 0.5*AKP*ABB9*CA*y(J2 + 17);
					if (J + 18 >= LC) DY(J1) += 0.5*AKP*ABB10*CA*y(J2 + 18);
					if (J + 19 >= LC) DY(J1) += 0.5*AKP*ABB11*CA*y(J2 + 19);




					DY(J2) = AKP * (1. / 2.*SOM6 + SOM7
						+ 0.5*(ABB1*(SY2_12)
							+ABB2 * (SY2_13)
							+ABB3 * (SY2_14)
							+ABB4 * (SY2_15)
							+ABB5 * (SY2_16)
							+ABB6 * (SY2_17)
							+ABB7 * (SY2_18)
							+ABB8 * (SY2_19)
							+ABB9 * (SY2_20)
							+ABB10 * (SY2_21)
							+ABB11 * (SY2_22))
						+ 0.5*FR*(ABB1*SY2_1 + ABB2 * SY2_2 + ABB3 * SY2_3
							+ ABB4 * SY2_4 + ABB5 * SY2_5 + ABB6 * SY2_6
							+ ABB7 * SY2_7 + ABB8 * SY2_8 + ABB9 * SY2_9
							+ ABB10 * SY2_10 + ABB11 * SY2_11));


					if (J >= LC) DY(J2) = DY(J2) - AKP * ((J - 8) + CA + CD)*y(J2);
					if (J + 3 >= LC) DY(J2) += AKP * CD*y(J2 + 3);
					if (J + 6 >= LC) DY(J2) += 0.5*AKP*ABB1*CA*y(J2 + 6);
					if (J + 7 >= LC) DY(J2) += 0.5*AKP*ABB2*CA*y(J2 + 7);
					if (J + 8 >= LC) DY(J2) += 0.5*AKP*ABB3*CA*y(J2 + 8);
					if (J + 9 >= LC) DY(J2) += 0.5*AKP*ABB4*CA*y(J2 + 9);
					if (J + 10 >= LC) DY(J2) += 0.5*AKP*ABB5*CA*y(J2 + 10);
					if (J + 11 >= LC) DY(J2) += 0.5*AKP*ABB6*CA*y(J2 + 11);
					if (J + 12 >= LC) DY(J2) += 0.5*AKP*ABB7*CA*y(J2 + 12);
					if (J + 13 >= LC) DY(J2) += 0.5*AKP*ABB8*CA*y(J2 + 13);
					if (J + 14 >= LC) DY(J2) += 0.5*AKP*ABB9*CA*y(J2 + 14);
					if (J + 15 >= LC) DY(J2) += 0.5*AKP*ABB10*CA*y(J2 + 15);
					if (J + 16 >= LC) DY(J2) += 0.5*AKP*ABB11*CA*y(J2 + 16);

				}

				else if (J == 7)
				{
					DY(J - 1) = AKP * (
						SOM1*(1. - (ABB1 + ABB2 + ABB3))
						+ (1. - FR)*(SOM2)*(1. - (ABB1 + ABB2 + ABB3))
						+ 0.5*(ABB1*SY_1 + ABB2 * SY_2 + ABB3 * SY_3
							+ ABB4 * SY_4 + ABB5 * SY_5
							+ ABB7 * SY_7 + ABB8 * SY_8 + ABB9 * SY_9
							+ ABB10 * SY_10 + ABB11 * SY_11)
						+ 0.5*FR*(ABB1*(SY_12)
							+ABB2 * (SY_13)
							+ABB3 * (SY_14)
							+ABB4 * (SY_15)
							+ABB5 * (SY_16)
							+ABB7 * (+SY_18)
							+ ABB8 * (SY_19)
							+ABB9 * (SY_20)
							+ABB10 * (SY_21)
							+ABB11 * (SY_22))
						+ ABB6 * CORRP(3 + J - 1));


					if (J >= LC) DY(J - 1) = DY(J - 1) - AKP * (J - 4)*y(J - 1);
					if (J + 4 >= LC) DY(J - 1) += AKP * (1. - FR)*CA*y(J1 + 4)* (1 - (ABB1 + ABB2 + ABB3));
					if (J + 9 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB1*CA*y(J1 + 9);
					if (J + 10 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB2*CA*y(J1 + 10);
					if (J + 11 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB3*CA*y(J1 + 11);
					if (J + 12 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB4*CA*y(J1 + 12);
					if (J + 13 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB5*CA*y(J1 + 13);
					if (J + 15 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB7*CA*y(J1 + 15);
					if (J + 16 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB8*CA*y(J1 + 16);
					if (J + 17 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB9*CA*y(J1 + 17);
					if (J + 18 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB10*CA*y(J1 + 18);
					if (J + 19 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB11*CA*y(J1 + 19);


					DY(J1) = AKP * (SOM3 +
						FR * (SOM4)
						+(SOM5)*(1. - (ABB1 + ABB2))
						+ FR * (SOM4)*(1. - (ABB1 + ABB2))
						+ 0.5*(ABB1*SY1_1 + ABB2 * SY1_2 + ABB4 * SY1_4
							+ ABB5 * SY1_5
							+ ABB6 * SY1_6 + ABB7 * SY1_7 + ABB8 * SY1_8
							+ ABB9 * SY1_9 + ABB10 * SY1_10 + ABB11 * SY1_11)
						+ 0.5*FR*(ABB1*(SY1_12 + SY12_1)
							+ ABB2 * (SY1_13 + SY12_2)
							+ ABB4 * (SY1_15 + SY12_4)
							+ ABB5 * (SY1_16 + SY12_5)
							+ ABB6 * (SY1_17 + SY12_6)
							+ ABB7 * (SY1_18 + SY12_7)
							+ ABB8 * (SY1_19 + SY12_8)
							+ ABB9 * (SY1_20 + SY12_9)
							+ ABB10 * (SY1_21 + SY12_10)
							+ ABB11 * (SY1_22 + SY12_11))
						+ 0.5*(ABB1*(SY1_23)
							+ABB2 * (SY1_24)
							+ABB4 * (SY1_26)
							+ABB5 * (SY1_27)
							+ABB6 * (SY1_28)
							+ABB7 * (SY1_29)
							+ABB8 * (SY1_30)
							+ABB9 * (SY1_31)
							+ABB10 * (SY1_32)
							+ABB11 * (SY1_33))
						+ ABB3 * CORRO(J - 1));


					if (J >= LC) DY(J1) = DY(J1) - AKP * ((J - 6) + FR * CD + (1. - FR) * CA)*y(J1);
					if (J + 3 >= LC) DY(J1) += AKP * FR*(CD*y(J1 + 3) + y(J1 + 3)*(1. - (ABB1 + ABB2)));
					if (J + 4 >= LC) DY(J1) += AKP * CA*y(J2 + 4) *(1. - (ABB1 + ABB2));
					if (J + 6 >= LC) DY(J1) += 0.5*AKP*FR*ABB1*CA*y(J1 + 6);
					if (J + 7 >= LC) DY(J1) += 0.5*AKP*FR*ABB2*CA*y(J1 + 7);
					if (J + 9 >= LC) DY(J1) += 0.5*AKP*FR*ABB4*CA*y(J1 + 9);
					if (J + 10 >= LC) DY(J1) += 0.5*AKP*FR*ABB5*CA*y(J1 + 10);
					if (J + 11 >= LC) DY(J1) += 0.5*AKP*FR*ABB6*CA *y(J1 + 11);
					if (J + 12 >= LC) DY(J1) += 0.5*AKP*FR*ABB7*CA * y(J1 + 12);
					if (J + 13 >= LC) DY(J1) += 0.5*AKP*FR*ABB8*CA *y(J1 + 13);
					if (J + 14 >= LC) DY(J1) += 0.5*AKP*FR*ABB9*CA * y(J1 + 14);
					if (J + 15 >= LC) DY(J1) += 0.5*AKP*FR*ABB10*CA * y(J1 + 15);
					if (J + 16 >= LC) DY(J1) += 0.5*AKP*FR*ABB11*CA * y(J1 + 16);
					if (J + 9 >= LC)  DY(J1) += 0.5*AKP*ABB1*CA*y(J2 + 9);
					if (J + 10 >= LC) DY(J1) += 0.5*AKP*ABB2*CA*y(J2 + 10);
					if (J + 12 >= LC) DY(J1) += 0.5*AKP*ABB4*CA*y(J2 + 12);
					if (J + 13 >= LC) DY(J1) += 0.5*AKP*ABB5*CA*y(J2 + 13);
					if (J + 14 >= LC) DY(J1) += 0.5*AKP*ABB6*CA*y(J2 + 14);
					if (J + 15 >= LC) DY(J1) += 0.5*AKP*ABB7*CA*y(J2 + 15);
					if (J + 16 >= LC) DY(J1) += 0.5*AKP*ABB8*CA*y(J2 + 16);
					if (J + 17 >= LC) DY(J1) += 0.5*AKP*ABB9*CA*y(J2 + 17);
					if (J + 18 >= LC) DY(J1) += 0.5*AKP*ABB10*CA*y(J2 + 18);
					if (J + 19 >= LC) DY(J1) += 0.5*AKP*ABB11*CA*y(J2 + 19);

					DY(J2) = AKP * (1. / 2.*SOM6 + SOM7
						+ 0.5*(ABB1*(SY2_12)
							+ABB2 * (SY2_13)
							+ABB3 * (SY2_14)
							+ABB4 * (SY2_15)
							+ABB5 * (SY2_16)
							+ABB6 * (SY2_17)
							+ABB7 * (SY2_18)
							+ABB8 * (SY2_19)
							+ABB9 * (SY2_20)
							+ABB10 * (SY2_21)
							+ABB11 * (SY2_22))
						+ 0.5*FR*(ABB1*SY2_1 + ABB2 * SY2_2 + ABB3 * SY2_3
							+ ABB4 * SY2_4 + ABB5 * SY2_5 + ABB6 * SY2_6
							+ ABB7 * SY2_7 + ABB8 * SY2_8 + ABB9 * SY2_9
							+ ABB10 * SY2_10 + ABB11 * SY2_11));

					if (J >= LC) DY(J2) = DY(J2) - AKP * CA*CD*y(J2);
					if (J + 3 >= LC) DY(J2) += AKP * CD*y(J2 + 3);
					if (J + 6 >= LC) DY(J2) += 0.5*AKP*ABB1*CA*y(J2 + 6);
					if (J + 7 >= LC) DY(J2) += 0.5*AKP*ABB2*CA*y(J2 + 7);
					if (J + 8 >= LC) DY(J2) += 0.5*AKP*ABB3*CA*y(J2 + 8);
					if (J + 9 >= LC) DY(J2) += 0.5*AKP*ABB4*CA*y(J2 + 9);
					if (J + 10 >= LC) DY(J2) += 0.5*AKP*ABB5*CA*y(J2 + 10);
					if (J + 11 >= LC) DY(J2) += 0.5*AKP*ABB6*CA*y(J2 + 11);
					if (J + 12 >= LC) DY(J2) += 0.5*AKP*ABB7*CA*y(J2 + 12);
					if (J + 13 >= LC) DY(J2) += 0.5*AKP*ABB8*CA*y(J2 + 13);
					if (J + 14 >= LC) DY(J2) += 0.5*AKP*ABB9*CA*y(J2 + 14);
					if (J + 15 >= LC) DY(J2) += 0.5*AKP*ABB10*CA*y(J2 + 15);
					if (J + 16 >= LC) DY(J2) += 0.5*AKP*ABB11*CA*y(J2 + 16);
				}

				else if (J == 6)
				{
					DY(J - 1) = AKP * (
						SOM1*(1. - (ABB1 + ABB2))
						+ (1. - FR)*(SOM2)*(1. - (ABB1 + ABB2))
						+ 0.5*(ABB1*SY_1 + ABB2 * SY_2 + ABB3 * SY_3
							+ ABB4 * SY_4 + ABB6 * SY_6
							+ ABB7 * SY_7 + ABB8 * SY_8 + ABB9 * SY_9
							+ ABB10 * SY_10 + ABB11 * SY_11)
						+ 0.5*FR*(ABB1*(SY_12)
							+ABB2 * (SY_13)
							+ABB3 * (SY_14)
							+ABB4 * (SY_15)
							+ABB6 * (SY_17)
							+ABB7 * (SY_18)
							+ABB8 * (SY_19)
							+ABB9 * (SY_20)
							+ABB10 * (SY_21)
							+ABB11 * (SY_22))
						+ ABB5 * CORRP(3 + J - 1));


					if (J >= LC) DY(J - 1) = DY(J - 1) - AKP * (J - 4)*y(J - 1);
					if (J + 4 >= LC) DY(J - 1) += AKP * (1. - FR)*CA*y(J1 + 4) * (1 - (ABB1 + ABB2));
					if (J + 9 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB1*CA*y(J1 + 9);
					if (J + 10 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB2*CA*y(J1 + 10);
					if (J + 11 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB3*CA*y(J1 + 11);
					if (J + 12 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB4*CA*y(J1 + 12);
					if (J + 14 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB6*CA*y(J1 + 14);
					if (J + 15 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB7*CA*y(J1 + 15);
					if (J + 16 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB8*CA*y(J1 + 16);
					if (J + 17 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB9*CA*y(J1 + 17);
					if (J + 18 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB10*CA*y(J1 + 18);
					if (J + 19 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB11*CA*y(J1 + 19);


					DY(J1) = AKP * (SOM3 +
						FR * (SOM4)
						+(SOM5)*(1. - ABB1)
						+ FR * (SOM4)*(1. - ABB1)
						+ 0.5*(ABB1*SY1_1 + ABB3 * SY1_3 + ABB4 * SY1_4
							+ ABB5 * SY1_5
							+ ABB6 * SY1_6 + ABB7 * SY1_7 + ABB8 * SY1_8
							+ ABB9 * SY1_9 + ABB10 * SY1_10 + ABB11 * SY1_11)
						+ 0.5*FR*(ABB1*(SY1_12 + SY12_1)
							+ ABB3 * (SY1_14 + SY12_3)
							+ ABB4 * (SY1_15 + SY12_4)
							+ ABB5 * (SY1_16 + SY12_5)
							+ ABB6 * (SY1_17 + SY12_6)
							+ ABB7 * (SY1_18 + SY12_7)
							+ ABB8 * (SY1_19 + SY12_8)
							+ ABB9 * (SY1_20 + SY12_9)
							+ ABB10 * (SY1_21 + SY12_10)
							+ ABB11 * (SY1_22 + SY12_11))
						+ 0.5*(ABB1*(SY1_23)
							+ABB3 * (SY1_25)
							+ABB4 * (SY1_26)
							+ABB5 * (SY1_27)
							+ABB6 * (SY1_28)
							+ABB7 * (SY1_29)
							+ABB8 * (SY1_30)
							+ABB9 * (SY1_31)
							+ABB10 * (SY1_32)
							+ABB11 * (SY1_33))
						+ ABB2 * CORRO(J - 1));


					if (J >= LC) DY(J1) = DY(J1) - AKP * ((J - 6) + FR * CD + (1. - FR) *CA)*y(J1);
					if (J + 3 >= LC) DY(J1) += AKP * FR*(CD*y(J1 + 3) + y(J1 + 3)*(1. - (ABB1)));
					if (J + 4 >= LC) DY(J1) += AKP * CA*y(J2 + 4) *(1. - (ABB1));
					if (J + 6 >= LC) DY(J1) += 0.5*AKP*FR*ABB1*CA*y(J1 + 6);
					if (J + 8 >= LC) DY(J1) += 0.5*AKP*FR*ABB3*CA*y(J1 + 8);
					if (J + 9 >= LC) DY(J1) += 0.5*AKP*FR*ABB4*CA*y(J1 + 9);
					if (J + 10 >= LC) DY(J1) += 0.5*AKP*FR*ABB5*CA*y(J1 + 10);
					if (J + 11 >= LC) DY(J1) += 0.5*AKP*FR*ABB6*CA *y(J1 + 11);
					if (J + 12 >= LC) DY(J1) += 0.5*AKP*FR*ABB7*CA *y(J1 + 12);
					if (J + 13 >= LC) DY(J1) += 0.5*AKP*FR*ABB8*CA * y(J1 + 13);
					if (J + 14 >= LC) DY(J1) += 0.5*AKP*FR*ABB9*CA * y(J1 + 14);
					if (J + 15 >= LC) DY(J1) += 0.5*AKP*FR*ABB10*CA * y(J1 + 15);
					if (J + 16 >= LC) DY(J1) += 0.5*AKP*FR*ABB11*CA * y(J1 + 16);
					if (J + 9 >= LC)  DY(J1) += 0.5*AKP*ABB1*CA*y(J2 + 9);
					if (J + 11 >= LC) DY(J1) += 0.5*AKP*ABB3*CA*y(J2 + 11);
					if (J + 12 >= LC) DY(J1) += 0.5*AKP*ABB4*CA*y(J2 + 12);
					if (J + 13 >= LC) DY(J1) += 0.5*AKP*ABB5*CA*y(J2 + 13);
					if (J + 14 >= LC) DY(J1) += 0.5*AKP*ABB6*CA*y(J2 + 14);
					if (J + 15 >= LC) DY(J1) += 0.5*AKP*ABB7*CA*y(J2 + 15);
					if (J + 16 >= LC) DY(J1) += 0.5*AKP*ABB8*CA*y(J2 + 16);
					if (J + 17 >= LC) DY(J1) += 0.5*AKP*ABB9*CA*y(J2 + 17);
					if (J + 18 >= LC) DY(J1) += 0.5*AKP*ABB10*CA*y(J2 + 18);
					if (J + 19 >= LC) DY(J1) += 0.5*AKP*ABB11*CA*y(J2 + 19);

					DY(J2) = AKP * (1. / 2.*SOM6 + SOM7
						+ 1. / (1. + CD)*(ABB1*(SY2_12)
							+ABB2 * (SY2_13)
							+ABB3 * (SY2_14)
							+ABB4 * (SY2_15)
							+ABB5 * (SY2_16)
							+ABB6 * (SY2_17)
							+ABB7 * (SY2_18)
							+ABB8 * (SY2_19)
							+ABB9 * (SY2_20)
							+ABB10 * (SY2_21)
							+ABB11 * (SY2_22))
						+ 1. / (1. + CD)*FR*(ABB1*SY2_1 + ABB2 * SY2_2 + ABB3 * SY2_3
							+ ABB4 * SY2_4 + ABB5 * SY2_5 + ABB6 * SY2_6
							+ ABB7 * SY2_7 + ABB8 * SY2_8 + ABB9 * SY2_9
							+ ABB10 * SY2_10 + ABB11 * SY2_11));


					if (J + 3 >= LC) DY(J2) += AKP * CD*y(J2 + 3);
					if (J + 6 >= LC) DY(J2) += 1. / (1. + CD)*AKP*ABB1*CA *y(J2 + 6);
					if (J + 7 >= LC) DY(J2) += 1. / (1. + CD)*AKP*ABB2*CA *y(J2 + 7);
					if (J + 8 >= LC) DY(J2) += 1. / (1. + CD)*AKP*ABB3*CA * y(J2 + 8);
					if (J + 9 >= LC) DY(J2) += 1. / (1. + CD)*AKP*ABB4*CA * y(J2 + 9);
					if (J + 10 >= LC) DY(J2) += 1. / (1. + CD)*AKP*ABB5*CA *y(J2 + 10);
					if (J + 11 >= LC) DY(J2) += 1. / (1. + CD)*AKP*ABB6*CA * y(J2 + 11);
					if (J + 12 >= LC) DY(J2) += 1. / (1. + CD)*AKP*ABB7*CA * y(J2 + 12);
					if (J + 13 >= LC) DY(J2) += 1. / (1. + CD)*AKP*ABB8*CA *y(J2 + 13);
					if (J + 14 >= LC) DY(J2) += 1. / (1. + CD)*AKP*ABB9*CA *y(J2 + 14);
					if (J + 15 >= LC) DY(J2) += 1. / (1. + CD)*AKP*ABB10*CA * y(J2 + 15);
					if (J + 16 >= LC) DY(J2) += 1. / (1. + CD)*AKP*ABB11*CA * y(J2 + 16);

				}

				else if (J == 5)
				{
					DY(J - 1) = AKP * (
						SOM1*(1. - ABB1)
						+ (1. - FR)*(SOM2)*(1. - ABB1)
						+ 0.5*(ABB1*SY_1 + ABB2 * SY_2 + ABB3 * SY_3
							+ ABB5 * SY_5 + ABB6 * SY_6
							+ ABB7 * SY_7 + ABB8 * SY_8 + ABB9 * SY_9
							+ ABB10 * SY_10 + ABB11 * SY_11)
						+ 0.5*FR*(ABB1*(SY_12)
							+ABB2 * (SY_13)
							+ABB3 * (SY_14)
							+ABB5 * (SY_16)
							+ABB6 * (SY_17)
							+ABB7 * (SY_18)
							+ABB8 * (SY_19)
							+ABB9 * (SY_20)
							+ABB10 * (SY_21)
							+ABB11 * (SY_22))
						+ ABB4 * CORRP(3 + J - 1));

					if (J >= LC) DY(J - 1) = DY(J - 1) - AKP * (J - 4)*y(J - 1);
					if (J + 4 >= LC) DY(J - 1) += AKP * (1. - FR)*CA*y(J1 + 4)* (1 - ABB1);
					if (J + 9 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB1*CA*y(J1 + 9);
					if (J + 10 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB2*CA*y(J1 + 10);
					if (J + 11 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB3*CA*y(J1 + 11);
					if (J + 13 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB5*CA*y(J1 + 13);
					if (J + 14 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB6*CA*y(J1 + 14);
					if (J + 15 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB7*CA*y(J1 + 15);
					if (J + 16 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB8*CA*y(J1 + 16);
					if (J + 17 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB9*CA*y(J1 + 17);
					if (J + 18 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB10*CA*y(J1 + 18);
					if (J + 19 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB11*CA*y(J1 + 19);


					DY(J1) = AKP * (SOM3 + FR * (SOM4)
							+(SOM5)+FR * (SOM4)
							+0.5*(ABB2*SY1_2 + ABB3 * SY1_3 + ABB4 * SY1_4
							+ ABB5 * SY1_5
							+ ABB6 * SY1_6 + ABB7 * SY1_7 + ABB8 * SY1_8
							+ ABB9 * SY1_9 + ABB10 * SY1_10 + ABB11 * SY1_11)
							+ 0.5*FR*(ABB2*(SY1_13 + SY12_2)
							+ ABB3 * (SY1_14 + SY12_3)
							+ ABB4 * (SY1_15 + SY12_4)
							+ ABB5 * (SY1_16 + SY12_5)
							+ ABB6 * (SY1_17 + SY12_6)
							+ ABB7 * (SY1_18 + SY12_7)
							+ ABB8 * (SY1_19 + SY12_8)
							+ ABB9 * (SY1_20 + SY12_9)
							+ ABB10 * (SY1_21 + SY12_10)
							+ ABB11 * (SY1_22 + SY12_11))
							+ 0.5*(ABB2*(SY1_24)
							+ABB3 * (SY1_25)
							+ABB4 * (SY1_26)
							+ABB5 * (SY1_27)
							+ABB6 * (SY1_28)
							+ABB7 * (SY1_29)
							+ABB8 * (SY1_30)
							+ABB9 * (SY1_31)
							+ABB10 * (SY1_32)
							+ABB11 * (SY1_33))
						+ ABB1 * CORRO(J - 1));

					if (J + 3 >= LC) DY(J1) += AKP * FR*(CD*y(J1 + 3) + y(J1 + 3));
					if (J + 4 >= LC) DY(J1) += AKP * CA*y(J2 + 4);
					if (J + 7 >= LC) DY(J1) += 0.5*AKP*FR*ABB2*CA*y(J1 + 7);
					if (J + 8 >= LC) DY(J1) += 0.5*AKP*FR*ABB3*CA*y(J1 + 8);
					if (J + 9 >= LC) DY(J1) += 0.5*AKP*FR*ABB4*CA*y(J1 + 9);
					if (J + 10 >= LC) DY(J1) += 0.5*AKP*FR*ABB5*CA*y(J1 + 10);
					if (J + 11 >= LC) DY(J1) += 0.5*AKP*FR*ABB6*CA * y(J1 + 11);
					if (J + 12 >= LC) DY(J1) += 0.5*AKP*FR*ABB7*CA *y(J1 + 12);
					if (J + 13 >= LC) DY(J1) += 0.5*AKP*FR*ABB8*CA * y(J1 + 13);
					if (J + 14 >= LC) DY(J1) += 0.5*AKP*FR*ABB9*CA *  y(J1 + 14);
					if (J + 15 >= LC) DY(J1) += 0.5*AKP*FR*ABB10*CA * y(J1 + 15);
					if (J + 16 >= LC) DY(J1) += 0.5*AKP*FR*ABB11*CA *   y(J1 + 16);
					if (J + 10 >= LC)  DY(J1) += 0.5*AKP*ABB2*CA*y(J2 + 10);
					if (J + 11 >= LC) DY(J1) += 0.5*AKP*ABB3*CA*y(J2 + 11);
					if (J + 12 >= LC) DY(J1) += 0.5*AKP*ABB4*CA*y(J2 + 12);
					if (J + 13 >= LC) DY(J1) += 0.5*AKP*ABB5*CA*y(J2 + 13);
					if (J + 14 >= LC) DY(J1) += 0.5*AKP*ABB6*CA*y(J2 + 14);
					if (J + 15 >= LC) DY(J1) += 0.5*AKP*ABB7*CA*y(J2 + 15);
					if (J + 16 >= LC) DY(J1) += 0.5*AKP*ABB8*CA*y(J2 + 16);
					if (J + 17 >= LC) DY(J1) += 0.5*AKP*ABB9*CA*y(J2 + 17);
					if (J + 18 >= LC) DY(J1) += 0.5*AKP*ABB10*CA*y(J2 + 18);
					if (J + 19 >= LC) DY(J1) += 0.5*AKP*ABB11*CA*y(J2 + 19);


					DY(J2) = AKP * (1. / 2.*SOM6 + SOM7
							+ (ABB1*(SY2_12)
							+ABB2 * (SY2_13)
							+ABB3 * (SY2_14)
							+ABB4 * (SY2_15)
							+ABB5 * (SY2_16)
							+ABB6 * (SY2_17)
							+ABB7 * (SY2_18)
							+ABB8 * (SY2_19)
							+ABB9 * (SY2_20)
							+ABB10 * (SY2_21)
							+ABB11 * (SY2_22))
							+ FR * (ABB1*SY2_1 + ABB2 * SY2_2 + ABB3 * SY2_3
							+ ABB4 * SY2_4 + ABB5 * SY2_5 + ABB6 * SY2_6
							+ ABB7 * SY2_7 + ABB8 * SY2_8 + ABB9 * SY2_9
							+ ABB10 * SY2_10 + ABB11 * SY2_11));

					if (J + 3 >= LC) DY(J2) += AKP * CD*y(J2 + 3);
					if (J + 6 >= LC) DY(J2) += AKP * ABB1*CA*y(J2 + 6);
					if (J + 7 >= LC) DY(J2) += AKP * ABB2*CA*y(J2 + 7);
					if (J + 8 >= LC) DY(J2) += AKP * ABB3*CA*y(J2 + 8);
					if (J + 9 >= LC) DY(J2) += AKP * ABB4*CA*y(J2 + 9);
					if (J + 10 >= LC) DY(J2) += AKP * ABB5*CA*y(J2 + 10);
					if (J + 11 >= LC) DY(J2) += AKP * ABB6*CA*y(J2 + 11);
					if (J + 12 >= LC) DY(J2) += AKP * ABB7*CA*y(J2 + 12);
					if (J + 13 >= LC) DY(J2) += AKP * ABB8*CA*y(J2 + 13);
					if (J + 14 >= LC) DY(J2) += AKP * ABB9*CA*y(J2 + 14);
					if (J + 15 >= LC) DY(J2) += AKP * ABB10*CA*y(J2 + 15);
					if (J + 16 >= LC) DY(J2) += AKP * ABB11*CA*y(J2 + 16);

				}

				else if (J == 4)
				{
					DY(J - 1) = AKP * (SOM1 + (1. - FR)*(SOM2)
						+0.5*(ABB1*SY_1 + ABB2 * SY_2
							+ ABB4 * SY_4 + ABB5 * SY_5 + ABB6 * SY_6
							+ ABB7 * SY_7 + ABB8 * SY_8 + ABB9 * SY_9
							+ ABB10 * SY_10 + ABB11 * SY_11)
						+ 0.5*FR*(ABB1*(SY_12)
							+ABB2 * (SY_13)
							+ABB4 * (SY_15)
							+ABB5 * (SY_16)
							+ABB6 * (SY_17)
							+ABB7 * (SY_18)
							+ABB8 * (SY_19)
							+ABB9 * (SY_20)
							+ABB10 * (SY_21)
							+ABB11 * (SY_22))
						+ ABB3 * CORRP(3 + J - 1));

					if (J + 4 >= LC) DY(J - 1) += AKP * (1. - FR)*CA*y(J1 + 4);
					if (J + 9 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB1*CA*y(J1 + 9);
					if (J + 10 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB2*CA*y(J1 + 10);
					if (J + 12 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB4*CA*y(J1 + 12);
					if (J + 13 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB5*CA*y(J1 + 13);
					if (J + 14 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB6*CA*y(J1 + 14);
					if (J + 15 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB7*CA*y(J1 + 15);
					if (J + 16 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB8*CA*y(J1 + 16);
					if (J + 17 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB9*CA*y(J1 + 17);
					if (J + 18 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB10*CA*y(J1 + 18);
					if (J + 19 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB11*CA*y(J1 + 19);



					DY(J1) = AKP * (SOM3 + FR * (SOM4)
						+(SOM5)+FR * (SOM4)
						+(ABB1*SY1_1 + ABB2 * SY1_2 + ABB3 * SY1_3 + ABB4 * SY1_4
							+ ABB5 * SY1_5 + ABB6 * SY1_6 + ABB7 * SY1_7
							+ ABB8 * SY1_8 + ABB9 * SY1_9 + ABB10 * SY1_10
							+ ABB11 * SY1_11)
						+ FR * (ABB1*(SY1_12 + 0.5*SY12_1)
							+ ABB2 * (SY1_13 + 0.5*SY12_2)
							+ ABB3 * (SY1_14 + 0.5*SY12_3)
							+ ABB4 * (SY1_15 + 0.5*SY12_4)
							+ ABB5 * (SY1_16 + 0.5*SY12_5)
							+ ABB6 * (SY1_17 + 0.5*SY12_6)
							+ ABB7 * (SY1_18 + 0.5*SY12_7)
							+ ABB8 * (SY1_19 + 0.5*SY12_8)
							+ ABB9 * (SY1_20 + 0.5*SY12_9)
							+ ABB10 * (SY1_21 + 0.5*SY12_10)
							+ ABB11 * (SY1_22 + 0.5*SY12_11))
						+ 0.5*(ABB1*(SY1_23)
							+ABB2 * (SY1_24)
							+ABB3 * (SY1_25)
							+ABB4 * (SY1_26)
							+ABB5 * (SY1_27)
							+ABB6 * (SY1_28)
							+ABB7 * (SY1_29)
							+ABB8 * (SY1_30)
							+ABB9 * (SY1_31)
							+ABB10 * (SY1_32)
							+ABB11 * (SY1_33)));


					if (J + 3 >= LC) DY(J1) += AKP * FR*(CD*y(J1 + 3) + y(J1 + 3));
					if (J + 4 >= LC)  DY(J1) += AKP * CA*y(J2 + 4);
					if (J + 6 >= LC)  DY(J1) += AKP * FR*ABB1*CA*y(J1 + 6);
					if (J + 7 >= LC)  DY(J1) += AKP * FR*ABB2*CA*y(J1 + 7);
					if (J + 8 >= LC)  DY(J1) += AKP * FR*ABB3*CA*y(J1 + 8);
					if (J + 9 >= LC)  DY(J1) += AKP * FR*ABB4*CA*y(J1 + 9);
					if (J + 10 >= LC) DY(J1) += AKP * FR*ABB5*CA*y(J1 + 10);
					if (J + 11 >= LC) DY(J1) += AKP * FR*ABB6*CA*y(J1 + 11);
					if (J + 12 >= LC) DY(J1) += AKP * FR*ABB7*CA*y(J1 + 12);
					if (J + 13 >= LC) DY(J1) += AKP * FR*ABB8*CA*y(J1 + 13);
					if (J + 14 >= LC) DY(J1) += AKP * FR*ABB9*CA*y(J1 + 14);
					if (J + 15 >= LC) DY(J1) += AKP * FR*ABB10*CA*y(J1 + 15);
					if (J + 16 >= LC) DY(J1) += AKP * FR*ABB11*CA*y(J1 + 16);
					if (J + 9 >= LC)  DY(J1) += 0.5*AKP*ABB1*CA*y(J2 + 9);
					if (J + 10 >= LC) DY(J1) += 0.5*AKP*ABB2*CA*y(J2 + 10);
					if (J + 11 >= LC) DY(J1) += 0.5*AKP*ABB3*CA*y(J2 + 11);
					if (J + 12 >= LC) DY(J1) += 0.5*AKP*ABB4*CA*y(J2 + 12);
					if (J + 13 >= LC) DY(J1) += 0.5*AKP*ABB5*CA*y(J2 + 13);
					if (J + 14 >= LC) DY(J1) += 0.5*AKP*ABB6*CA*y(J2 + 14);
					if (J + 15 >= LC) DY(J1) += 0.5*AKP*ABB7*CA*y(J2 + 15);
					if (J + 16 >= LC) DY(J1) += 0.5*AKP*ABB8*CA*y(J2 + 16);
					if (J + 17 >= LC) DY(J1) += 0.5*AKP*ABB9*CA*y(J2 + 17);
					if (J + 18 >= LC) DY(J1) += 0.5*AKP*ABB10*CA*y(J2 + 18);
					if (J + 19 >= LC) DY(J1) += 0.5*AKP*ABB11*CA*y(J2 + 19);



					DY(J2) = AKP * (1. / 2.*SOM6*CA + SOM7 * CA
						+ (ABB1*(SY2_12)
							+ABB2 * (SY2_13)
							+ABB3 * (SY2_14)
							+ABB4 * (SY2_15)
							+ABB5 * (SY2_16)
							+ABB6 * (SY2_17)
							+ABB7 * (SY2_18)
							+ABB8 * (SY2_19)
							+ABB9 * (SY2_20)
							+ABB10 * (SY2_21)
							+ABB11 * (SY2_22))
						+ FR * (ABB1*SY2_1 + ABB2 * SY2_2 + ABB3 * SY2_3
							+ ABB4 * SY2_4 + ABB5 * SY2_5 + ABB6 * SY2_6
							+ ABB7 * SY2_7 + ABB8 * SY2_8 + ABB9 * SY2_9
							+ ABB10 * SY2_10 + ABB11 * SY2_11));


					if (J + 3 >= LC) DY(J2) += AKP * CA*CD*y(J2 + 3);
					if (J + 6 >= LC) DY(J2) += AKP * ABB1*CA*y(J2 + 6);
					if (J + 7 >= LC) DY(J2) += AKP * ABB2*CA*y(J2 + 7);
					if (J + 8 >= LC) DY(J2) += AKP * ABB3*CA*y(J2 + 8);
					if (J + 9 >= LC) DY(J2) += AKP * ABB4*CA*y(J2 + 9);
					if (J + 10 >= LC) DY(J2) += AKP * ABB5*CA*y(J2 + 10);
					if (J + 11 >= LC) DY(J2) += AKP * ABB6*CA*y(J2 + 11);
					if (J + 12 >= LC) DY(J2) += AKP * ABB7*CA*y(J2 + 12);
					if (J + 13 >= LC) DY(J2) += AKP * ABB8*CA*y(J2 + 13);
					if (J + 14 >= LC) DY(J2) += AKP * ABB9*CA*y(J2 + 14);
					if (J + 15 >= LC) DY(J2) += AKP * ABB10*CA*y(J2 + 15);
					if (J + 16 >= LC) DY(J2) += AKP * ABB11*CA*y(J2 + 16);

				}

				else if (J == 3)
				{
					DY(J - 1) = AKP * (SOM1 + (1. - FR)*(SOM2)
						+0.5*(ABB1*SY_1 + ABB3 * SY_3
							+ ABB4 * SY_4 + ABB5 * SY_5 + ABB6 * SY_6
							+ ABB7 * SY_7 + ABB8 * SY_8 + ABB9 * SY_9
							+ ABB10 * SY_10 + ABB11 * SY_11)
						+ 0.5*FR*(ABB1*(SY_12)
							+ABB3 * (SY_14)
							+ABB4 * (SY_15)
							+ABB5 * (SY_16)
							+ABB6 * (SY_17)
							+ABB7 * (SY_18)
							+ABB8 * (SY_19)
							+ABB9 * (SY_20)
							+ABB10 * (SY_21)
							+ABB11 * (SY_22))
						+ ABB2 * CORRP(3 + J - 1));

					if (J + 4 >= LC) DY(J - 1) += AKP * (1. - FR)*CA*y(J1 + 4);
					if (J + 9 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB1*CA*y(J1 + 9);
					if (J + 11 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB3*CA*y(J1 + 11);
					if (J + 12 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB4*CA*y(J1 + 12);
					if (J + 13 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB5*CA*y(J1 + 13);
					if (J + 14 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB6*CA*y(J1 + 14);
					if (J + 15 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB7*CA*y(J1 + 15);
					if (J + 16 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB8*CA*y(J1 + 16);
					if (J + 17 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB9*CA*y(J1 + 17);
					if (J + 18 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB10*CA*y(J1 + 18);
					if (J + 19 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB11*CA*y(J1 + 19);


					DY(J1) = AKP * (SOM3 + FR * (CD + 1.)*SOM4 +
						(SOM5)*CD
						+ (ABB1*SY1_1 + ABB2 * SY1_2 + ABB3 * SY1_3 + ABB4 * SY1_4
							+ ABB5 * SY1_5 + ABB6 * SY1_6 + ABB7 * SY1_7
							+ ABB8 * SY1_8 + ABB9 * SY1_9 + ABB10 * SY1_10
							+ ABB11 * SY1_11)
						+ FR * (ABB1*(SY1_12 + CD / (1. + CD)*SY12_1)
							+ ABB2 * (SY1_13 + CD / (1. + CD)*SY12_2)
							+ ABB3 * (SY1_14 + CD / (1. + CD)*SY12_3)
							+ ABB4 * (SY1_15 + CD / (1. + CD)*SY12_4)
							+ ABB5 * (SY1_16 + CD / (1. + CD)*SY12_5)
							+ ABB6 * (SY1_17 + CD / (1. + CD)*SY12_6)
							+ ABB7 * (SY1_18 + CD / (1. + CD)*SY12_7)
							+ ABB8 * (SY1_19 + CD / (1. + CD)*SY12_8)
							+ ABB9 * (SY1_20 + CD / (1. + CD)*SY12_9)
							+ ABB10 * (SY1_21 + CD / (1. + CD)*SY12_10)
							+ ABB11 * (SY1_22 + CD / (1. + CD)*SY12_11))
						+ CD / (1. + CD)*(ABB1*(SY1_23)
							+ABB2 * (SY1_24)
							+ABB3 * (SY1_25)
							+ABB4 * (SY1_26)
							+ABB5 * (SY1_27)
							+ABB6 * (SY1_28)
							+ABB7 * (SY1_29)
							+ABB8 * (SY1_30)
							+ABB9 * (SY1_31)
							+ABB10 * (SY1_32)
							+ABB11 * (SY1_33)));


					if (J + 3 >= LC) DY(J1) += AKP * (CD*y(J1 + 3));
					if (J + 4 >= LC) DY(J1) += AKP * CA*CD*y(J2 + 4);
					if (J + 6 >= LC) DY(J1) += AKP * FR*ABB1*CA*y(J1 + 6);
					if (J + 7 >= LC) DY(J1) += AKP * FR*ABB2*CA*y(J1 + 7);
					if (J + 8 >= LC) DY(J1) += AKP * FR*ABB3*CA*y(J1 + 8);
					if (J + 9 >= LC) DY(J1) += AKP * FR*ABB4*CA*y(J1 + 9);
					if (J + 10 >= LC) DY(J1) += AKP * FR*ABB5*CA*y(J1 + 10);
					if (J + 11 >= LC) DY(J1) += AKP * FR*ABB6*CA*y(J1 + 11);
					if (J + 12 >= LC) DY(J1) += AKP * FR*ABB7*CA*y(J1 + 12);
					if (J + 13 >= LC) DY(J1) += AKP * FR*ABB8*CA*y(J1 + 13);
					if (J + 14 >= LC) DY(J1) += AKP * FR*ABB9*CA*y(J1 + 14);
					if (J + 15 >= LC) DY(J1) += AKP * FR*ABB10*CA*y(J1 + 15);
					if (J + 16 >= LC) DY(J1) += AKP * FR*ABB11*CA*y(J1 + 16);
					if (J + 9 >= LC) DY(J1) += CD / (1. + CD)*AKP*ABB1*CA* y(J2 + 9);
					if (J + 10 >= LC) DY(J1) += CD / (1. + CD)*AKP*ABB2*CA* y(J2 + 10);
					if (J + 11 >= LC) DY(J1) += CD / (1. + CD)*AKP*ABB3*CA * y(J2 + 11);
					if (J + 12 >= LC) DY(J1) += CD / (1. + CD)*AKP*ABB4*CA* y(J2 + 12);
					if (J + 13 >= LC) DY(J1) += CD / (1. + CD)*AKP*ABB5*CA * y(J2 + 13);
					if (J + 14 >= LC) DY(J1) += CD / (1. + CD)*AKP*ABB6*CA* y(J2 + 14);
					if (J + 15 >= LC) DY(J1) += CD / (1. + CD)*AKP*ABB7*CA* y(J2 + 15);
					if (J + 16 >= LC) DY(J1) += CD / (1. + CD)*AKP*ABB8*CA* y(J2 + 16);
					if (J + 17 >= LC) DY(J1) += CD / (1. + CD)*AKP*ABB9*CA* y(J2 + 17);
					if (J + 18 >= LC) DY(J1) += CD / (1. + CD)*AKP*ABB10*CA* y(J2 + 18);
					if (J + 19 >= LC) DY(J1) += CD / (1. + CD)*AKP*ABB11*CA* y(J2 + 19);

				}

				else if (J == 2)
				{
					DY(J - 1) = AKP * (SOM1 + (1. - FR)*(SOM2)
							+0.5*(ABB2*SY_2 + ABB3 * SY_3
							+ ABB4 * SY_4 + ABB5 * SY_5 + ABB6 * SY_6
							+ ABB7 * SY_7 + ABB8 * SY_8 + ABB9 * SY_9
							+ ABB10 * SY_10 + ABB11 * SY_11)
							+ 0.5*FR*(ABB2*(SY_13)
							+ABB3 * (SY_14)
							+ABB4 * (SY_15)
							+ABB5 * (SY_16)
							+ABB6 * (SY_17)
							+ABB7 * (SY_18)
							+ABB8 * (SY_19)
							+ABB9 * (SY_20)
							+ABB10 * (SY_21)
							+ABB11 * (SY_22))
						+ ABB1 * CORRP(3 + J - 1));


					if (J + 4 >= LC) DY(J - 1) += AKP * (1. - FR)*CA*y(J1 + 4);
					if (J + 10 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB2*CA*y(J1 + 10);
					if (J + 11 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB3*CA*y(J1 + 11);
					if (J + 12 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB4*CA*y(J1 + 12);
					if (J + 13 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB5*CA*y(J1 + 13);
					if (J + 14 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB6*CA*y(J1 + 14);
					if (J + 15 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB7*CA*y(J1 + 15);
					if (J + 16 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB8*CA*y(J1 + 16);
					if (J + 17 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB9*CA*y(J1 + 17);
					if (J + 18 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB10*CA*y(J1 + 18);
					if (J + 19 >= LC) DY(J - 1) += 0.5*AKP*FR*ABB11*CA*y(J1 + 19);

				}


			} // end for gas

		}

		// Printo
		{
			// Gaseous phase mass
			{
				double SOMMAPL = 0.;
				double SOMMAOL = 0.;
				double SOMMADL = 0.;

				for (int I = 1; I <= LC - 1; I++)
				{
					SOMMAPL += y(I-1)*(I*14. + 2.);
					SOMMAOL += y(I-1 + NPA)*(I*14.);
					SOMMADL += y(I-1 + NPA + NOL)*(I*14. - 2.);
				}

				const double GAS = SOMMAPL + SOMMAOL + SOMMADL;
			}

			// Liquid phase mass
			{
				double RESPL = 0.;
				double RESOL = 0.;
				double RESDL = 0.;

				double PMOL = 0.;
				double DMOL = 0.;
				double AOMOL = 0.;

				for (int I = NPA; I >= LC; I--)
				{
					RESPL += y(I-1)*(I*14. + 2.);
					PMOL += y(I-1);

					RESOL += y(I-1 + NPA)*(I*14.);
					AOMOL += y(I-1 + NPA);

					RESDL += y(I-1 + NPA + NOL)*(I*14. - 2.);
					DMOL += y(I-1 + NPA + NOL);
				}
					
				const double RES = RESPL + RESOL + RESDL;

				//const double TG = RES / WG;
			}

			// Gaseous-phase distribution
			{
				Eigen::VectorXd Pmass(NPA);
				Eigen::VectorXd Omass(NPA);
				Eigen::VectorXd Dmass(NPA);
				for (int I = 1; I <= LC - 1; I++)
				{
					Pmass(I-1) = y(I-1)*(I*14. + 2);
					Omass(I-1) = y(I-1 + NPA)*(I*14.);
					Dmass(I-1) = y(I-1 + 2 * NPA)*(I*14. - 2);
				}
			}

			// Overall distribution
			{
				Eigen::VectorXd mass_distribution(ne/3);
				for (int i = 1; i <= ne/3; i++)
					mass_distribution(i - 1) = y(i-1)*(i*14. + 2.) + 
											   y(i-1 + NPA)*(i*14.) + 
											   y(i-1 + 2*NPA)*(i*14. - 2.);

				std::ofstream fOutput("DistributionPP.out", std::ios::out);
				fOutput.setf(std::ios::scientific);

				const double sum = mass_distribution.sum();
				std::cout << sum << std::endl;
				for (int i = 1; i < ne/3; i++)
					fOutput << i << " " << mass_distribution(i-1) << " " << mass_distribution(i-1) / sum << std::endl;

				fOutput.close();
			}
		}

		std::ofstream fLiq("Liq.out", std::ios::out);
		fLiq.setf(std::ios::scientific);
		for (int J = NPA; J >= LC; J--)
		{
			fLiq << std::setw(12) << std::right << J;
			fLiq << std::setw(24) << std::right << std::setprecision(15) << y(J - 1);
			fLiq << std::setw(24) << std::right << std::setprecision(15) << DY(J - 1);
			fLiq << std::endl;
		}
		fLiq.close();

		std::ofstream fGas("Gas.out", std::ios::out);
		fGas.setf(std::ios::scientific);
		for (int J = LC - 1; J >= 2; J--)
		{
			fGas << std::setw(12) << std::right << J;
			fGas << std::setw(24) << std::right << std::setprecision(15) << y(J - 1);
			fGas << std::setw(24) << std::right << std::setprecision(15) << DY(J - 1);
			fGas << std::endl;
		}
		fGas.close();
		
	}

	}



	getchar();
	return 1;
}


void InitialDistribution(	Eigen::VectorXd& y, const double epsi, const double MWm, const double MWp, 
							int& N, double& wg)
{
	y.setZero();

	const int ls = static_cast<int>( (y.size() - 1) / 3.);

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
			FatalErrorMessage("Max dimension error");

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

	std::ofstream fDistribution("DistributionPP.out", std::ios::out);
	fDistribution.setf(std::ios::scientific);

	std::ofstream fFractions("Fractions.out", std::ios::out);
	fFractions.setf(std::ios::scientific);

	double t = 0.;
	for (int i = 1; i <= N; i++)
	{
		const double s = y(i-1)*(i*MWm + 2.);
		t += s;

		fDistribution << i * MWm << " " << s << " " << t << std::endl;
		fFractions << i * MWm << " " << y(i) << std::endl;
	}

	wg = t;	// TODO

	fDistribution.close();
	fFractions.close();

	int NPA = N;
	int NOL = N;
	int NDO = N;

	y(NPA + NOL - 1) = 0.;
	y(ne - 1 - 1) = 0.;
	y(ne - 2 - 1) = 0.;
}

void Corrections(const Eigen::VectorXd& Y, Eigen::VectorXd& CORRO, Eigen::VectorXd& CORRP,
				 const int NPA, const int LC, const double CA, const double CD)
{
	double SOMXX = Y(NPA-1);
	double SOMZZ1 = Y(2 * NPA - 1);
	double SOMZZ2 = Y(2 * NPA - 1);
	double SOMDD = Y(3 * NPA - 1);
	double ASOMXX = SOMXX;
	double ASOMZZ1 = SOMZZ1;
	double ASOMZZ2 = SOMZZ2;
	double ASOMDD = SOMDD;

	int IVAL = 15;

	for (int jj = NPA - 4; jj >= IVAL + 2; jj--)
	{
		if (jj + 3 >= LC)
			ASOMXX += Y(jj + 3 - 1);
		SOMXX += ASOMXX;
	}

	for (int jj = NPA - 6; jj >= IVAL + 2; jj--)
	{
		if (jj + 5 >= LC)
			ASOMZZ1 += Y(jj + 5 + NPA - 1);
		SOMZZ1 += ASOMZZ1;
	}

	for (int jj = NPA - 4; jj >= IVAL + 4; jj--)
	{
		if (jj + 3 >= LC)
			ASOMZZ2 += Y(jj + 3 + NPA);
		SOMZZ2 += ASOMZZ2;
	}
		
	for (int jj = NPA - 6; jj >= IVAL + 4; jj--)
	{
		if (jj + 5 >= LC)
			ASOMDD += Y(2 * NPA + jj + 5 - 1);
		SOMDD += ASOMDD;
	}


	double SOMP1 = 0.0;
	for (int j = 2 * IVAL; j <= NPA; j++)
	{
		if (j >= LC)
			SOMP1 += Y(j - 1);
	}

	double SOMD1 = 0.0;
	for (int j = IVAL+1; j <= NPA; j++)
	{
		if (j >= LC)
			SOMD1 += Y(j + 2*NPA-1);
	}
		
	double SOMD2 = 0.0;
	for (int j = 2*IVAL + 5; j <= NPA; j++)
	{
		if (j >= LC)
			SOMD2 += Y(j + 2*NPA-1);
	}

	double SOMO1 = 0.0;
	for (int j = IVAL + 6; j <= NPA; j++)
	{
		if (j >= LC)
			SOMO1 += Y(j + NPA - 1);
	}

	double SOMO2 = 0.0;
	for (int j = 2 * IVAL + 3; j <= NPA; j++)
	{
		if (j >= LC)
			SOMO2 += Y(j + NPA - 1);
	}

	double SOMP2 = 0.0;
	for (int j = IVAL + 4; j <= NPA; j++)
	{
		if (j >= LC)
			SOMP2 += Y(j - 1);
	}

	CORRO(IVAL - 1) = 0.5*(SOMXX + SOMP1) + 0.5*(SOMDD + CA * SOMD1) +
		CD / (1. + CD)*SOMD1 + 0.5*SOMD2 +
		0.5*(0.5*(SOMZZ1 + CA * SOMO1) + 0.5*SOMO2) +
		0.5*(0.5*SOMZZ2 + CD / (1. + CD)*SOMO1 + 0.5*SOMO2);


	if (IVAL + 7 >= LC)
		CORRO(IVAL - 1) += CD / (1. + CD)*CA*Y(2 * NPA + IVAL + 7 - 1);

	if (2 * IVAL + 4 >= LC)
		CORRO(IVAL - 1) += 0.5*CA*Y(2 * NPA + 2 * IVAL + 4 - 1);

	if (2 * IVAL + 1 >= LC)
		CORRO(IVAL - 1) += 0.25*CA*Y(2 * IVAL + 1 + NPA - 1);

	if (2 * IVAL + 2 >= LC)
		CORRO(IVAL - 1) += 0.25*Y(2 * IVAL + 2 + NPA - 1);

	CORRP(IVAL - 1) = 0.5*SOMXX + (SOMP2 + SOMP2 + 0.5*SOMP1) +
		0.5*(SOMDD + CA * SOMD1) + 1. / (1. + CD) * SOMD1 + SOMD1 + SOMD1 +
		0.5*(0.5*(SOMZZ1 + CA * SOMO1) + 0.5*SOMO2 + SOMO1 + SOMO1) +
		0.5*(0.5*SOMZZ2 + 1. / (1. + CD)*SOMO1 + SOMO1 + SOMO1);

	if (IVAL + 3 >= LC)
		CORRP(IVAL - 1) += Y(IVAL + 3 - 1);

	if (IVAL + 7 >= LC)
		CORRP(IVAL - 1) += CA / (1. + CD)*Y(IVAL + 7 + 2 * NPA - 1);

	if (IVAL + 6 >= LC)
		CORRP(IVAL - 1) += CA * Y(2 * NPA + IVAL + 6 - 1);

	if (IVAL + 7 >= LC)
		CORRP(IVAL - 1) += Y(2 * NPA + IVAL + 7 - 1);

	if (IVAL + 5 >= LC)
		CORRP(IVAL - 1) += CA * Y(2 * NPA + IVAL + 5 - 1);

	if (IVAL + 7 >= LC)
		CORRP(IVAL - 1) += Y(2 * NPA + IVAL + 7 - 1);

	if (IVAL + 6 >= LC)
		CORRP(IVAL - 1) += Y(2 * NPA + IVAL + 6 - 1);

	if (2 * IVAL + 1 >= LC)
		CORRP(IVAL - 1) += 0.25*CA*Y(NPA + 2 * IVAL + 1 - 1);

	if (2 * IVAL + 2 >= LC)
		CORRP(IVAL - 1) += 0.25*Y(NPA + 2 * IVAL + 2 - 1);

	if (IVAL + 5 >= LC)
		CORRP(IVAL - 1) += 0.5*CA*Y(NPA + IVAL + 5 - 1);

	if (IVAL + 4 >= LC)
		CORRP(IVAL - 1) += 0.5*CA*Y(NPA + IVAL + 4 - 1);

	if (IVAL + 5 >= LC)
		CORRP(IVAL - 1) += 0.5*Y(NPA + IVAL + 5 - 1);

	if (IVAL + 5 >= LC)
		CORRP(IVAL - 1) += 0.5*Y(NPA + IVAL + 5 - 1);

	if (IVAL + 4 >= LC)
		CORRP(IVAL - 1) += 0.5*Y(NPA + IVAL + 4 - 1);

	if (IVAL + 5 >= LC)
		CORRP(IVAL - 1) += 0.5*Y(NPA + IVAL + 5 - 1);

	for (int IVAL = 14; IVAL >= 5; IVAL--)
	{
		if (2 * IVAL + 1 >= LC)
			SOMP1 += Y(2 * IVAL + 1 - 1);

		if (2 * IVAL >= LC)
			SOMP1 += Y(2 * IVAL - 1);

		if (IVAL + 4 >= LC)
			SOMP2 += Y(IVAL + 4 - 1);

		if (IVAL + 8 >= LC)
			SOMD1 += Y(2 * NPA + IVAL + 8 - 1);

		if (2 * IVAL + 6 >= LC)
			SOMD2 += Y(2 * NPA + 2 * IVAL + 5 + 1 - 1);

		if (2 * IVAL + 5 >= LC)
			SOMD2 += Y(2 * NPA + 2 * IVAL + 5 - 1);

		if (IVAL + 6 >= LC)
			SOMO1 += Y(NPA + IVAL + 6 - 1);

		if (2 * IVAL + 4 >= LC)
			SOMO2 += Y(NPA + 2 * IVAL + 3 + 1 - 1);

		if (2 * IVAL + 3 >= LC)
			SOMO2 += Y(NPA + 2 * IVAL + 3 - 1);

		int JJ = IVAL + 2;

		if (JJ + 3 >= LC)
			ASOMXX += Y(JJ + 3 - 1);
		SOMXX += ASOMXX;

		if (JJ + 5 >= LC)
			ASOMZZ1 += Y(JJ + NPA + 5 - 1);
		SOMZZ1 += ASOMZZ1;

		JJ = IVAL + 4;

		if (JJ + 3 >= LC)
			ASOMZZ2 += Y(JJ + NPA + 3 - 1);
		SOMZZ2 += ASOMZZ2;

		if (JJ + 5 >= LC)
			ASOMDD += Y(JJ + 2 * NPA + 5 - 1);
		SOMDD += ASOMDD;

		CORRO(IVAL - 1) = 0.5*(SOMXX + SOMP1) + 0.5*(SOMDD + CA * SOMD1) +
			CD / (1. + CD)*SOMD1 + 0.5*SOMD2 + 0.5*(0.5*(SOMZZ1 + CA * SOMO1) +
				0.5*SOMO2) + 0.5*(0.5*SOMZZ2 + CD / (1. + CD)*SOMO1 + 0.5*SOMO2);

		if (IVAL + 7 >= LC)
			CORRO(IVAL - 1) += CD / (1. + CD)*CA*Y(2 * NPA + IVAL + 7 - 1);

		if (2 * IVAL + 4 >= LC)
			CORRO(IVAL - 1) += 0.5*CA*Y(2 * NPA + 2 * IVAL + 4 - 1);

		if (2 * IVAL + 1 >= LC)
			CORRO(IVAL - 1) += 0.25*CA*Y(2 * IVAL + 1 + NPA - 1);


		if (2 * IVAL + 2 >= LC)
			CORRO(IVAL - 1) += 0.25*Y(2 * IVAL + 2 + NPA - 1);


		CORRP(IVAL - 1) = 0.5*SOMXX + (SOMP2 + SOMP2 + 0.5*SOMP1) +
			0.5*(SOMDD + CA * SOMD1) + 1. / (1. + CD) *
			SOMD1 + SOMD1 + SOMD1 +
			0.5*(0.5*(SOMZZ1 + CA * SOMO1) + 0.5*SOMO2 + SOMO1 + SOMO1) +
			0.5*(0.5*SOMZZ2 + 1. / (1. + CD)*SOMO1 + SOMO1 + SOMO1);


		if (IVAL + 3 >= LC)
			CORRP(IVAL - 1) += Y(IVAL + 3 - 1);

		if (IVAL + 7 >= LC)
			CORRP(IVAL - 1) += CA / (1. + CD)*Y(IVAL + 7 + 2 * NPA - 1);

		if (IVAL + 6 >= LC)
			CORRP(IVAL - 1) += CA * Y(2 * NPA + IVAL + 6 - 1);

		if (IVAL + 7 >= LC)
			CORRP(IVAL - 1) += Y(2 * NPA + IVAL + 7 - 1);

		if (IVAL + 5 >= LC)
			CORRP(IVAL - 1) += CA * Y(2 * NPA + IVAL + 5 - 1);

		if (IVAL + 7 >= LC)
			CORRP(IVAL - 1) += Y(2 * NPA + IVAL + 7 - 1);

		if (IVAL + 6 >= LC)
			CORRP(IVAL - 1) += Y(2 * NPA + IVAL + 6 - 1);

		if (2 * IVAL + 1 >= LC)
			CORRP(IVAL - 1) += 0.25*CA*Y(NPA + 2 * IVAL + 1 - 1);

		if (2 * IVAL + 2 >= LC)
			CORRP(IVAL - 1) += 0.25*Y(NPA + 2 * IVAL + 2 - 1);

		if (IVAL + 5 >= LC)
			CORRP(IVAL - 1) += 0.5*CA*Y(NPA + IVAL + 5 - 1);

		if (IVAL + 4 >= LC)
			CORRP(IVAL - 1) += 0.5*CA*Y(NPA + IVAL + 4 - 1);

		if (IVAL + 5 >= LC)
			CORRP(IVAL - 1) += 0.5*Y(NPA + IVAL + 5 - 1);

		if (IVAL + 5 >= LC)
			CORRP(IVAL - 1) += 0.5*Y(NPA + IVAL + 5 - 1);

		if (IVAL + 4 >= LC)
			CORRP(IVAL - 1) += 0.5*Y(NPA + IVAL + 4 - 1);

		if (IVAL + 5 >= LC)
			CORRP(IVAL - 1) += 0.5*Y(NPA + IVAL + 5 - 1);

	}
}

void Aggregation(const Eigen::VectorXd& Y, const int J, 
	double& SOM1, double& SOM2, double& SOM3, double& SOM4, double& SOM5, double& SOM6, double& SOM7,
	Eigen::VectorXd& SY, Eigen::VectorXd& SY1, Eigen::VectorXd& SY12, Eigen::VectorXd& SY2, 
	const int NPA, const int LC)
{

	if ((NPA - 2 >= J) && (J + 2 >= LC))
	{
		const int M0 = J + 2;
		const int M1 = M0 + NPA;
		//M2 = M0 + NPA * 2;
		SOM3 += Y(M0-1);
		SOM6 += Y(M1-1);
	}

	if ((NPA - 3 >= J) && (J + 3 >= LC))
	{
		const int M0 = J + 3;
		//M1 = M0 + NPA;
		//M2 = M0 + NPA * 2;
		SOM1 += Y(M0 - 1);
	}

	if ((NPA - 4 >= J) && (J + 4 >= LC))
	{
		const int M0 = J + 4;
		const int M1 = M0 + NPA;
		const int M2 = M0 + NPA * 2;
		SOM4 += Y(M1 - 1);
		SOM7 += Y(M2 - 1);
	}

	if ((NPA - 5 >= J) && (J + 5 >= LC))
	{
		const int M0 = J + 5;
		const int M1 = M0 + NPA;
		const int M2 = M0 + NPA * 2;
		SY1(1 - 1) = SY1(1) + Y(M0 - 1);
		SY2(1 - 1) = SY2(1) + Y(M1 - 1);
		SOM2 += Y(M1 - 1);
		SOM5 += Y(M2 - 1);
	}

	if ((NPA - 6 >= J) && (J + 6 >= LC))
	{
		const int M0 = J + 6;
		const int M1 = M0 + NPA;
		SY2(2 - 1) += Y(M1 - 1);
		SY1(2 - 1) += Y(M0 - 1);
	}

	if ((NPA - 7 >= J) && (J + 7 >= LC))
	{
		const int M0 = J + 7;
		const int M1 = M0 + NPA;
		const int M2 = M0 + NPA * 2;
		SY1(3 - 1) += Y(M0 - 1);
		SY1(12 - 1) += Y(M1 - 1);
		SY2(3 - 1) += Y(M1 - 1);
		SY2(12 - 1) += Y(M2 - 1);
	}

	if ((NPA - 8 >= J) && (J + 8 >= LC))
	{
		const int M0 = J + 8;
		const int M1 = M0 + NPA;
		const int M2 = M0 + NPA * 2;
		SY(1 - 1) += Y(M0 - 1);
		SY1(4 - 1) += Y(M0 - 1);
		SY2(4 - 1) += Y(M1 - 1);
		SY12(1 - 1) += Y(M1 - 1);
		SY1(13 - 1) += Y(M1 - 1);
		SY2(13 - 1) += Y(M2 - 1);
	}

	if ((NPA - 9 >= J) && (J + 9 >= LC))
	{
		const int M0 = J + 9;
		const int M1 = M0 + NPA;
		const int M2 = M0 + NPA * 2;
		SY2(14 - 1) += Y(M2 - 1);
		SY1(14 - 1) += Y(M1 - 1);
		SY12(2 - 1) += Y(M1 - 1);
		SY(2 - 1) += Y(M0 - 1);
		SY1(5 - 1) += Y(M0 - 1);
		SY2(5 - 1) += Y(M1 - 1);
	}

	if ((NPA - 10 >= J) && (J + 10 >= LC))
	{
		const int M0 = J + 10;
		const int M1 = M0 + NPA;
		const int M2 = M0 + NPA * 2;
		SY(3 - 1) += Y(M0 - 1);
		SY(12 - 1) += Y(M1 - 1);
		SY1(15 - 1) += Y(M1 - 1);
		SY12(3 - 1) += Y(M1 - 1);
		SY1(23 - 1) += Y(M2 - 1);
		SY2(15 - 1) += Y(M2 - 1);
		SY1(6 - 1) += Y(M0 - 1);
		SY2(6 - 1) += Y(M1 - 1);
	}

	if ((NPA - 11 >= J) && (J + 11 >= LC))
	{
		const int M0 = J + 11;
		const int M1 = M0 + NPA;
		const int M2 = M0 + NPA * 2;
		SY(4 - 1) += Y(M0 - 1);
		SY12(4 - 1) += Y(M1 - 1);
		SY(13 - 1) += Y(M1 - 1);
		SY1(24 - 1) += Y(M2 - 1);
		SY2(16 - 1) += Y(M2 - 1);
		SY1(16 - 1) += Y(M1 - 1);
		SY1(7 - 1) += Y(M0 - 1);
		SY2(7 - 1) += Y(M1 - 1);
	}

	if ((NPA - 12 >= J) && (J + 12 >= LC))
	{
		const int M0 = J + 12;
		const int M1 = M0 + NPA;
		const int M2 = M0 + NPA * 2;
		SY(5 - 1) += Y(M0 - 1);
		SY12(5 - 1) += Y(M1 - 1);
		SY(14 - 1) += Y(M1 - 1);
		SY1(25 - 1) += Y(M2 - 1);
		SY2(17 - 1) += Y(M2 - 1);
		SY1(17 - 1) += Y(M1 - 1);
		SY1(8 - 1) += Y(M0 - 1);
		SY2(8 - 1) += Y(M1 - 1);
	}

	if ((NPA - 13 >= J) && (J + 13 >= LC))
	{
		const int M0 = J + 13;
		const int M1 = M0 + NPA;
		const int M2 = M0 + NPA * 2;
		SY(6 - 1) += Y(M0 - 1);
		SY12(6 - 1) += Y(M1 - 1);
		SY(15 - 1) += Y(M1 - 1);
		SY1(26 - 1) += Y(M2 - 1);
		SY2(18 - 1) += Y(M2 - 1);
		SY1(18 - 1) += Y(M1 - 1);
		SY1(9 - 1) += Y(M0 - 1);
		SY2(9 - 1) += Y(M1 - 1);
	}

	if ((NPA - 14 >= J) && (J + 14 >= LC))
	{
		const int M0 = J + 14;
		const int M1 = M0 + NPA;
		const int M2 = M0 + NPA * 2;
		SY1(27 - 1) += Y(M2 - 1);
		SY(16 - 1) += Y(M1 - 1);
		SY(7 - 1) += Y(M0 - 1);
		SY12(7 - 1) += Y(M1 - 1);
		SY2(19 - 1) += Y(M2 - 1);
		SY1(19 - 1) += Y(M1 - 1);
		SY1(10 - 1) += Y(M0 - 1);
		SY2(10 - 1) += Y(M1 - 1);
	}

	if ((NPA - 15 >= J) && (J + 15 >= LC))
	{
		const int M0 = J + 15;
		const int M1 = M0 + NPA;
		const int M2 = M0 + NPA * 2;
		SY(17 - 1) += Y(M1 - 1);
		SY1(28 - 1) += Y(M2 - 1);
		SY(8 - 1) += Y(M0 - 1);
		SY12(8 - 1) += Y(M1 - 1);
		SY2(20 - 1) += Y(M2 - 1);
		SY1(20 - 1) += Y(M1 - 1);
		SY1(11 - 1) += Y(M0 - 1);
		SY2(11 - 1) += Y(M1 - 1);
	}

	if ((NPA - 16 >= J) && (J + 16 >= LC))
	{
		const int M0 = J + 16;
		const int M1 = M0 + NPA;
		const int M2 = M0 + NPA * 2;
		SY(18 - 1) += Y(M1 - 1);
		SY1(29 - 1) += Y(M2 - 1);
		SY(9 - 1) += Y(M0 - 1);
		SY12(9 - 1) += Y(M1 - 1);
		SY1(21 - 1) += Y(M1 - 1);
		SY2(21 - 1) += Y(M2 - 1);
	}

	if ((NPA - 17 >= J) && (J + 17 >= LC))
	{
		const int M0 = J + 17;
		const int M1 = M0 + NPA;
		const int M2 = M0 + NPA * 2;
		SY(19 - 1) += Y(M1 - 1);
		SY1(30 - 1) += Y(M2 - 1);
		SY(10 - 1) += Y(M0 - 1);
		SY12(10 - 1) += Y(M1 - 1);
		SY2(22 - 1) += Y(M2 - 1);
		SY1(22 - 1) += Y(M1 - 1);
	}

	if ((NPA - 18 >= J) && (J + 18 >= LC))
	{
		const int M0 = J + 18;
		const int M1 = M0 + NPA;
		const int M2 = M0 + NPA * 2;
		SY(20 - 1) += Y(M1 - 1);
		SY1(31 - 1) += Y(M2 - 1);
		SY(11 - 1) += Y(M0 - 1);
		SY12(11 - 1) += Y(M1 - 1);
	}

	if ((NPA - 19 >= J) && (J + 19 >= LC))
	{
		const int M0 = J + 19;
		const int M1 = M0 + NPA;
		const int M2 = M0 + NPA * 2;
		SY(21 - 1) += Y(M1 - 1);
		SY1(32 - 1) += Y(M2 - 1);
	}

	if ((NPA - 20 >= J) && (J + 20 >= LC))
	{
		const int M0 = J + 20;
		const int M1 = M0 + NPA;
		const int M2 = M0 + NPA * 2;
		SY(22 - 1) += Y(M1 - 1);
		SY1(33 - 1) += Y(M2 - 1);
	}
}