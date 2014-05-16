#include <iostream>
#include <iomanip>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <complex>
#include <fftw3.h>
#include <vector>
#include <memory>
#include <string>
#include <gsl/gsl_sf_coupling.h>

#include "mkl.h"

#include <boost/algorithm/string.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/lexical_cast.hpp>

#include "array_structs.h"
#include "numerics.h"
#include "quantum_dynamics.h"

using namespace std;

typedef std::complex<double> cplx;
typedef std::vector< std::vector<int> > basis_states;

double gaussian(double x, double y, double xcen, double ycen);
void print_matrix(char* desc, int m, int n, double *a, int lda);
void check(int n);

int main (int argc, char *argv[])
{
	int ii, jj, kk;

/********Atomic Units Conversions***************/
	double HBAR = 1.0;
	double EMASS = 1.0;
	double ECHARGE = 1.0;
	double VACPERM = (1.0/(4.0*M_PI));
	double LEN = 0.0529177; // nm
	double VEL = 2.18e8; // cm/s
	double EN = 27.21; // eV
	double TIME = 2.42e-17; // s
	double FREQ = 4.13e16; // Hz
	double ELECFIELD = 5.14e9; // V/cm^2
	double LASERINTEN = 3.51e16; // W/cm^2 (0.5*vacuum_permitivity*speed_of_light*EField^2)
	double C = 2.998e10/VEL; //speed of light

/***************************
Get Data from the inputfile
***************************/

//	First copy input file and strip comments
	fstream inputfile;
	if (argc < 2)
	{
		cout << "Please provide the name of the json inputfile" << endl;
		exit(1);
	}
	else
	{
		inputfile.open(argv[1]);

	}

	ofstream inputcopy;
	inputcopy.open("output_data/inputs.json");
	string line;

	int found, found2;
	while (getline(inputfile,line))
	{
		found = line.find('/');
		found2 = line.find('/', found+1);
		if (found != line.npos && found2 == found+1)
		{
			inputcopy << line.erase(found, line.length()) << endl;
		}
		else
		{
			inputcopy << line << endl;
		}
	}
	inputcopy.close();

	boost::property_tree::ptree pt;
	boost::property_tree::json_parser::read_json("output_data/inputs.json",pt);

//	[FDTDField]
    string FieldFile 	= pt.get<string>("FDTDField.File","XAgFilmX-3col.txt");
    cout << "Reading file " << FieldFile << endl;
    double laser_power 	= pt.get<double>("FDTDField.laser_power",10.0);
    double laser_focus 	= pt.get<double>("FDTDField.laser_focus",5.0);
    double Field_a = pt.get<double>("FDTDField.meep_a",100.0);
    double Field_res = pt.get<double>("MolFile.meep_res",20.0);

//	[MoleculeParameters]

    string MolFile 	= pt.get<string>("MoleculeParameters.ParamFile","Molecules.json");
    string Molecule = pt.get<string>("MoleculeParameters.Molecule","N2");

//	[CalculationParameters]
	int Procs 	= pt.get<int>("CalculationParameters.Procs",1);
	int JMax 	= pt.get<int>("CalculationParameters.JStates",10);
	int UseM	= pt.get<int>("CalculationParameters.UseM",0);
	int UseOdd	= pt.get<int>("CalculationParameters.UseOdd",0);
	int AllField = pt.get<int>("CalculationParameters.AllField",0);
	int num_procs = pt.get<int>("CalculationParameters.Threads",1);

//	[TrajectoryParameters]
	double InitialX = pt.get<double>("TrajectoryParameters.InitialX",3.0);
    double InitialY = pt.get<double>("TrajectoryParameters.InitialY",25.0);
    double ttotal 	= pt.get<double>("TrajectoryParameters.ttotal",2.0e8);
    int tsteps = pt.get<int>("TrajectoryParameters.tsteps",2000);
    int dt_image_out = pt.get<int>("TrajectoryParameters",100);

//	[Outputs]
	int JBasisOut = pt.get<int>("Outputs.JBasis",0);
	int EigenVectorsOut = pt.get<int>("Outputs.EigenVectors",0);

/***************************
Get molecule data
***************************/

	boost::property_tree::ptree mol;
	boost::property_tree::json_parser::read_json(MolFile,mol);

	double amass1	= mol.get<double>(Molecule + ".Mass1",14.0);
	double amass2	= mol.get<double>(Molecule + ".Mass2",14.0);

	double bond 	= mol.get<double>(Molecule + ".BondLength",1.0);

	boost::property_tree::ptree pol;
	pol = mol.get_child(Molecule + ".Polarizability");
	double apara 	= pol.get<double>("ZZ",12.49002644);
	double aperp 	= pol.get<double>("XX",7.06531958);

/*******************
Set up file buffers
*******************/

	streambuf *psbuf, *backup;
	ofstream JBasisBuffer;

/**********************************
Get FDTD Array Dimensions
Allocate Space for 2D Array
Fill Array with data from file
***********************************/

	int xsize = 0;
	int ysize = 0;

	ifstream data(FieldFile.c_str());

	time_t t = time(0);   // get time now
    struct tm * now = localtime( & t );

	cout << "\n*************************************************************" << endl;
	cout << "Beginning calculation of molecular alignment in an FDTD field" << endl;
	cout << "Time: " << asctime(now);
	cout << "*************************************************************" << endl;

	while(getline(data,line))
	{
		vector <string> tokens;
		boost::split(tokens,line,boost::is_any_of( " ,\t\n" ));
		if (tokens.size() == 3)
		{
			int xtest = boost::lexical_cast<int>(tokens[0]);
			int ytest = boost::lexical_cast<int>(tokens[1]);
			if (xtest > xsize) xsize = xtest;
			if (ytest > ysize) ysize = ytest;
		}
	}
	xsize += 1;
	ysize += 1;
	cout << "The field file contains a " << xsize << " by " << ysize << " data file." << endl;
	cout << "(" << xsize*ysize << " grid points)" << endl;
	cout << "Loading EMF file... ";
	Array_2D <double> EMField(xsize,ysize);
	data.clear(); //Needed to clear the EOF indicator, else seekg function fails
	data.seekg(0,data.beg);
	double efield;

//	Convert Field to atomic units
	laser_focus /= 1.0e4; //Convert to cm
	double power_factor = sqrt(sqrt(10.0));
	double laser_intensity;
	laser_intensity = laser_power / (M_PI * pow(laser_focus,2)); // W/cm^2
	laser_intensity /= LASERINTEN; // to au

//	double field_E0 = sqrt(2.0*laser_intensity/(C*VACPERM)); // E0 value in atomic units

	while(getline(data,line))
	{
		vector <string> tokens;
		boost::split(tokens,line,boost::is_any_of( " ,\t\n" ));
		if (tokens.size() == 3)
		{
			ii = boost::lexical_cast<int>(tokens[0]);
			jj = boost::lexical_cast<int>(tokens[1]);
			efield = boost::lexical_cast<double>(tokens[2]);
//	File holds |E|^2, just want |E|
			EMField.set_elem(ii,jj,efield*laser_intensity);
		}
	}
	cout << "Done." << endl;
//	cout << "The laser intensity is " << scientific << setprecision(2) << laser_intensity*LASERINTEN << "W/cm^2"<< endl;

/**************************************
Define field iterator array for loop
***************************************/

	Array_2D <double> *EMF;
	Array_2D <double> Powers(13*4,1);
	if (AllField != 0)
	{
		EMF = &EMField;
	}
	else
	{
		ii = 0;
		for (laser_power = 1.0e3; laser_power < 1.0e16; laser_power *= power_factor,ii++)
		{
			//cout << laser_power << endl;
			Powers.grid[ii] = laser_power/LASERINTEN;
		}
		EMF = &Powers;
	}

/**************************************
Define the basis set using |JMK> states
***************************************/

	int NEq;
	int factor; //Used to specify only even states
	if (UseM == 0 && UseOdd == 0)
	{
		factor = 2;
	}
	else
	{
		factor = 1;
	}

	if (UseM == 0)
	{
		NEq = JMax+1;
	}
	else
	{
		NEq = (JMax+1)*(JMax+1);
	}
	basis_states RB(NEq, std::vector<int>(3));  //Element of vector represent J,K,M (RB = "RotationalBasis")
	if (UseM == 0)
	{
		kk = 0;
		for (ii = 0; ii < NEq; ii++)
		{
			RB[kk][0] = factor*ii;
			RB[kk][1] = 0; //Linear Molecules, K = 0
			RB[kk][2] = 0;
			kk++;
		}
	}
	else
	{
		kk = 0;
		for (ii = 0; ii <= JMax; ii++)
		{
			for (jj = -1*ii; jj <= ii; jj++)
			{
				RB[kk][0] = ii;
				RB[kk][1] = 0; //Linear Molecules, K = 0
				RB[kk][2] = jj;
				kk++;
			}
		}
	}



	if (JBasisOut == 1  || JBasisOut == 2)
	{
		if (JBasisOut == 2)
		{
			JBasisBuffer.open("output_data/JBasis.txt");
			backup = cout.rdbuf(); 			// backup cout's stream buffer
			psbuf = JBasisBuffer.rdbuf(); 	// get file's stream buffer
			cout.rdbuf(psbuf); 				//assign to cout
		}
//		Print to the desired output
		for (ii = 0; ii < RB.size(); ii++)
		{
			for (std::vector<int>::const_iterator j = RB[ii].begin(); j != RB[ii].end(); ++j)
			{
				cout << *j << ' ';
			}
			cout << endl;
		}
		if (JBasisOut == 2)
		{
			cout.rdbuf(backup);
			JBasisBuffer.close();
		}
	}


/*********************************
Define the Cos^2 Matrix
**********************************/

	Array_2D <double> cossq(NEq,NEq);

	for (ii = 0; ii < cossq.Nx; ii++)
	{
		for (jj = 0; jj < cossq.Ny; jj++)
		{
			cossq.set_elem(ii,jj,(1.0/3.0)*kron_delta(RB[ii].at(0),RB[jj].at(0))*kron_delta(RB[ii].at(2),RB[jj].at(2))
				+(2.0/3.0)*FMIME(RB[ii].at(0),RB[ii].at(1),RB[ii].at(2),0,0,RB[jj].at(0),RB[jj].at(1),RB[jj].at(2)));
		}
	}

	int n = NEq;
	int lda = NEq;
	int info, lwork;
	double wkopt;
	double *work;

// local arrays
	double *EigenValues, *EigenVectors;
	double E_JK, coeff, coupling;

/*********************************
Define the Hamiltonian
**********************************/

	amass1 *= 1823.0; 	//convert amu to au
	amass2 *= 1823.0;
	bond /= 0.528; 		//Angstroms to au
//	Define rotational constant
	double Be = pow(HBAR,2)*0.5*(1.0/pow(bond,2))*(amass1+amass2)/(amass2*amass1);
	double percent_complete,prev_kk;

	Array_2D <double> H_total_FFB(NEq,NEq);
	Array_2D <double> COSARRAY(EMF->Nx,EMF->Ny);
	Array_2D <double> ENERGYARRAY(EMF->Nx,EMF->Ny);

	ofstream EigenOut;
	if (EigenVectorsOut != 0)
	{
		EigenOut.open("output_data/EigenVectorPopulations.txt");
	}
//
	cout << "Beginning Calculation..." << endl << fixed;
	prev_kk = 0;
	for (kk = 0; kk < EMF->Nx*EMF->Ny; kk++)
	{
		efield = EMF->grid[kk];
		percent_complete = 100 * (float)kk/(float)(xsize*ysize);
		if ((int)percent_complete != prev_kk) cout << percent_complete << "\r";
		prev_kk = percent_complete;

		for (ii = 0; ii < H_total_FFB.Nx; ii++)
		{
			for (jj = 0; jj < H_total_FFB.Ny; jj++)
			{
				E_JK =  kron_delta(RB[ii].at(2),RB[jj].at(2))*kron_delta(RB[ii].at(0),RB[jj].at(0))*(Be*(RB[ii].at(0))*(RB[ii].at(0)+1) - 0.25*efield*aperp); //Be(J(J+1) - (1/4)|E|^2*a_perp
				coeff = (-1.0/6.0)*abs(apara-aperp)*efield; //EField already squared here
				coupling = FMIME(RB[ii].at(0),RB[ii].at(1),RB[ii].at(2),0,0,RB[jj].at(0),RB[jj].at(1),RB[jj].at(2));
				H_total_FFB.set_elem(ii,jj,kron_delta(RB[ii].at(1),RB[jj].at(1))*(E_JK + coeff*coupling));
			}
			if (ii < 9 && jj < 9) cout << endl;
		}
		cout.precision(1);

	/*************************************
	Find Eigenvectors of Total Hamiltonian
	**************************************/

		EigenValues = new double[n];
		EigenVectors = new double[n*n];
		for (ii = 0; ii < n*n; ii++)
		{
			EigenVectors[ii] = H_total_FFB.grid[ii];
		}

		lwork = -1;
	// get eigenvalues and eigenvectors
		dsyev(	"V",			// 'N' for eigenvalues, 'V' for eigenvalues and eigenvectors
				"U", 			// 'U'pper or 'L'ower part of matrix used
				&n, 			// matrix order
				EigenVectors, 	// matrix
				&lda,			// lead dimension of matrix
				EigenValues, 	// array for eigenvalues
				&wkopt,			//
				&lwork,			//
				&info); 		//
	// Note: dsyev gives needed workspace as output for lwork = -1
	// Modify and run again for actual eigenvector calculation
		lwork = (int)wkopt;
		work = new double[lwork];
		dsyev("V","U",&n,EigenVectors,&lda,EigenValues,work,&lwork,&info);
		if (info > 0)
		{
			cout << "DSYEV failed to find the eigenvalues." << endl;
			exit(0);
		}
	/****************************
	Get cos^2 of the ground state
	*****************************/

	//	Place Eigenvector ground state in it's own array
		Array_1D <double> GroundState(n);
		for (ii=0; ii<n; ii++)
		{
			GroundState.set_elem(ii,EigenVectors[ii]);
		}

	//	Matrix multiply, C := alpha*A*B + beta*C
		double alpha = 1;
		double beta = 0;
		int ONE = 1;
		Array_1D <double> Product(n); //C Matrix

		dgemm(	"N", 				//Don't transpose A
				"N", 				//Don't transpose B
				&n, 				//Rows of A
				&ONE, 				//Columns of B (or C)
				&n, 				//Columns of A (or Rows of B)
				&alpha,
				cossq.grid, 		//A matrix
				&n, 				//Leading dimension of A
				GroundState.grid, 	//B Martrix
				&n, 				//Leading dimension of B
				&beta,
				Product.grid, 		//C matrix
				&n); 				//Leading dimension of C

		double cossq_alignment = ddot(&n, Product.grid, &ONE, GroundState.grid, &ONE);
	//	print_matrix("cos^2 theta",10,10,cossq.grid,cossq.Nx);

		if (EigenVectorsOut != 0)
		{
			for (int i = 0; i < n; ++i)
			{
				EigenOut << efield*LASERINTEN << " " << 2*i << " " << pow(abs(GroundState.get_elem(i)),2) << endl;
			}
			EigenOut << endl;
		}

		COSARRAY.grid[kk] = cossq_alignment;
		ENERGYARRAY.grid[kk] = EigenValues[0];

		delete [] EigenValues;
		delete [] EigenVectors;
}
	EigenOut.close();

	cout << endl << "Calculation Complete" << endl;
	ofstream AlignmentFile;
	AlignmentFile.open("output_data/AlignmentData.txt");
	for (ii=0; ii<ENERGYARRAY.Nx; ii++)
	{
		for (jj=0; jj<ENERGYARRAY.Ny; jj++)
		{
			AlignmentFile << ii << " " << jj << " " << EMF->get_elem(ii,jj) << " " << COSARRAY.get_elem(ii,jj) << " " << ENERGYARRAY.get_elem(ii,jj) << endl;
		}
		AlignmentFile << endl;
	}

	AlignmentFile.close();

/*******************************************
Calculate QM Trajectory using calculated PES
********************************************/

	double dt = ttotal/tsteps;
	dt /= 2.41e-2; //fs -> au
	ttotal /= 2.41e-2;
	cout << "Considering trajectory for " << dt*tsteps/1.89e7 << " ns." << endl;
	cout << endl << "Beginning quantum mechanical trajectory..." << endl;

	//Initialize fftw threading
	fftw_init_threads();
	fftw_plan_with_nthreads(Procs);

	//define step sizes and calculation arrays
	double xstep = Field_a/(LEN*Field_res); //Field Spacing in au
	double ystep = xstep;
	double pstep = M_PI*2.0/(xsize*xstep);
	double qstep = M_PI*2.0/(ysize*ystep);

	Array_1D <double> Xgrid(xsize);
	Xgrid.xinit = 0.0;
	Xgrid.xstep = xstep;
	Xgrid.fill_array();

	Array_1D <double> Ygrid(ysize);
	Ygrid.xinit = 0.0;
	Ygrid.xstep = ystep;
	Ygrid.fill_array();

	Array_1D <double> Pgrid(xsize);
	Pgrid.xinit = 0.0;
	Pgrid.xstep = pstep;
	Pgrid.fill_array();

	Array_1D <double> Qgrid(ysize);
	Qgrid.xinit = 0.0;
	Qgrid.xstep = qstep;
	Qgrid.fill_array();

	//Correct array ordering
	for (ii = xsize/2; ii < xsize; ii++)
	{
		Pgrid.grid[ii] = -1.0*Pgrid.grid[xsize - ii];
	}

	for (jj = ysize/2; jj < ysize; jj++)
	{
		Qgrid.grid[jj] = -1.0*Qgrid.grid[ysize - jj];
	}

	Array_2D <double> KE(xsize, ysize);
	KE.xstep = pstep;
	KE.ystep = qstep;
	for (ii = 0; ii < KE.Nx; ii++)
	{
		for (jj = 0; jj < KE.Ny; jj++)
		{
			KE.set_elem(ii,jj, (0.5/(amass1+amass2))*(pow(Pgrid.grid[ii],2) + pow(Qgrid.grid[jj],2)));
		}
	}

	cplx comp;

	Array_2D <cplx> KinetOp(KE.Nx,KE.Ny);
	for (ii=0; ii < KE.Nx*KE.Ny; ii++)
	{
		comp				= -0.5j * dt * KE.grid[ii];
		KinetOp.grid[ii] 	= exp(comp);
	}

	Array_2D <cplx> PotenOp(ENERGYARRAY.Nx,ENERGYARRAY.Ny);
	for (ii=0; ii < ENERGYARRAY.Nx*ENERGYARRAY.Ny; ii++)
	{
		comp 				= -1.0j * dt * ENERGYARRAY.grid[ii];
		PotenOp.grid[ii] 	= exp(comp);
	}

	//Set up the wavefunction array
	Array_2D <cplx> Wvfxn(xsize, ysize);
	Wvfxn.xstep = xstep;
	Wvfxn.ystep = ystep;

	fftw_plan forplan, backplan;
	forplan = fftw_plan_dft_2d(Wvfxn.Nx, Wvfxn.Ny, (fftw_complex *)Wvfxn.grid, (fftw_complex *)Wvfxn.grid, FFTW_FORWARD, FFTW_MEASURE);
	backplan = fftw_plan_dft_2d(Wvfxn.Nx, Wvfxn.Ny, (fftw_complex *)Wvfxn.grid, (fftw_complex *)Wvfxn.grid, FFTW_BACKWARD, FFTW_MEASURE);

	//Set initial wavefunction
	// for (ii = 0; ii < Wvfxn.Nx*Wvfxn.Ny; ii++)
	// {
	// 	Wvfxn.grid[ii] = 0.0;
	// }
	for (ii = 0; ii < Wvfxn.Nx; ii++)
	{
		for (jj = 0; jj < Wvfxn.Ny; jj++)
		{
			Wvfxn.set_elem(ii,jj,gaussian(Xgrid.grid[ii],Ygrid.grid[jj],xstep*InitialX*Field_res,ystep*InitialY*Field_res));
		}
	}
//	normalize_wxfxn_2D(Wvfxn);

	//Place the Particle at an initial location
	cout << (int)(InitialY*Field_res) << " " << (int)(InitialX*Field_res) << endl;
 	//Wvfxn.set_elem((int)(InitialY*Field_res),(int)(InitialX*Field_res), 1.0);
	int tindex = 0;
	do
	{
		if (tindex % dt_image_out == 0)
		{
			print_wvfxn(Wvfxn, Xgrid, Ygrid, tindex);
		}
		SplitOp2D_Step(Wvfxn,KinetOp,PotenOp,forplan,backplan);
		Wvfxn.time += dt;
		tindex++;
	} while (Wvfxn.time < ttotal);


	cout << endl << "Done. " << tindex << " total time steps." << endl;

/*********************
Print final commments
**********************/

	t = time(0);   // get time now
    now = localtime( & t );
	cout << "\n*********************************************" << endl;
	cout << "Code complete as of "<< asctime(now);
	cout << "*********************************************" << endl << endl;

	exit(0);
}

double gaussian(double x, double y, double xcen, double ycen)
{
	return exp(-1.0*(pow(x-xcen,2.0) + pow(y-ycen,2.0))/pow(50.0,2.0));
};


void print_matrix(char* desc, int m, int n, double* a, int lda)
{
	int ii, jj;
	cout << scientific << endl << desc << endl;
	for( ii = 0; ii < m; ii++ )
	{
		for( jj = 0; jj < n; jj++ )
		{
			cout << a[ii+jj*lda] << " ";
		}
		cout << endl;
	}
}

void check(int n)
{
	cout << "Check " << n << "." << endl;
}


