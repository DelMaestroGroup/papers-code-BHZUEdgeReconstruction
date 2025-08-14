#ifndef Parameters_BHZ_class
#define Parameters_BHZ_class
#include <iostream>
#include <fstream>
#include "tensors.h"

//This class will read from the **input file**

class Parameters_BHZ{

	public:
		//Declaring system params:---->
		bool PBC_X, PBC_Y;
		int Lx_, Ly_, N_Orbs, W_;
		int Total_Sites, Total_Cells, Ham_Size;

		//Declaring model params:----->
		double A_val, B_val, M_val, D_val, C_val;
		double no_val, Vo_val, fill_;
		double Hfield;
		int total_particles_;
		double mu_fixed, temp_, beta;
		//double mu_canonical, mu_;
		double U_, JHbyU_, JH_, Up_;
		double V_nn;
		bool canonical_;
		int numberofSoftedges_;

		double helical_edge_pot;

		//Declaring self-consistency params:---->
		bool Simple_mixing, Broyden_mixing;
		double Conv_err, alpha_;
		int Seed_, Max_iters;

		//Observables details:--->
		bool get_dos, get_ldoe, get_wave_fn, get_SScorr;
		bool get_sc, get_Akxw, get_Akyw, get_ldos;
		double w_min, w_max, dw_, eta_;
		int state_, npthreads_;
		//----------------------->

        bool read_OP_;
		string input_OP_file;

		bool restricted_updn_symm=true;
		bool sitewiseIonic;
        
		void Initialize(std::string input_file);
		double matchstring(std::string file, std::string match);
		std::string matchstring2(std::string file, std::string match);
};


void Parameters_BHZ::Initialize(string input_file){

	string PBC_X_string,PBC_Y_string;
	string pbc_x_out, pbc_y_out;

	cout<<"-------------------------------------------\n"
		<<"Reading the input file = "<<input_file<<"\n"
		<<"-------------------------------------------\n"<<endl;

	PBC_X_string=matchstring2(input_file,"PBC_X");
	if(PBC_X_string=="True")
	{				PBC_X=true;		pbc_x_out="PBC";	}
	else if(PBC_X_string == "False")
	{				PBC_X=false;	pbc_x_out="OBC";	}
	
	PBC_Y_string=matchstring2(input_file,"PBC_Y");
	if(PBC_Y_string=="True")
	{				PBC_Y=true;		pbc_y_out="PBC";	}
	else if(PBC_Y_string == "False")
	{				PBC_Y=false;	pbc_y_out="OBC";	}

    Lx_ = int(matchstring(input_file, "Cells_X"));
    Ly_ = int(matchstring(input_file, "Cells_Y"));
    N_Orbs = int(matchstring(input_file, "Total_Orbs"));
    W_ = int(matchstring(input_file, "Softwall_Width"));

	Total_Cells = Lx_*Ly_;
	Total_Sites = N_Orbs*Total_Cells;
	Ham_Size = 2*Total_Sites;

	cout<<"Boundary conditions = "<<pbc_x_out<<"x"<<pbc_y_out<<"\n"
			<<"Total size of the Hamiltonian = "<<Ham_Size<<"x"<<Ham_Size<<endl;
	//----------------------------------------------------------------------//

	A_val = matchstring(input_file, "Inter_Orb_Hopping_A");
	B_val = matchstring(input_file, "Intra_Orb_Hopping_B");
	D_val = matchstring(input_file, "Intra_Orb_Hopping_D");
	M_val = matchstring(input_file, "Gap_Parameter_M");
	C_val = matchstring(input_file, "Gap_Parameter_C");
	Hfield = matchstring(input_file, "Mag_Field");

	U_ = matchstring(input_file, "U_Onsite");
	JHbyU_ = matchstring(input_file, "JHbyU_Ratio");
	JH_ = 1.0*U_*JHbyU_;
	Up_ = U_ - 2.0 * JH_;
	V_nn = matchstring(input_file, "VNN_Hartree");

	string Ensemble_str, Ens_out;
	Ensemble_str = matchstring2(input_file, "Ensemble");
	if (Ensemble_str == "CE")
	{       canonical_ = true;              Ens_out = "Canonical";		}
	else
	{       canonical_ = false;             Ens_out = "Grand-Canonical";}

	fill_ = matchstring(input_file, "Filling");
	temp_ = matchstring(input_file, "Temperature");
	beta = 1.0 / (1.0 * temp_);

	numberofSoftedges_ = int(matchstring(input_file, "Number_of_Softwalls"));

	W_ = int(matchstring(input_file, "Softwall_Width"));
	Vo_val = matchstring(input_file, "Confining_Potential_Vo");
	mu_fixed = matchstring(input_file, "Fixed_mu");
	no_val = (3.0*U_-5.0*JH_)/2.0;
	helical_edge_pot = matchstring(input_file, "Helical_Edge_Pot");
	
	string ionic_str;
	ionic_str = matchstring2(input_file, "siteBasedIonic");
	if(ionic_str=="True"){	sitewiseIonic=true;	}
	else{	sitewiseIonic=false;	}
	

	cout << "(A, B, D, M) = (" << A_val << " , " << B_val << " , " << D_val << " , " << M_val << ")" << "\n"
		 << "(U, Up, JH) = (" << U_ << " , " << Up_ << " , " << JH_ << ")" << "\n"
                 << "Performing " << Ens_out << "-ensemble" << endl;
	
	if (canonical_ == true){
		total_particles_ = (int) (fill_*Ham_Size + 0.01);
		cout << "filling = " << fill_ << endl;
	}
	else if (canonical_ == false){
		cout << "Confining-width = " << W_ << endl;
		cout << "Confining-potential = " << Vo_val << endl;
		cout << "mu fixed = " << mu_fixed << endl;
		if(sitewiseIonic==false){
			cout << "site-based Ionic is OFF" << endl;
		}
		else{
			cout << "site-based Ionic is ON" << endl;
		}
	}
	//----------------------------------------------------------------------//
	
	string mixing_, mixing_out;
	mixing_ = matchstring2(input_file,"Simple_Mixing");
	if(mixing_=="True"){
		Simple_mixing=true;		Broyden_mixing=false;
		mixing_out="Simple_Mixing";
	}
	else{
		Simple_mixing=false;	Broyden_mixing=true;
		mixing_out="Broyden_Mixing";
	}

	Conv_err = matchstring(input_file, "Convergence_Error");
	alpha_ = matchstring(input_file, "alpha_OP");

	Seed_ = int(matchstring(input_file, "Random_Seed"));
	Max_iters = int(matchstring(input_file, "Max_Iterations"));
	cout<<"Self-consistency solver = "<<mixing_out<<endl;
	//----------------------------------------------------------------------//
	
	string ldoe_str, dos_str, sc_str, wave_fn_str, Akxw_str, Akyw_str, ss_corr_str, ldos_str;
	ldoe_str = matchstring2(input_file, "Calculate_LDOE");
	if (ldoe_str == "True")
	{
		get_ldoe = true;
		cout << "Measuring LDOE" << endl;
	}
	else {  get_ldoe = false;       }

	sc_str = matchstring2(input_file, "Calculate_SC");
	if (sc_str == "True")
	{
		get_sc = true;
		cout << "Measuring Spin Currents" << endl;
	}
	else {  get_sc = false;         }

	wave_fn_str = matchstring2(input_file, "Calculate_wave_fn");
	if (wave_fn_str == "True")
	{
		get_wave_fn = true;
		state_ = int(matchstring(input_file, "State_val"));
		cout << "Measuring wave functions" << endl;
	}
	else {  get_wave_fn = false;    }

	dos_str = matchstring2(input_file, "Calculate_DOS");
	if (dos_str == "True")
	{
		get_dos = true;
		cout << "Measuring DOS" << endl;
	}
	else {  get_dos = false;        }

	Akxw_str = matchstring2(input_file, "Spectral_Akxw");
	if (Akxw_str == "True")
	{
		get_Akxw = true;
		cout << "Measuing spectral function A(kx,w)" << endl;
	}
	else {  get_Akxw = false;       }

	Akyw_str = matchstring2(input_file, "Spectral_Akyw");
    if (Akyw_str == "True")
	{
		get_Akyw = true;
		cout << "Measuing spectral function A(ky,w)" << endl;
	}
	else {  get_Akyw = false;       }

	ss_corr_str = matchstring2(input_file, "Calculate_SS_corr");
	if (ss_corr_str == "True")
	{
		get_SScorr = true;
		cout << "Measuing spin-spin correlations" << endl;
	}
	else {  get_SScorr = false;       }

	ldos_str = matchstring2(input_file, "Calculate_LDOS_along_x");
	if (ldos_str == "True")
	{
		get_ldos = true;
		cout << "Measuring LDOS along-x" << endl;
	}
	else {  get_ldos = false;       }

	w_min = matchstring(input_file, "omega_min");
	w_max = matchstring(input_file, "omega_max");
	dw_ = matchstring(input_file, "d_omega");
	eta_ = matchstring(input_file, "broadening");

	//----------------------------------------------------------------------//
	string read_op_str;
	read_op_str = matchstring2(input_file, "Read_OP_from_file");
	if (read_op_str == "True"){
		read_OP_ = true;
		input_OP_file = matchstring2(input_file, "OP_Filename");
	}
	else{
		read_OP_ = false;
	}
	//----------------------------------------------------------------------//
	

	cout << "-------------------------------------------\n"
         << "Finish Reading the input file" << "\n"
		 << "-------------------------------------------\n"<< endl;				 
}

double Parameters_BHZ::matchstring(string file, string match){

	std::string test, line;
    std::ifstream inputfile(file);
    double amount;

    bool pass = false;

    while (std::getline(inputfile, line)) {
        std::istringstream iss(line);

        if (std::getline(iss, test, '=') && !pass) {
            if (iss >> amount && test == match) {
                pass = true;
            }
            if (pass) {
                break;
            }
        }
    }

    if (!pass) {
        throw std::invalid_argument("Missing argument in the input file: " + match);
    }

    return amount;
}


string Parameters_BHZ::matchstring2(string file, string match){
	std::string line;
    std::ifstream inputfile(file);
    std::string amount;
    int offset;

	/* Referenced from https://cplusplus.com/forum/beginner/121556/
			 * https://stackoverflow.com/questions/12463750/ */
	if(inputfile.is_open()) {
        while(std::getline(inputfile, line)) {
            if((offset = line.find(match, 0)) != std::string::npos) {
                amount = line.substr(offset + match.length() + 1);
                break; // Break early if the match is found
            }
        }
        inputfile.close();
    } else {
        std::cerr << "Unable to open input file while in the Parameter class." << std::endl;
        return "";
    }

    //std::cout << match << " = " << amount << std::endl;
    return amount;
}


#endif
