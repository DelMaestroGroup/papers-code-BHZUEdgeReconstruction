#include <sstream>
#include <cassert>
#include "tensors.h"
#include "Parameters_BHZ.h"
#include "random"

#ifndef MFParam_BHZ_class
#define MFParam_BHZ_class

class MFParam_BHZ
{

public:
	// Constructor for the SC class
	MFParam_BHZ(Parameters_BHZ &Parameters_BHZ__, mt19937_64 &Generator1__)
		: Parameters_BHZ_(Parameters_BHZ__), Generator1_(Generator1__)
	{
		Initialize();
		initialOrderParams();
		initialIonicPotential();
	}

	Parameters_BHZ &Parameters_BHZ_;

	void Initialize();
	double random1();
	void initialOrderParams();
	void initialIonicPotential();

	uniform_real_distribution<double> dis1_;
	mt19937_64 &Generator1_;

	int lx_, ly_, N_orbs_, N_spin_ = 2, N_cells_, wx_;
	int H_size_, half_size_;
	Mat_1_Complex_doub initialOPs_;
	Mat_1_doub initial_ionic;
	int ionic_size;
};


void MFParam_BHZ::Initialize(){
	lx_ = Parameters_BHZ_.Lx_;
	ly_ = Parameters_BHZ_.Ly_;
	wx_ = Parameters_BHZ_.W_;
	N_orbs_ = Parameters_BHZ_.N_Orbs;

	N_cells_ = Parameters_BHZ_.Total_Cells;
	half_size_ = lx_ * ly_ * N_orbs_;
	H_size_ = 2 * half_size_;

	initialOPs_.resize(H_size_);
	//dis1_ = std::uniform_real_distribution<double>(0.0, 1.0);
	if(Parameters_BHZ_.numberofSoftedges_==1){	ionic_size = lx_ - wx_ + 1;		}
    else{	ionic_size = lx_ - 2*wx_ + 2;	}
	initial_ionic.resize(ionic_size);
}

double MFParam_BHZ::random1(){
	return dis1_(Generator1_);
}

void MFParam_BHZ::initialOrderParams(){

	if(!Parameters_BHZ_.read_OP_){
		Mat_4_Complex_doub OP_dummy;
		OP_dummy.resize(lx_);
		for(int rx=0;rx<lx_;rx++){
			OP_dummy[rx].resize(ly_);
			for(int ry=0;ry<ly_;ry++){
				OP_dummy[rx][ry].resize(N_orbs_);
				for(int orb=0;orb<N_orbs_;orb++){
					OP_dummy[rx][ry][orb].resize(N_spin_);
					for(int spin=0;spin<N_spin_;spin++){
						OP_dummy[rx][ry][orb][spin]=Zero_Complex;
					}
				}
			}
		}

		string OP_file_ = "Initial_Order_Parameter_file.txt";
		ofstream OP_file_out(OP_file_.c_str());

		Mat_2_Complex_doub PBC_OPs;
		PBC_OPs.resize(N_orbs_);
		for(int orb=0;orb<N_orbs_;orb++){
			PBC_OPs[orb].resize(N_spin_);
			for(int spin=0;spin<N_spin_;spin++){
				PBC_OPs[orb][spin] = random1()*One_Complex;
			}
		}

		int r;
		for (int rx = 0; rx < lx_; rx++){
			for (int ry = 0; ry < ly_; ry++){
				for (int orb = 0; orb < N_orbs_; orb++){
					for(int spin=0;spin<N_spin_;spin++){

						if(Parameters_BHZ_.PBC_X==true && Parameters_BHZ_.PBC_Y==true){
							//Translation symmetry along-x and -y
							OP_dummy[rx][ry][orb][spin] = PBC_OPs[orb][spin];
						}
						else if(Parameters_BHZ_.canonical_==true && Parameters_BHZ_.restricted_updn_symm==true){
							if(spin==0){
								OP_dummy[rx][ry][orb][0] = random1()*One_Complex;
								OP_dummy[rx][ry][orb][1] = OP_dummy[rx][ry][orb][0];
							}
						}
/*						else if(Parameters_BHZ_.PBC_X==false && Parameters_BHZ_.PBC_Y==true){
							if(Parameters_BHZ_.canonical_==false && Parameters_BHZ_.numberofSoftedges_==1){
								OP_dummy[rx][ry][orb][spin] = random1()*One_Complex;
							}
							else{
								//Mirror symmetry along-x:
								if(lx_%2==0 && rx<lx_/2){
									OP_dummy[rx][ry][orb][spin] = random1()*One_Complex;
									OP_dummy[lx_-1-rx][ry][orb][spin] = OP_dummy[rx][ry][orb][spin];		
								}
								else if(lx_%2==1 && rx<(lx_+1)/2){
									OP_dummy[rx][ry][orb][spin] = random1()*One_Complex;
									OP_dummy[lx_-1-rx][ry][orb][spin] = OP_dummy[rx][ry][orb][spin];
								}
							}
						}
*/						else{
							OP_dummy[rx][ry][orb][spin] = random1()*One_Complex;
						}

						r = spin * half_size_ + (orb + 2 * ry + 2 * ly_ * rx);
						initialOPs_[r] = OP_dummy[rx][ry][orb][spin];

						OP_file_out<<r<<"  "<< initialOPs_[r].real()<<"	"<<initialOPs_[r].imag()<< endl;
					}
				}
			}
		}
		
	}
	else{
		string OP_file_name=Parameters_BHZ_.input_OP_file;
		ifstream OP_file_in(OP_file_name.c_str());
		if (!OP_file_in.is_open()) {
			cerr << "Error: Could not open file " << OP_file_name << endl;
			return;
		}

		int temp_r;
		double temp_OP_real, temp_OP_imag;
		Mat_1_doub input_OP_real,input_OP_imag;

		string line;
		while(getline(OP_file_in,line) ){
			stringstream line_ss;
			line_ss << line;
			line_ss>>temp_r>>temp_OP_real>>temp_OP_imag;
			input_OP_real.push_back(temp_OP_real);
			input_OP_imag.push_back(temp_OP_imag);
		}

		assert(input_OP_real.size()==H_size_);

		int r,cell;
		for (int ind=0; ind<H_size_;ind++){
			initialOPs_[ind].real(input_OP_real[ind]);
			initialOPs_[ind].imag(input_OP_imag[ind]);
		}
		input_OP_real.clear();
		input_OP_imag.clear();

	}
	
}

void MFParam_BHZ::initialIonicPotential(){
	for(int rx=0;rx<ionic_size;rx++){
		initial_ionic[rx] = Parameters_BHZ_.no_val;
		if(Parameters_BHZ_.numberofSoftedges_==1 && rx==ionic_size-1){
			initial_ionic[rx] += -Parameters_BHZ_.helical_edge_pot;
		}
		//initial_ionic[rx] = 0.0;
	}
}

#endif
