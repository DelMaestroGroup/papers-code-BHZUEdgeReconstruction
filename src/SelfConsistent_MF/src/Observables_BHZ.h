#include "tensors.h"
#include "Parameters_BHZ.h"
#include "Hamiltonian_BHZ.h"

#ifndef Observables_BHZ_class
#define Observables_BHZ_class
#define PI acos(-1.0)

class Observables_BHZ{

    public:
    Observables_BHZ(Parameters_BHZ &Parameters_BHZ__, Hamiltonian_BHZ &Hamiltonian_BHZ__): 
        Parameters_BHZ_(Parameters_BHZ__),Hamiltonian_BHZ_(Hamiltonian_BHZ__){
            Initialize();
    }

    Parameters_BHZ &Parameters_BHZ_;
    Hamiltonian_BHZ &Hamiltonian_BHZ_;

    void Initialize();
    double fermifunction(double en_, double mu_val);
    double chemicalpotential(int particles_, double muin);
    double chemicalpotential2(int particles_);

    void getNewOPs();
    pair<double,double> getEnergies();
    double getOPError();
    void updateOrderParams(int iter);
    void updateIonicPotential();

    void calculateDOS();
    complex<double> calculateAmplitude(int state1, int pos1, int state2, int pos2);
    void calculateLDOE();
    void calculateSpinCurrents();
    void calculateAkxw();
    void calculateAkyw();
    void calculateWaveFunctions();
    void calculateSpinSpinCorr();
    void calculateLDOS();
	
    int lx_, ly_, H_size, MHS;
    int N_cells_,N_orbs_, N_spin_;
    
    double one_by_PI_=1.0/(1.0*PI);
    double w_min,w_max,dw,eta;
    int w_size;

    double mu;

    Mat_1_Complex_doub new_OPs_;
    Mat_1_doub new_ionic;
    int ionic_size;
    double ionic_relaxation;

    //For Broyden:
    int OP_size;
    Mat_2_doub B_m, B_mp1;
    Mat_1_doub del_X,del_F;
    Mat_1_doub F_m,F_mp1;
};

void Observables_BHZ::Initialize(){

    lx_ = Parameters_BHZ_.Lx_;
    ly_ = Parameters_BHZ_.Ly_;
    
    N_cells_ = Parameters_BHZ_.Total_Cells;
    N_orbs_ = Parameters_BHZ_.N_Orbs;
    N_spin_ = 2;
    H_size = Parameters_BHZ_.Ham_Size;
    MHS = Hamiltonian_BHZ_.half_size_;

    w_min = Parameters_BHZ_.w_min;
    w_max = Parameters_BHZ_.w_max;

    dw = Parameters_BHZ_.dw_;
    eta = Parameters_BHZ_.eta_;

    w_size = (int) ( (w_max - w_min)/dw );

/*    mu = Parameters_BHZ_.mu_;
    cout<< mu <<endl;

    if(canonical_==true){
        mu_=mu_canonical;
        cout<<"Canonically calculated mu = "<<mu_<<endl;
    }
    else{
        mu_=mu_fixed;
        cout<<"Grand-canonical fixed mu = "<<mu_<<endl;
    }
*/

    new_OPs_.resize(H_size);

    OP_size = Hamiltonian_BHZ_.OPs_.size();
    B_m.resize(OP_size);                B_mp1.resize(OP_size);
    for(int i=0;i<OP_size;i++){
        B_m[i].resize(OP_size);         B_mp1[i].resize(OP_size);
    }
    del_X.resize(OP_size);              del_F.resize(OP_size);
    F_m.resize(OP_size);                F_mp1.resize(OP_size);

    ionic_size = Hamiltonian_BHZ_.ionic_size;
    new_ionic.resize(ionic_size);
    ionic_relaxation = 0.01;
}

double Observables_BHZ::chemicalpotential2(int particles_){
    double mu_temp, eps_, dmu_by_dN, N_temp, dmu_by_dN_min, Ne_;
    eps_=1e-2;
    mu_temp=Hamiltonian_BHZ_.evals_[0];
    N_temp=100000;
    Ne_=1.0*particles_;

    int iters=0;

    dmu_by_dN = 0.01*( Hamiltonian_BHZ_.evals_[H_size-1] - Hamiltonian_BHZ_.evals_[0] )*( 1.0/(1.0*H_size) );
    dmu_by_dN_min = 0.0001*( Hamiltonian_BHZ_.evals_[H_size-1] - Hamiltonian_BHZ_.evals_[0] )*( 1.0/(1.0*H_size) );

    while( abs(N_temp - Ne_) > eps_){
        N_temp=0;
        for(int i=0;i<H_size;i++){
            N_temp += fermifunction(Hamiltonian_BHZ_.evals_[i],mu_temp);
        }
        mu_temp = mu_temp + dmu_by_dN*(Ne_-N_temp);
        iters++;
        if(iters%1000==0){
            dmu_by_dN = (1000.0/(10.0*iters))*dmu_by_dN;
            dmu_by_dN = max(dmu_by_dN, dmu_by_dN_min);
        }
    }
    
    cout<<"Calculated mu = "<<mu_temp<<endl;
    cout<<"Calculated # of particles = "<<N_temp<<endl;
    return mu_temp;

}

double Observables_BHZ::chemicalpotential(int particles_, double muin){

    double mu_out;
    double n1,N;

    double dMubydN;
    double nstate = Hamiltonian_BHZ_.evals_.size();
    dMubydN = 0.0005*(Hamiltonian_BHZ_.evals_[nstate-1] - Hamiltonian_BHZ_.evals_[0])/nstate;
    N=particles_;
    mu_out = muin;
    bool converged=false;
    int final_i;

    double mu1, mu2;
    double mu_temp = muin;

    mu1=Hamiltonian_BHZ_.evals_[0]- (5.0/Parameters_BHZ_.beta);
    mu2=Hamiltonian_BHZ_.evals_[nstate-1] + (5.0/Parameters_BHZ_.beta);
    for(int i=0;i<400000;i++){
        n1=0.0;
        for(int j=0;j<nstate;j++){
            n1+=double(1.0/( exp( (Hamiltonian_BHZ_.evals_[j]-mu_temp)*Parameters_BHZ_.beta ) + 1.0));
        }
        //cout <<"i  "<< i << "  n1  " << n1 << "  mu  " << mu_out<< endl;
        if(abs(N-n1)<double(0.000001)){
            //cout<<abs(N-n1)<<endl;
            converged=true;
            break;
        }
        else {
            if(n1 >N){
                mu2=mu_temp;
                mu_temp=0.5*(mu1 + mu_temp);
            }
            else{
                mu1=mu_temp;
                mu_temp=0.5*(mu2 + mu_temp);
            }
        }
        //cout<<"mu_temp = "<<mu_temp<<"   "<<mu1<<"    "<<mu2<<"   "<<eigs_[nstate-1]<<"  "<<eigs_[0]<<"  "<<n1<<endl;
    }
    if(!converged){
        cout<<"mu_not_converged, N = "<<n1<<endl;
    }
    else{
        cout<<"mu converged, N = "<<n1<<endl;
    }
    
    mu_out = mu_temp;

    return mu_out;
}

double Observables_BHZ::fermifunction(double en_, double mu_val){
	double ffn, temp;
	temp = Parameters_BHZ_.temp_;
	ffn = ((1.0) / (1.0 + exp((en_ - mu_val) / temp)));
	return ffn;
}

void Observables_BHZ::getNewOPs(){
	//Calculating new order-parameters:
    int r;
    double mu_val=mu;
	complex<double> value;

	for (int rx = 0; rx < lx_; rx++){
		for (int ry = 0; ry < ly_; ry++){
			for (int orb = 0; orb < N_orbs_; orb++){
				for (int spin = 0; spin < N_spin_; spin++){
                    r = Hamiltonian_BHZ_.makeIndex(rx, ry, orb, spin);

                    value = Zero_Complex;

					for (int n = 0; n < H_size; n++){
						value += One_Complex * calculateAmplitude(n, r, n, r) * fermifunction(Hamiltonian_BHZ_.evals_[n],mu_val);
					}
                    //cout<<mu_val<<endl;

					new_OPs_[r] = value;
				}
			}
		}
	}
}

double Observables_BHZ::getOPError(){
    // Obtaining the classical and quantum energies:
	double error_OP;
	error_OP = 0.0;

    int r;
    for (int rx = 0; rx < lx_; rx++){
		for (int ry = 0; ry < ly_; ry++){
            for (int orb = 0; orb < N_orbs_; orb++){
                for (int spin = 0; spin < N_spin_; spin++){

                    r = Hamiltonian_BHZ_.makeIndex(rx,ry,orb,spin);

                    error_OP += abs(new_OPs_[r] - Hamiltonian_BHZ_.OPs_[r]) *
							abs(new_OPs_[r] - Hamiltonian_BHZ_.OPs_[r]);
                }
			}
		}
	}
	error_OP = sqrt(error_OP);

	return error_OP;
}

pair<double, double> Observables_BHZ::getEnergies(){
    //Calculating classical and quantum energies:
	double E_class, E_quant;
	E_class = 0.0;
	E_quant = 0.0;

    //Classical
    int r_up, r_dn;
    for (int rx = 0; rx < lx_; rx++){
		for (int ry = 0; ry < ly_; ry++){
            for (int orb = 0; orb < N_orbs_; orb++){

                r_up = Hamiltonian_BHZ_.makeIndex(rx, ry, orb, 0);
                r_dn = Hamiltonian_BHZ_.makeIndex(rx, ry, orb, 1);

                // Hartree hubbard classical energy:"-U<n_{r,a,up}><n_{r,a,dn}>"
                E_class += -1.0 * Parameters_BHZ_.U_ * Hamiltonian_BHZ_.OPs_[r_up].real() * Hamiltonian_BHZ_.OPs_[r_dn].real();

            }
		}
    }
    int r_s_up, r_s_dn, r_p_up, r_p_dn;
    for (int rx = 0; rx < lx_; rx++){
		for (int ry = 0; ry < ly_; ry++){
            r_s_up = Hamiltonian_BHZ_.makeIndex(rx, ry, 0, 0);
            r_s_dn = Hamiltonian_BHZ_.makeIndex(rx, ry, 0, 1);
            r_p_up = Hamiltonian_BHZ_.makeIndex(rx, ry, 1, 0);
            r_p_dn = Hamiltonian_BHZ_.makeIndex(rx, ry, 1, 1);
            
            // Hartree orbital-repulsion classical energy:"-(U'-0.5J)<n_{r,s}><n_{r,p}>"
            E_class += -1.0 * (Parameters_BHZ_.Up_ - 0.5 * Parameters_BHZ_.JH_) 
                            * (Hamiltonian_BHZ_.OPs_[r_s_up].real() + Hamiltonian_BHZ_.OPs_[r_s_dn].real()) 
                            * (Hamiltonian_BHZ_.OPs_[r_p_up].real() + Hamiltonian_BHZ_.OPs_[r_p_dn].real());
        }
	}

    int rpx_s_up,rpx_s_dn,rpx_p_up,rpx_p_dn;
    int rpy_s_up,rpy_s_dn,rpy_p_up,rpy_p_dn;

    for (int rx = 0; rx < lx_; rx++){
		for (int ry = 0; ry < ly_; ry++){
            r_s_up = Hamiltonian_BHZ_.makeIndex(rx, ry, 0, 0);
            r_s_dn = Hamiltonian_BHZ_.makeIndex(rx, ry, 0, 1);
            r_p_up = Hamiltonian_BHZ_.makeIndex(rx, ry, 1, 0);
            r_p_dn = Hamiltonian_BHZ_.makeIndex(rx, ry, 1, 1);

            double n_r = Hamiltonian_BHZ_.OPs_[r_s_up].real() + Hamiltonian_BHZ_.OPs_[r_s_dn].real()
                        + Hamiltonian_BHZ_.OPs_[r_p_up].real() + Hamiltonian_BHZ_.OPs_[r_p_dn].real();

            // Hartree nearest-neighbour classical energy:"-V(<n_{r}n_{r+x}>+<n_r><n_{r+y}>)"
            if(rx <lx_-1){
                rpx_s_up = Hamiltonian_BHZ_.makeIndex(rx+1, ry, 0, 0);
                rpx_s_dn = Hamiltonian_BHZ_.makeIndex(rx+1, ry, 0, 1);
                rpx_p_up = Hamiltonian_BHZ_.makeIndex(rx+1, ry, 1, 0);
                rpx_p_dn = Hamiltonian_BHZ_.makeIndex(rx+1, ry, 1, 1);

                double n_rpx = Hamiltonian_BHZ_.OPs_[rpx_s_up].real() + Hamiltonian_BHZ_.OPs_[rpx_s_dn].real()
                                + Hamiltonian_BHZ_.OPs_[rpx_p_up].real() + Hamiltonian_BHZ_.OPs_[rpx_p_dn].real();

                E_class += -1.0 * Parameters_BHZ_.V_nn * n_r * n_rpx;
            }
            else if(rx == lx_-1 && Parameters_BHZ_.PBC_X){
                rpx_s_up = Hamiltonian_BHZ_.makeIndex(0, ry, 0, 0);
                rpx_s_dn = Hamiltonian_BHZ_.makeIndex(0, ry, 0, 1);
                rpx_p_up = Hamiltonian_BHZ_.makeIndex(0, ry, 1, 0);
                rpx_p_dn = Hamiltonian_BHZ_.makeIndex(0, ry, 1, 1);

                double n_rpx = Hamiltonian_BHZ_.OPs_[rpx_s_up].real() + Hamiltonian_BHZ_.OPs_[rpx_s_dn].real()
                                + Hamiltonian_BHZ_.OPs_[rpx_p_up].real() + Hamiltonian_BHZ_.OPs_[rpx_p_dn].real();

                E_class += -1.0 * Parameters_BHZ_.V_nn * n_r * n_rpx;
            }
            if(ry<ly_-1){
                rpy_s_up = Hamiltonian_BHZ_.makeIndex(rx, ry+1, 0, 0);
                rpy_s_dn = Hamiltonian_BHZ_.makeIndex(rx, ry+1, 0, 1);
                rpy_p_up = Hamiltonian_BHZ_.makeIndex(rx, ry+1, 1, 0);
                rpy_p_dn = Hamiltonian_BHZ_.makeIndex(rx, ry+1, 1, 1);

                double n_rpy = Hamiltonian_BHZ_.OPs_[rpy_s_up].real() + Hamiltonian_BHZ_.OPs_[rpy_s_dn].real()
                                + Hamiltonian_BHZ_.OPs_[rpy_p_up].real() + Hamiltonian_BHZ_.OPs_[rpy_p_dn].real();

                E_class += -1.0 * Parameters_BHZ_.V_nn * n_r * n_rpy;
            }
            else if(ry==ly_-1 && Parameters_BHZ_.PBC_Y){
                rpy_s_up = Hamiltonian_BHZ_.makeIndex(rx, 0, 0, 0);
                rpy_s_dn = Hamiltonian_BHZ_.makeIndex(rx, 0, 0, 1);
                rpy_p_up = Hamiltonian_BHZ_.makeIndex(rx, 0, 1, 0);
                rpy_p_dn = Hamiltonian_BHZ_.makeIndex(rx, 0, 1, 1);

                double n_rpy = Hamiltonian_BHZ_.OPs_[rpy_s_up].real() + Hamiltonian_BHZ_.OPs_[rpy_s_dn].real()
                                + Hamiltonian_BHZ_.OPs_[rpy_p_up].real() + Hamiltonian_BHZ_.OPs_[rpy_p_dn].real();

                E_class += -1.0 * Parameters_BHZ_.V_nn * n_r * n_rpy;
            }
        }
    }

    //Quantum
	for (int n = 0; n < Hamiltonian_BHZ_.evals_.size(); n++){
		E_quant += Hamiltonian_BHZ_.evals_[n] * fermifunction(Hamiltonian_BHZ_.evals_[n], mu);
	}

	return make_pair(E_class, E_quant);
}

void Observables_BHZ::updateOrderParams(int iter){
	double alpha_OP = Parameters_BHZ_.alpha_;

	// Updating order parameters based on Simple mixing:----->
	if (Parameters_BHZ_.Simple_mixing == true){
		//      For simple (linear) mixing, we use:----->
		//      X^{m+1}_{in} = (1-alpha)*X^{m}_{in} + alpha*X^{m}_{out}
        //or    X^{m+1}_{in} = X^{m}_{in} + alpha*F^{m}
        //where, F^{m} = X^{m}_{out} - X^{m}_{in} 

        int r;        
		for (int rx = 0; rx < lx_; rx++){
            for (int ry = 0; ry < ly_; ry++){
                for (int orb = 0; orb < N_orbs_; orb++){
                    for (int spin = 0; spin < N_spin_; spin++){
                        r = Hamiltonian_BHZ_.makeIndex(rx, ry, orb, spin);

                        new_OPs_[r] = (1.0 - alpha_OP) * Hamiltonian_BHZ_.OPs_[r] + alpha_OP * new_OPs_[r];
                    }
				}
			}
		}

	}

	//  Updating order parameters based on Broyden mixing:----->
    //  Details of this algorithm can be found in: arXiv:0805.4446 (2008)!!
	if (Parameters_BHZ_.Broyden_mixing == true){

		/*For Broyden mixing, we use the multidimensional quasi-Newton-Raphson method:->    
            X^{m+1}_{in} = X^{m}_{in} - B^{m}*F^{m};
    
            where, based on Sherman-Morrison formula:
            B^{m+1} = B^{m} + [ |d(X_{in})> - B^{m}|dF> ] * <d(X_{in})|B^{m};
                                 ----------------------
                                  <d(X_{in})|B^{m}|dF>

            with B^{0} = alpha*I;   and, |dF> = F^{m+1} - F^{m};
            |d(X_{in})> = X^{m+1}_{in} - X^{m}_{in};
        */

        Mat_1_doub Left_vec,Right_vec;
        Left_vec.resize(OP_size);       Right_vec.resize(OP_size);
        double amplitude;

		if (iter == 0){
			// Initializing the approximated inverse Jacobian matrix (B):->
			for(int i = 0; i < OP_size; i++){
				for (int j = 0; j < OP_size; j++){
					if (i == j){
						B_mp1[i][j] = -alpha_OP;
					}
                    else{
                        B_mp1[i][j] = 0.0;
                    }
				}
                //Obtaining the F-vector for 0th-iteration:->
                F_mp1[i] = (new_OPs_[i] - Hamiltonian_BHZ_.OPs_[i]).real();
			}
            
            //Calculating difference in OPs:->
            for(int i=0;i<OP_size;i++){
                del_X[i]=0.0;
                for (int j = 0; j < OP_size; j++){
                    del_X[i] += -B_mp1[i][j]*F_mp1[j];
                }
                //Getting new estimates of order parameters:->
                new_OPs_[i] = Hamiltonian_BHZ_.OPs_[i] + del_X[i]*One_Complex;
            }

            //Saving F_m and B_m for next iteration:->         
            F_m=F_mp1;
            B_m=B_mp1;
            
		}
        else{
            
            //Obtaining new F-vetor:->
            for(int i=0;i<OP_size;i++){
                F_mp1[i] = (new_OPs_[i] - Hamiltonian_BHZ_.OPs_[i]).real();
            }
            

            //Calculating del_F vector for calculating new B:->
            for(int i=0;i<OP_size;i++){
                del_F[i] = F_mp1[i] - F_m[i];
            }

            //Calculating left and right vectors for getting the new B:->
            for(int i=0;i<OP_size;i++){
                Left_vec[i]=0.0;
                Right_vec[i]=0.0;
                for(int j=0; j<OP_size;j++){
                    Left_vec[i] += -B_m[i][j]*del_F[j];
                    Right_vec[i] += B_m[j][i]*del_X[j];
                }
            }
/*            for(int i=0;i<OP_size;i++){
                cout<<Left_vec[i]<<endl;
            }
            cout<<endl;
            for(int i=0;i<OP_size;i++){
                cout<<Right_vec[i]<<endl;
            }
            cout<<endl;
*/
            amplitude=0.0;
            for(int i=0;i<OP_size;i++){
                amplitude +=-del_X[i]*Left_vec[i];
            //    cout<<amplitude<<endl;
            }

            for(int i=0;i<OP_size;i++){
                Left_vec[i] = del_X[i]+Left_vec[i];
            }

            for(int i=0;i<OP_size;i++){
                for(int j=0;j<OP_size;j++){
                    B_mp1[i][j] = B_m[i][j] + (Left_vec[i]*Right_vec[j])/amplitude;
                }
            }

            //Calculating difference in OPs:->
            for(int i=0;i<OP_size;i++){
                del_X[i]=0.0;
                for (int j = 0; j < OP_size; j++){
                    del_X[i] += -B_mp1[i][j]*F_mp1[j];
                }
                //Getting new estimates of order parameters:->
                new_OPs_[i] = Hamiltonian_BHZ_.OPs_[i] + del_X[i]*One_Complex;
            }

            //Saving F_m and B_m for next iteration:->         
            F_m=F_mp1;
            B_m=B_mp1;
        }
	}

	
}

void Observables_BHZ::updateIonicPotential(){
    for(int rx=0;rx<ionic_size;rx++){
        new_ionic[rx] = 0.0;
    }

    Mat_1_doub avg_density(ionic_size,0.0);
    complex<double> temp_den;

    //For numberofSoftedges==1: rx_bulk = [0, ... lx_-wx_]
    //For numberofSoftedges==2: rx_bulk = [wx_-1, ... lx_-wx_]
    int rx_start, rx_end;
    if(Parameters_BHZ_.numberofSoftedges_ == 1){
      rx_start = 0;
      rx_end   = lx_ - Parameters_BHZ_.W_;
    } else {
      rx_start = Parameters_BHZ_.W_ - 1;
      rx_end   = lx_ - Parameters_BHZ_.W_;
    }
    //Sanity check:
    int expected_size = rx_end - rx_start + 1;
    assert(expected_size == ionic_size);

    for(int ind=0;ind<ionic_size;ind++){
        int rx = rx_start + ind;
        temp_den = Zero_Complex;
        for(int ry=0;ry<ly_;ry++){
            for(int orb=0;orb<N_orbs_;orb++){
                for(int spin=0;spin<N_spin_;spin++){
                    int r = Hamiltonian_BHZ_.makeIndex(rx, ry, orb, spin);

                    for(int n=0;n<H_size;n++){
                        temp_den += (1.0/(1.0*ly_))*calculateAmplitude(n,r,n,r)*fermifunction(Hamiltonian_BHZ_.evals_[n],mu);
                    }
                }
            }
        }
        avg_density[rx] = temp_den.real();
        if(ind>ionic_size-3){
            cout <<" avg_density["<<rx<<"]="<<avg_density[ind]<<endl;
        }
    }

    for(int rx=0;rx<ionic_size;rx++){
        double delta_n = avg_density[rx] - 2.0;
        new_ionic[rx] = Hamiltonian_BHZ_.ionic[rx] + ionic_relaxation * delta_n;
        if(rx>ionic_size-3){
            cout<<"delta_n="<<delta_n<<endl;
            cout<<"new_ionic["<<rx<<"]="<<new_ionic[rx]<<endl;
        }
    }
}

void Observables_BHZ::calculateDOS(){
    double value, w;

    string DOS_out="Density_of_states.txt";
    ofstream DOS_file_out(DOS_out.c_str());

    w=w_min;
    while(w<=w_max){
        value=0.0;
        DOS_file_out<< w<<"     ";

        for(int n=0;n<H_size;n++){
            value=value+(1.0/(4.0*lx_*ly_))*(one_by_PI_)*((eta)/((w-Hamiltonian_BHZ_.evals_[n])*(w-Hamiltonian_BHZ_.evals_[n])+(eta*eta)));
        }
        DOS_file_out<<value<<endl;
        w=w+dw;
    }
}

complex<double> Observables_BHZ::calculateAmplitude(int state1, int pos1, int state2, int pos2){
    complex<double> value;
    value = (conj(Hamiltonian_BHZ_.evecs_[state1][pos1]))*(Hamiltonian_BHZ_.evecs_[state2][pos2]);
    return value;
}

void Observables_BHZ::calculateLDOE(){

    Mat_4_Complex_doub ChargeDensity;
    ChargeDensity.resize(lx_);
    for(int rx=0;rx<lx_;rx++){
        ChargeDensity[rx].resize(ly_);
        for(int ry=0;ry<ly_;ry++){
            ChargeDensity[rx][ry].resize(N_orbs_);
            for(int orb=0;orb<N_orbs_;orb++){
                ChargeDensity[rx][ry][orb].resize(N_spin_);
                for(int spin=0;spin<N_spin_;spin++){
                    ChargeDensity[rx][ry][orb][spin]=Zero_Complex;
                }
            }
        }
    }

    complex<double> Tot_num_ples, local_charge_density;
    Tot_num_ples=Zero_Complex;

    int r;
//    #pragma omp parallel for collapse(4) private(r,local_charge_density) reduction(+:Tot_num_ples)
    for(int rx=0;rx<lx_;rx++){
        for(int ry=0;ry<ly_;ry++){
            for(int orb=0;orb<N_orbs_;orb++){
                for(int spin=0;spin<N_spin_;spin++){
                    r = Hamiltonian_BHZ_.makeIndex(rx, ry, orb, spin);

                    local_charge_density = Zero_Complex;
                    for(int n=0;n<H_size;n++){
                        local_charge_density += calculateAmplitude(n, r, n, r)*fermifunction(Hamiltonian_BHZ_.evals_[n],mu);
                    }
                    ChargeDensity[rx][ry][orb][spin] = local_charge_density;
                    Tot_num_ples += ChargeDensity[rx][ry][orb][spin];
                }
            }
        }
    }
    cout<<"Tot_num_ples_LDOE="<<Tot_num_ples.real()<<endl;

    double inv_ly=(1.0/(1.0*ly_));

    string file_Avg_LDOE="LDOE_plus_sz_along_x.txt";
    ofstream file_Avg_CD_out(file_Avg_LDOE.c_str());
    if (!file_Avg_CD_out.is_open()){
        cerr << "Error: Could not open file " << file_Avg_LDOE << endl;
        return;
    }

    file_Avg_CD_out <<"#rx       LDOE_Orb-0    LDOE_Orb-1    Sz_Orb-0    Sz_Orb-1"<< endl;

    complex<double> CD_s_up,CD_s_dn,CD_p_up,CD_p_dn;
    for(int rx=0;rx<lx_;rx++){
        CD_s_up = Zero_Complex;
        CD_s_dn = Zero_Complex;
        CD_p_up = Zero_Complex;
        CD_p_dn = Zero_Complex;

        for(int ry=0;ry<ly_;ry++){
            CD_s_up += inv_ly*ChargeDensity[rx][ry][0][0];
            CD_s_dn += inv_ly*ChargeDensity[rx][ry][0][1];
            CD_p_up += inv_ly*ChargeDensity[rx][ry][1][0];
            CD_p_dn += inv_ly*ChargeDensity[rx][ry][1][1];
        }
        file_Avg_CD_out<<rx<<"  "<<CD_s_up.real()+CD_s_dn.real()<<"    "<<CD_p_up.real()+CD_p_dn.real()<<"      "<<
                                    (CD_s_up.real()-CD_s_dn.real())*0.5<<"       "<<0.5*(CD_p_up.real()-CD_p_dn.real())<<endl;
    }

    if(!Parameters_BHZ_.canonical_){
    string file_Bulk_Edge_LDOE="Bulk_edge_particles.txt";
    ofstream file_BE_out(file_Bulk_Edge_LDOE.c_str());
    if (!file_BE_out.is_open()) {
        cerr << "Error: Could not open file " << file_Bulk_Edge_LDOE << endl;
        return;
    }

    complex<double> Tot_edge_ples_orb_s,Tot_edge_ples_orb_p;
    Tot_edge_ples_orb_s=Zero_Complex;
    Tot_edge_ples_orb_p=Zero_Complex;

    complex<double> Tot_edge_ples, Tot_bulk_ples;
    Tot_edge_ples=Zero_Complex;
    Tot_bulk_ples=Zero_Complex;

    complex<double> Hard_edge_ples, SoftHard_edge_ples, Bulk_ples;
    Hard_edge_ples=Zero_Complex;
    SoftHard_edge_ples=Zero_Complex;
    Bulk_ples = Zero_Complex;

    for(int rx=0;rx<lx_;rx++){
        for(int ry=0;ry<ly_;ry++){
            if(Parameters_BHZ_.numberofSoftedges_==1){
                if(rx>lx_-Hamiltonian_BHZ_.wx_){
                    Tot_edge_ples_orb_s += ChargeDensity[rx][ry][0][0] + ChargeDensity[rx][ry][0][1];
                    Tot_edge_ples_orb_p += ChargeDensity[rx][ry][1][0] + ChargeDensity[rx][ry][1][1];
                    Tot_edge_ples=Tot_edge_ples_orb_s + Tot_edge_ples_orb_p;
                }
                else{
                    Tot_bulk_ples+=ChargeDensity[rx][ry][0][0] + ChargeDensity[rx][ry][0][1] + ChargeDensity[rx][ry][1][0] + ChargeDensity[rx][ry][1][1];
                }
                if(rx>=0 && rx<=2){
                    Hard_edge_ples += ChargeDensity[rx][ry][0][0] + ChargeDensity[rx][ry][0][1] + ChargeDensity[rx][ry][1][0] + ChargeDensity[rx][ry][1][1];
                }
                else if(rx<=lx_-Hamiltonian_BHZ_.wx_ && rx>=lx_-Hamiltonian_BHZ_.wx_-2){
                    SoftHard_edge_ples += ChargeDensity[rx][ry][0][0] + ChargeDensity[rx][ry][0][1] + ChargeDensity[rx][ry][1][0] + ChargeDensity[rx][ry][1][1];
                }
                else if(rx>2 && rx<lx_-Hamiltonian_BHZ_.wx_-2){
                    Bulk_ples += ChargeDensity[rx][ry][0][0] + ChargeDensity[rx][ry][0][1] + ChargeDensity[rx][ry][1][0] + ChargeDensity[rx][ry][1][1];
                }
                else{

                }
            }
            else{
                if(rx<Hamiltonian_BHZ_.wx_-1 || rx>lx_-Hamiltonian_BHZ_.wx_){
                    Tot_edge_ples_orb_s += ChargeDensity[rx][ry][0][0] + ChargeDensity[rx][ry][0][1];
                    Tot_edge_ples_orb_p += ChargeDensity[rx][ry][1][0] + ChargeDensity[rx][ry][1][1];
                    Tot_edge_ples=Tot_edge_ples_orb_s + Tot_edge_ples_orb_p;
                }
                else{
                    Tot_bulk_ples+=ChargeDensity[rx][ry][0][0] + ChargeDensity[rx][ry][0][1] + ChargeDensity[rx][ry][1][0] + ChargeDensity[rx][ry][1][1];
                }
            }
        }
    }

    file_BE_out<<"Edge_ples_orb_s= "<<Tot_edge_ples_orb_s.real()<<endl;
    file_BE_out<<"Edge_ples_orb_p= "<<Tot_edge_ples_orb_p.real()<<endl;
    file_BE_out<<"Total_edge_ples= "<<Tot_edge_ples.real()<<endl;
    file_BE_out<<"Total_bulk_ples= "<<Tot_bulk_ples.real()<<endl;
    file_BE_out<<"Total_num_ples= "<<Tot_num_ples.real()<<endl;
    file_BE_out<<"Hard_edge_ples= "<<Hard_edge_ples.real()<<endl;
    file_BE_out<<"HardSoft_edge_ples= "<<SoftHard_edge_ples.real()<<endl;
    file_BE_out<<"Bulk_ples= "<<Bulk_ples.real()<<endl;
    }

}

void Observables_BHZ::calculateSpinCurrents(){
    Mat_4_Complex_doub SpinCurrent_UP, SpinCurrent_DN;

    SpinCurrent_UP.resize(N_cells_);
    SpinCurrent_DN.resize(N_cells_);
    for(int cell1=0;cell1<N_cells_;cell1++){
        SpinCurrent_UP[cell1].resize(N_cells_);
        SpinCurrent_DN[cell1].resize(N_cells_);

        for(int cell2=0;cell2<N_cells_;cell2++){
            SpinCurrent_UP[cell1][cell2].resize(N_orbs_);
            SpinCurrent_DN[cell1][cell2].resize(N_orbs_);

            for(int orb1=0;orb1<N_orbs_;orb1++){
                SpinCurrent_UP[cell1][cell2][orb1].resize(N_orbs_);
                SpinCurrent_DN[cell1][cell2][orb1].resize(N_orbs_);

                for(int orb2=0;orb2<N_orbs_;orb2++){
                    SpinCurrent_UP[cell1][cell2][orb1][orb2]=Zero_Complex;
                    SpinCurrent_DN[cell1][cell2][orb1][orb2]=Zero_Complex;
                }
            }
        }
    }

    int r1, r2, cell1, cell2;
//    complex<double> local_spinup_curr, local_spindn_curr;
//    #pragma omp parallel for collapse(6) private(r1, r2, cell1, cell2)
    for(int r1x=0;r1x<lx_;r1x++){
        for(int r1y=0;r1y<ly_;r1y++){
            cell1 = r1y + ly_*r1x;

            for(int orb1=0;orb1<N_orbs_;orb1++){
                r1 = Hamiltonian_BHZ_.makeIndex(r1x, r1y, orb1, 0);

                for(int r2x=0;r2x<lx_;r2x++){
                    for(int r2y=0;r2y<ly_;r2y++){
                        cell2 = r2y + ly_*r2x;

                        for(int orb2=0;orb2<N_orbs_;orb2++){
                            r2 = Hamiltonian_BHZ_.makeIndex(r2x, r2y, orb2, 0);

                            if(cell1 < cell2){
                                if(Hamiltonian_BHZ_.C_mat[r1][r2].real()!=0 || Hamiltonian_BHZ_.C_mat[r1][r2].imag()!=0){
                                    for(int n=0;n<H_size;n++){
                                        SpinCurrent_UP[cell1][cell2][orb1][orb2] += 
                                        0.5*Iota_Complex*( ( Hamiltonian_BHZ_.C_mat[r2][r1]*calculateAmplitude(n, r2, n, r1) ) - 
                                        ( Hamiltonian_BHZ_.C_mat[r1][r2]*calculateAmplitude(n, r1, n, r2) ) )*(fermifunction(Hamiltonian_BHZ_.evals_[n],mu) );
                                    }
                                }
                                if(Hamiltonian_BHZ_.C_mat[r1+MHS][r2+MHS].real()!=0 || Hamiltonian_BHZ_.C_mat[r1+MHS][r2+MHS].imag()!=0){
                                    for(int n=0;n<H_size;n++){
                                        SpinCurrent_DN[cell1][cell2][orb1][orb2] 
                                        += 0.5*Iota_Complex*( ( Hamiltonian_BHZ_.C_mat[r2+MHS][r1+MHS]*calculateAmplitude(n, r2+MHS, n, r1+MHS) )- 
                                        ( Hamiltonian_BHZ_.C_mat[r1+MHS][r2+MHS]*calculateAmplitude(n, r1+MHS, n, r2+MHS) ) )*(fermifunction(Hamiltonian_BHZ_.evals_[n],mu) );
                                    }
                                }
                            //    cout<<cell1<<"  "<<cell2<<" "<<orb1<<"  "<<orb2<<"  "<<SpinCurrent_UP[cell1][cell2][orb1][orb2].real()<<"   "
                            //    <<SpinCurrent_DN[cell1][cell2][orb1][orb2].real()<<endl;
                            }
                        }
                    }
                }
            }
        }
    }

    string file_spin_current="Total_spin_current_at_each_cell_link.txt";
    ofstream file_spin_current_out(file_spin_current.c_str());
    if (!file_spin_current_out.is_open()) {
        cerr << "Error: Could not open file " << file_spin_current << endl;
        return;
    }

    for(int r1x=0;r1x<lx_;r1x++){
        for(int r1y=0;r1y<ly_;r1y++){
            cell1 = r1y + ly_*r1x;

            for(int r2x=0;r2x<lx_;r2x++){
                for(int r2y=0;r2y<ly_;r2y++){
                    cell2 = r2y + ly_*r2x;

                    if(cell1<cell2){
                        if(SpinCurrent_UP[cell1][cell2][0][0].real()!=0 || SpinCurrent_UP[cell1][cell2][1][1].real()!=0){
                            file_spin_current_out<<cell1<<" "<<cell2<<" "
                            <<SpinCurrent_UP[cell1][cell2][0][0].real()-SpinCurrent_DN[cell1][cell2][0][0].real()<<"  "
                            <<SpinCurrent_UP[cell1][cell2][1][1].real()-SpinCurrent_DN[cell1][cell2][1][1].real()<<endl;
                        }
                    }
                }
            }
        }
    }

    string file_avg_spin_current="Avg_spin_current_along_ry_vs_rx.txt";
    ofstream file_avg_spin_current_out(file_avg_spin_current.c_str());
    if (!file_avg_spin_current_out.is_open()) {
        cerr << "Error: Could not open file " << file_avg_spin_current << endl;
        return;
    }

    file_avg_spin_current_out<<"#rx     #(s->s),up      #(p->p),up      #(s->s),dn      #(p->p),dn"<<endl;
    
    complex<double> SC_ss_up,SC_pp_up,SC_ss_dn,SC_pp_dn;
    for(int r1x=0;r1x<lx_;r1x++){
        SC_ss_up=Zero_Complex;  SC_pp_up=Zero_Complex;
        SC_ss_dn=Zero_Complex;  SC_pp_dn=Zero_Complex;

        for(int r1y=0;r1y<ly_;r1y++){
            cell1 = r1y + ly_*r1x;

            for(int r2y=0;r2y<ly_;r2y++){
                cell2 = r2y + ly_*r1x;

                if(cell1<cell2){
                    SC_ss_up += SpinCurrent_UP[cell1][cell2][0][0];
                    SC_pp_up += SpinCurrent_UP[cell1][cell2][1][1];

                    SC_ss_dn += SpinCurrent_DN[cell1][cell2][0][0];
                    SC_pp_dn += SpinCurrent_DN[cell1][cell2][1][1];
                }
            }
        }
        file_avg_spin_current_out<<r1x<<"   "<<SC_ss_up.real()<<"     "<<SC_pp_up.real()<<"     "<<SC_ss_dn.real()<<"   "<<SC_pp_dn.real()<<endl;
    }
    
}

void Observables_BHZ::calculateLDOS(){
    string file_ldos("LDOS_along_x.txt");
    ofstream file_ldos_out(file_ldos.c_str());

    double w=0.0;
    int r0_up, r0_dn, r1_up, r1_dn;
    complex<double> ldos_0_up, ldos_0_dn, ldos_1_up, ldos_1_dn;

    for(int rx=0;rx<lx_;rx++){
        for(int om=0;om<w_size;om++){
            w=w_min+om*dw;
            ldos_0_up=Zero_Complex;
            ldos_0_dn=Zero_Complex;
            ldos_1_up=Zero_Complex;
            ldos_1_dn=Zero_Complex;

            for(int ry=0;ry<ly_;ry++){
                r0_up = Hamiltonian_BHZ_.makeIndex(rx, ry, 0, 0);
                r0_dn = Hamiltonian_BHZ_.makeIndex(rx, ry, 0, 1);
                r1_up = Hamiltonian_BHZ_.makeIndex(rx, ry, 1, 0);
                r1_dn = Hamiltonian_BHZ_.makeIndex(rx, ry, 1, 1);

                for(int n=0;n<H_size;n++){
                    ldos_0_up += (one_by_PI_)*(calculateAmplitude(n, r0_up, n, r0_up))*
                        ( (eta)/((w-Hamiltonian_BHZ_.evals_[n])*(w-Hamiltonian_BHZ_.evals_[n])+(eta*eta)) );
                    ldos_0_dn += (one_by_PI_)*(calculateAmplitude(n, r0_dn, n, r0_dn))*
                        ( (eta)/((w-Hamiltonian_BHZ_.evals_[n])*(w-Hamiltonian_BHZ_.evals_[n])+(eta*eta)) );
                    ldos_1_up += (one_by_PI_)*(calculateAmplitude(n, r1_up, n, r1_up))*
                        ( (eta)/((w-Hamiltonian_BHZ_.evals_[n])*(w-Hamiltonian_BHZ_.evals_[n])+(eta*eta)) );
                    ldos_1_dn += (one_by_PI_)*(calculateAmplitude(n, r1_dn, n, r1_dn))*
                        ( (eta)/((w-Hamiltonian_BHZ_.evals_[n])*(w-Hamiltonian_BHZ_.evals_[n])+(eta*eta)) );
                }

            }
            file_ldos_out<<rx<<"    "<<w<<"     "<<ldos_0_up.real()<<"   "<<ldos_0_dn.real()<<"     "<<ldos_1_up.real()<<"  "<<ldos_1_dn.real()<<endl;
        }
        file_ldos_out<<endl;
    }

}

void Observables_BHZ::calculateAkxw(){
    
    Mat_3_Complex_doub B_mat;
    B_mat.resize(lx_);
    for(int r1x=0;r1x<lx_;r1x++){
        B_mat[r1x].resize(lx_);
        for(int r2x=0;r2x<lx_;r2x++){
            B_mat[r1x][r2x].resize(w_size);
            for(int om=0;om<w_size;om++){
                B_mat[r1x][r2x][om] = Zero_Complex;
            }
        }
    }
    
    int r1,r2;
    double w;
    w=0.0;
    for(int r1x=0;r1x<lx_;r1x++){
        for(int r2x=0;r2x<lx_;r2x++){
            for(int om=0;om<w_size;om++){
                //w=w_min+om*dw;
                w=w_min + mu + om*dw;

                for(int ry=0;ry<ly_;ry++){
                    for(int orb1=0;orb1<N_orbs_;orb1++){
                        for(int orb2=0;orb2<N_orbs_;orb2++){
                            for(int spin=0;spin<N_spin_;spin++){

                                r1 = Hamiltonian_BHZ_.makeIndex(r1x, ry, orb1, spin);
                                r2 = Hamiltonian_BHZ_.makeIndex(r2x, ry, orb2, spin);

                                for(int n=0;n<H_size;n++){
                                    B_mat[r1x][r2x][om] += (one_by_PI_)*(calculateAmplitude(n, r1, n, r2))*
                                        ( (eta)/((w-Hamiltonian_BHZ_.evals_[n])*(w-Hamiltonian_BHZ_.evals_[n])+(eta*eta)) );
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    int kx_ind,ky_ind,k_point;
    double kx,ky;
    complex<double> Akxw_;

    int kx_min,kx_max;
    double momentum_step;
    if(Parameters_BHZ_.PBC_X==true){
        kx_min=-lx_/2;       kx_max=lx_/2;
        momentum_step = ((2.0*PI)/(1.0*lx_));
    }
    else{
        kx_min=1;       kx_max=lx_;
        momentum_step = ((1.0*PI)/(1.0*lx_+1.0));
    }
    
    string file_Akxw("Akxw.txt");
    ofstream file_Akxw_out(file_Akxw.c_str());

    for(kx_ind=kx_min;kx_ind<=kx_max;kx_ind++){
        kx=kx_ind*momentum_step;

        for(int om=0;om<w_size;om++){
            //w=w_min+om*dw;
            w = w_min + mu + om*dw;
            Akxw_ = Zero_Complex;
            for(int r1x=0;r1x<lx_;r1x++){
                for(int r2x=0;r2x<lx_;r2x++){

                    if(Parameters_BHZ_.PBC_X==true){
                        Akxw_ += (1.0/(1.0*lx_))*( exp(Iota_Complex*kx*(1.0*(r1x-r2x)))*B_mat[r1x][r2x][om] );
                    }
                    else{
                        Akxw_ += (2.0/((lx_+1)*1.0))*( sin(kx*(1.0*r1x+1.0))*sin(kx*(1.0*r2x+1.0))*B_mat[r1x][r2x][om] );
                    }

                }
            }
            file_Akxw_out<<kx<<"    "<<w<<"         "<<Akxw_.real()<<endl;
        }
        file_Akxw_out<<endl;
    }

}

void Observables_BHZ::calculateAkyw(){

    Mat_4_Complex_doub B_s_mat,B_p_mat,B_mat;
    B_s_mat.resize(ly_);    B_p_mat.resize(ly_);
    B_mat.resize(ly_);

    for(int r1y=0;r1y<ly_;r1y++){
        B_s_mat[r1y].resize(ly_);    B_p_mat[r1y].resize(ly_);
        B_mat[r1y].resize(ly_);

        for(int r2y=0;r2y<ly_;r2y++){
            B_s_mat[r1y][r2y].resize(N_spin_);      B_p_mat[r1y][r2y].resize(N_spin_);
            B_mat[r1y][r2y].resize(N_spin_);

            for(int spin=0;spin<N_spin_;spin++){
                B_s_mat[r1y][r2y][spin].resize(w_size);      B_p_mat[r1y][r2y][spin].resize(w_size);
                B_mat[r1y][r2y][spin].resize(w_size);

                for(int om=0;om<w_size;om++){
                    B_s_mat[r1y][r2y][spin][om] = Zero_Complex;     B_p_mat[r1y][r2y][spin][om] = Zero_Complex;
                    B_mat[r1y][r2y][spin][om] = Zero_Complex;
                }
            }
        }
    }

    int r1,r2;
    double w=0.0;

    for(int r1y=0;r1y<ly_;r1y++){
        for(int r2y=0;r2y<ly_;r2y++){
            for(int spin=0;spin<N_spin_;spin++){
                for(int om=0;om<w_size;om++){
                    //w=w_min+om*dw;
                    w=w_min + mu + om*dw;

                    for(int orb=0;orb<N_orbs_;orb++){
                        for(int rx=0;rx<lx_;rx++){
                            
                            r1 = Hamiltonian_BHZ_.makeIndex(rx, r1y, orb, spin);
                            r2 = Hamiltonian_BHZ_.makeIndex(rx, r2y, orb, spin);

                            for(int n=0;n<H_size;n++){
                                B_mat[r1y][r2y][spin][om] += (one_by_PI_)*(calculateAmplitude(n, r1, n, r2))*
                                                            ( (eta)/((w-Hamiltonian_BHZ_.evals_[n])*(w-Hamiltonian_BHZ_.evals_[n])+(eta*eta)) );
                            }

                            if(orb==0){
                                for(int n=0;n<H_size;n++){
                                    B_s_mat[r1y][r2y][spin][om] += (one_by_PI_)*(calculateAmplitude(n, r1, n, r2))*
                                        ( (eta)/((w-Hamiltonian_BHZ_.evals_[n])*(w-Hamiltonian_BHZ_.evals_[n])+(eta*eta)) );
                                }
                            }

                            if(orb==1){
                                for(int n=0;n<H_size;n++){
                                    B_p_mat[r1y][r2y][spin][om] += (one_by_PI_)*(calculateAmplitude(n, r1, n, r2))*
                                        ( (eta)/((w-Hamiltonian_BHZ_.evals_[n])*(w-Hamiltonian_BHZ_.evals_[n])+(eta*eta)) );
                                }
                            }

                        }
                    }
                }
            }
        }
    }

    string file_Akyw("Akyw.txt");
    ofstream file_Akyw_out(file_Akyw.c_str());
    file_Akyw_out<<"#ky"<<" "<<"omega"<<"   "<<"Akyw_up"<<"     "<<"Akyw_dn"<<"     "<<"Akyw_tot"<<endl;

    string file_Akyw_s_spin_resolved("Akyw_s_spin_resolved.txt");
    ofstream file_Akyw_s_sr_out(file_Akyw_s_spin_resolved.c_str());
    file_Akyw_s_sr_out<<"#ky"<<" "<<"omega"<<"   "<<"Akyw_s_up"<<"  "<<"Akyw_s_dn"<<"  "<<"Akyw_s_tot"<<endl;

    string file_Akyw_p_spin_resolved("Akyw_p_spin_resolved.txt");
    ofstream file_Akyw_p_sr_out(file_Akyw_p_spin_resolved.c_str());
    file_Akyw_p_sr_out<<"#ky"<<" "<<"omega"<<"   "<<"Akyw_p_up"<<"  "<<"Akyw_p_dn"<<"  "<<"Akyw_p_tot"<<endl;

    complex<double> Akyw_s,Akyw_s_up,Akyw_s_dn;
    complex<double> Akyw_p,Akyw_p_up,Akyw_p_dn;
    complex<double> Akyw_,Akyw_up,Akyw_dn;

    int ky_ind,ky_min,ky_max;
    double momentum_step, ky;
    if(Parameters_BHZ_.PBC_Y){
        ky_min=-ly_/2;       ky_max=ly_/2;
        momentum_step = ((2.0*PI)/(1.0*ly_));
    }
    else{
        ky_min=1;       ky_max=ly_;
        momentum_step = ((1.0*PI)/(1.0*ly_+1.0));
    }

    for(ky_ind=ky_min;ky_ind<=ky_max;ky_ind++){
        ky=ky_ind*momentum_step;

        for(int om=0;om<w_size;om++){
            //w=w_min+om*dw;
            w=w_min + mu + om*dw;

            Akyw_s=Zero_Complex;Akyw_s_up=Zero_Complex;Akyw_s_dn=Zero_Complex;
            Akyw_p=Zero_Complex;Akyw_p_up=Zero_Complex;Akyw_p_dn=Zero_Complex;
            Akyw_=Zero_Complex;Akyw_up=Zero_Complex;Akyw_dn=Zero_Complex;

            for(int r1y=0;r1y<ly_;r1y++){
                for(int r2y=0;r2y<ly_;r2y++){

                        if(Parameters_BHZ_.PBC_Y){
                            Akyw_s_up += (1.0/(1.0*ly_*lx_))*( exp(Iota_Complex*ky*(1.0*(r1y-r2y)))*B_s_mat[r1y][r2y][0][om] );
                            Akyw_s_dn += (1.0/(1.0*ly_*lx_))*( exp(Iota_Complex*ky*(1.0*(r1y-r2y)))*B_s_mat[r1y][r2y][1][om] );
                            Akyw_s += (1.0/(2.0*ly_*lx_))*( exp(Iota_Complex*ky*(1.0*(r1y-r2y)))* (B_s_mat[r1y][r2y][0][om] + B_s_mat[r1y][r2y][1][om]) );

                            Akyw_p_up += (1.0/(1.0*ly_*lx_))*( exp(Iota_Complex*ky*(1.0*(r1y-r2y)))*B_p_mat[r1y][r2y][0][om] );
                            Akyw_p_dn += (1.0/(1.0*ly_*lx_))*( exp(Iota_Complex*ky*(1.0*(r1y-r2y)))*B_p_mat[r1y][r2y][1][om] );
                            Akyw_p += (1.0/(2.0*ly_*lx_))*( exp(Iota_Complex*ky*(1.0*(r1y-r2y)))* (B_p_mat[r1y][r2y][0][om] + B_p_mat[r1y][r2y][1][om]) );

                            Akyw_up += (1.0/(1.0*ly_*lx_))*( exp(Iota_Complex*ky*(1.0*(r1y-r2y)))*B_mat[r1y][r2y][0][om] );
                            Akyw_dn += (1.0/(1.0*ly_*lx_))*( exp(Iota_Complex*ky*(1.0*(r1y-r2y)))*B_mat[r1y][r2y][1][om] );
                            Akyw_ += (1.0/(2.0*ly_*lx_))*( exp(Iota_Complex*ky*(1.0*(r1y-r2y)))* (B_mat[r1y][r2y][0][om] + B_mat[r1y][r2y][1][om]) );
                        }
                        else{
                            Akyw_ += (2.0/((ly_+1)*1.0))*( sin(ky*(1.0*r1y+1.0))*sin(ky*(1.0*r2y+1.0))*(B_mat[r1y][r2y][0][om] + B_mat[r1y][r2y][1][om]) );
                        }

                }
            }
            file_Akyw_out<<ky<<"    "<<w<<"         "<<Akyw_up.real()<<"    "<<Akyw_dn.real()<<"    "<<Akyw_.real()<<endl;
            file_Akyw_p_sr_out<<ky<<"    "<<w<<"         "<<Akyw_p_up.real()<<"     "<<Akyw_p_dn.real()<<"      "<<Akyw_p.real()<<endl;
            file_Akyw_s_sr_out<<ky<<"    "<<w<<"         "<<Akyw_s_up.real()<<"     "<<Akyw_s_dn.real()<<"      "<<Akyw_s.real()<<endl;
        }
        file_Akyw_s_sr_out<<endl;
        file_Akyw_p_sr_out<<endl;
        file_Akyw_out<<endl;
    }
}

void Observables_BHZ::calculateWaveFunctions(){

    double E_min=-0.2, E_max=0.2, eval;
    eval = 0.0;
    int state_,  rs_up,rs_dn, rp_up,rp_dn;
    complex<double>val_s_up,val_s_dn, val_p_up,val_p_dn;
    
    string head = "wave_function_for_state_";
    string base = ".txt";

    for(int np=0;np<Hamiltonian_BHZ_.evals_.size();np++){
        eval = Hamiltonian_BHZ_.evals_[np];

        if(eval > E_min && eval < E_max){
            state_=np;
            val_s_up = Zero_Complex;    val_s_dn = Zero_Complex;
            val_p_up = Zero_Complex;    val_p_dn = Zero_Complex;

            string file_wave_fn_(head + to_string(state_) + base);
            ofstream file_out(file_wave_fn_.c_str());
            if (!file_out.is_open()) {
                cerr << "Error: Could not open file " << file_wave_fn_ << endl;
                continue;  // Skipping to next state when a file opening fails
            }

            file_out << "#r=ry+ly*rx    phi(r,0,0)    phi(r,0,1)  phi(r,1,0)  phi(r,1,1)" << endl;

            for(int rx=0;rx<lx_;rx++){
                for(int ry=0;ry<ly_;ry++){
                    rs_up = Hamiltonian_BHZ_.makeIndex(rx, ry, 0, 0);
                    rs_dn = Hamiltonian_BHZ_.makeIndex(rx, ry, 0, 1);
                    rp_up = Hamiltonian_BHZ_.makeIndex(rx, ry, 1, 0);
                    rp_dn = Hamiltonian_BHZ_.makeIndex(rx, ry, 1, 1);
                                
                    val_s_up = calculateAmplitude(state_, rs_up, state_, rs_up);
                    val_s_dn = calculateAmplitude(state_, rs_dn, state_, rs_dn);
                    val_p_up = calculateAmplitude(state_, rp_up, state_, rp_up);
                    val_p_dn = calculateAmplitude(state_, rp_dn, state_, rp_dn);
                    
                    file_out<<ry+ly_*rx<<"  "<<val_s_up.real()<<"   "<<val_s_dn.real()<<"   "<<val_p_up.real()<<
                    "   "<<val_p_dn.real()<<endl;
                     
                }
            }

        }
    }

}

void Observables_BHZ::calculateSpinSpinCorr(){

    Mat_6_Complex_doub SS_occ,SS_unocc;
    SS_occ.resize(N_cells_);
    SS_unocc.resize(N_cells_);
    for(int cell1=0;cell1<N_cells_;cell1++){
        SS_occ[cell1].resize(N_cells_);
        SS_unocc[cell1].resize(N_cells_);
        for(int cell2=0;cell2<N_cells_;cell2++){
            SS_occ[cell1][cell2].resize(N_orbs_);
            SS_unocc[cell1][cell2].resize(N_orbs_);
            for(int orb1=0;orb1<N_orbs_;orb1++){
                SS_occ[cell1][cell2][orb1].resize(N_orbs_);
                SS_unocc[cell1][cell2][orb1].resize(N_orbs_);
                for(int orb2=0;orb2<N_orbs_;orb2++){
                    SS_occ[cell1][cell2][orb1][orb2].resize(N_spin_);
                    SS_unocc[cell1][cell2][orb1][orb2].resize(N_spin_);
                    for(int spin1=0;spin1<N_spin_;spin1++){
                        SS_occ[cell1][cell2][orb1][orb2][spin1].resize(N_spin_);
                        SS_unocc[cell1][cell2][orb1][orb2][spin1].resize(N_spin_);
                    }
                }
            }
        }
    }

    int r1,r2;
    int r1x,r2x,r1y,r2y;
    for(int cell1=0;cell1<N_cells_;cell1++){
        r1x = cell1/ly_;
        r1y = cell1%ly_;
        for(int cell2=0;cell2<N_cells_;cell2++){
            r2x = cell2/ly_;
            r2y = cell2%ly_;
            for(int orb1=0;orb1<N_orbs_;orb1++){
                for(int orb2=0;orb2<N_orbs_;orb2++){
                    for(int spin1=0;spin1<N_spin_;spin1++){
                        for(int spin2=0;spin2<N_spin_;spin2++){

                            r1 = Hamiltonian_BHZ_.makeIndex(r1x, r1y, orb1, spin1);
                            r2 = Hamiltonian_BHZ_.makeIndex(r2x, r2y, orb2, spin2);

                            SS_occ[cell1][cell2][orb1][orb2][spin1][spin2] = Zero_Complex;
                            SS_unocc[cell1][cell2][orb1][orb2][spin1][spin2] = Zero_Complex;
                            for(int n=0;n<H_size;n++){
                                SS_occ[cell1][cell2][orb1][orb2][spin1][spin2] += calculateAmplitude(n, r1, n, r2)*fermifunction(Hamiltonian_BHZ_.evals_[n], mu);
                                SS_unocc[cell1][cell2][orb1][orb2][spin1][spin2] += calculateAmplitude(n, r1, n, r2)*(1.0-fermifunction(Hamiltonian_BHZ_.evals_[n], mu));
                            }

                        }
                    }
                }
            }
        }
    }

    Mat_4_Complex_doub SS_orb_resolved;
    SS_orb_resolved.resize(N_cells_);
    for(int cell1=0;cell1<N_cells_;cell1++){
        SS_orb_resolved[cell1].resize(N_cells_);
        for(int cell2=0;cell2<N_cells_;cell2++){
            SS_orb_resolved[cell1][cell2].resize(N_orbs_);
            for(int orb1=0;orb1<N_orbs_;orb1++){
                SS_orb_resolved[cell1][cell2][orb1].resize(N_orbs_);
                for(int orb2=0;orb2<N_orbs_;orb2++){
                    SS_orb_resolved[cell1][cell2][orb1][orb2]=Zero_Complex;
                }
            }
        }
    }

    for(int cell1=0;cell1<N_cells_;cell1++){
        for(int cell2=0;cell2<N_cells_;cell2++){
            for(int orb1=0;orb1<N_orbs_;orb1++){
                for(int orb2=0;orb2<N_orbs_;orb2++){

                    //0.5*(<S+S- + S-S+) contribution:
                    SS_orb_resolved[cell1][cell2][orb1][orb2] +=  
                    //Direct term:--->
                                0.5*(SS_occ[cell1][cell1][orb1][orb1][0][1]*SS_occ[cell2][cell2][orb2][orb2][1][0] + 
                                SS_occ[cell1][cell1][orb1][orb1][1][0]*SS_occ[cell2][cell2][orb2][orb2][0][1]) +
                    //Exchange term:--->
                                0.5*(SS_occ[cell1][cell2][orb1][orb2][0][0]*SS_unocc[cell2][cell1][orb2][orb1][1][1] + 
                                SS_occ[cell1][cell2][orb1][orb2][1][1]*SS_unocc[cell2][cell1][orb2][orb1][0][0]);

                    //SzSz contribution:
                    SS_orb_resolved[cell1][cell2][orb1][orb2] +=
                    //Direct term:--->
                                0.25*(SS_occ[cell1][cell1][orb1][orb1][0][0]*SS_occ[cell2][cell2][orb2][orb2][0][0] +
                                SS_occ[cell1][cell1][orb1][orb1][1][1]*SS_occ[cell2][cell2][orb2][orb2][1][1] -
                                SS_occ[cell1][cell1][orb1][orb1][0][0]*SS_occ[cell2][cell2][orb2][orb2][1][1] -
                                SS_occ[cell1][cell1][orb1][orb1][1][1]*SS_occ[cell2][cell2][orb2][orb2][0][0] ) + 
                    //Exchange term:--->
                                0.25*(SS_occ[cell1][cell2][orb1][orb2][0][0]*SS_unocc[cell2][cell1][orb2][orb1][0][0] +
                                SS_occ[cell1][cell2][orb1][orb2][1][1]*SS_unocc[cell2][cell1][orb2][orb1][1][1] - 
                                SS_occ[cell1][cell2][orb1][orb2][0][1]*SS_unocc[cell2][cell1][orb2][orb1][1][0] -
                                SS_occ[cell1][cell2][orb1][orb2][1][0]*SS_unocc[cell2][cell1][orb2][orb1][0][1]);
                }
            }
        }
    }

    string spin_spin_corr = "SiSj_corr.txt";
    ofstream spin_spin_corr_out(spin_spin_corr.c_str());
    spin_spin_corr_out<<"#cell1     r1x      r1y      cell2      r2x     r2y     S_{cell1}.S_{cell2}"<<endl;

    Mat_2_Complex_doub SS_mat;
    SS_mat.resize(N_cells_);
    for(int cell1=0;cell1<N_cells_;cell1++){
        SS_mat[cell1].resize(N_cells_);
        for(int cell2=0;cell2<N_cells_;cell2++){
            SS_mat[cell1][cell2]=Zero_Complex;
        }
    }

    for(int cell1=0;cell1<N_cells_;cell1++){
        r1x = cell1/ly_;
        r1y = cell1%ly_;
        for(int cell2=0;cell2<N_cells_;cell2++){
            r2x = cell2/ly_;
            r2y = cell2%ly_;
            for(int orb1=0;orb1<N_orbs_;orb1++){
                for(int orb2=0;orb2<N_orbs_;orb2++){
                    SS_mat[cell1][cell2] += SS_orb_resolved[cell1][cell2][orb1][orb2];
                }
            }
            spin_spin_corr_out<<cell1<<"        "<<r1x<<"      "<<r1y<<"        "<<cell2<<"       "<<r2x<<"       "<<r2y<<"        "
                                    <<SS_mat[cell1][cell2].real()<<"        "<<SS_mat[cell1][cell2].imag()<<endl;
        }
    }


    string local_moments = "Magnetic_moments_at_each_cell.txt";
    ofstream local_moments_out(local_moments.c_str());
    local_moments_out<<"#cell      rx      ry      S^2_{cell}"<<endl;

    int rx, ry;
    for(int cell=0;cell<N_cells_;cell++){
        rx = cell/ly_;
        ry = cell%ly_;
        local_moments_out<<cell<<"      "<<rx<<"        "<<ry<<"        "<<SS_mat[cell][cell].real()<<"     "<<SS_mat[cell][cell].imag()<<endl;
    }

    string avg_local_moments = "Avg_local_moments_vs_rx.txt";
    ofstream avg_moments_out(avg_local_moments.c_str());
    avg_moments_out<<"#rx       avg_S^2"<<endl;

    complex<double> avg_moment;
    int cell;
    for(int rx=0;rx<lx_;rx++){
        avg_moment=Zero_Complex;
        for(int ry=0;ry<ly_;ry++){
            cell = ry + ly_*rx;
            avg_moment += ((1.0)/(1.0*ly_))*SS_mat[cell][cell];
        }
        avg_moments_out<<rx<<"      "<<avg_moment.real()<<"     "<<avg_moment.imag()<<endl;
    }

}

#endif
