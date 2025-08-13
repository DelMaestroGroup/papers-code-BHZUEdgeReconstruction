#include <cassert>
#include "tensors.h"
#include "Parameters_BHZ.h"
extern "C" {
//    #include <lapacke.h>
//    #include <cblas.h>
    void zheev_(char *, char *, int *, std::complex<double> *, int *, double *, std::complex<double> *, int *, double *, int *);
}

#ifndef Hamiltonian_BHZ_class
#define Hamiltonian_BHZ_class
#define PI acos(-1.0)

class Hamiltonian_BHZ{

    public:
    //Constructor:
    Hamiltonian_BHZ(Parameters_BHZ &Parameters_BHZ__): Parameters_BHZ_(Parameters_BHZ__){
        Initialize();
    }
    void Initialize();
    int bar_function(int val);
    int makeIndex(int rx, int ry, int orb, int spin);

    void addOnsiteTerms(int r, int orb, int spin);
    void addHoppingXTerm(int r, int rpx, int orb);
    void addHoppingYTerm(int r, int rpy, int orb);
    void addOrbMixingXTerm(int r, int rpx, int orb, int spin);
    void addOrbMixingYTerm(int r, int rpy, int orb, int spin);
    void addOnsiteHartreeTerm(int rx, int ry, int orb, int spin);
    void addNNHartreeXTerm(int rx, int ry, int orb, int spin);
    void addNNHartreeYTerm(int rx, int ry, int orb, int spin);

    void connectionMatrix();
    void Diagonalizer();


    Parameters_BHZ &Parameters_BHZ_;
    
    int lx_, ly_, wx_;
    int N_cells_,N_orbs_,N_spin_;
    int N_particles_, H_size, half_size_;
    double A_, B_, M_, D_, C_, Vo_, no_;
    double U_, Up_, JH_, Vnn_;

    Mat_1_doub evals_;
    Mat_2_Complex_doub C_mat;
    Mat_2_Complex_doub evecs_;
    Mat_1_Complex_doub OPs_;
    Mat_1_doub ionic;
    int ionic_size;
};

void Hamiltonian_BHZ::Initialize(){

    A_ = Parameters_BHZ_.A_val;
    B_ = Parameters_BHZ_.B_val;
    M_ = Parameters_BHZ_.M_val;
    D_ = Parameters_BHZ_.D_val;
    C_ = Parameters_BHZ_.C_val;

    lx_ = Parameters_BHZ_.Lx_;
    ly_ = Parameters_BHZ_.Ly_;
    wx_ = Parameters_BHZ_.W_;
    Vo_ = Parameters_BHZ_.Vo_val;
    no_ = Parameters_BHZ_.no_val;

    U_ = Parameters_BHZ_.U_;
	JH_ = Parameters_BHZ_.JH_;
    Up_ = Parameters_BHZ_.Up_;
    Vnn_ = Parameters_BHZ_.V_nn;
    
    N_cells_ = Parameters_BHZ_.Total_Cells;
    N_orbs_ = Parameters_BHZ_.N_Orbs;
    N_spin_ = 2;
    H_size = Parameters_BHZ_.Ham_Size;
    half_size_=(int) (H_size/2);

    evals_.resize(H_size);
    evecs_.resize(H_size);
    C_mat.resize(H_size);
    for(int i=0;i<H_size;i++){
        evecs_[i].resize(H_size);
        C_mat[i].resize(H_size);
    }

//cout << "A: " << A_ << ", B: " << B_ << ", D: " << D_ << ", M: " << M_ << endl;
//cout << "lx: " << lx_ << ", ly: " << ly_ << ", wx: " << wx_ << endl;
//cout << "N_cells: " << N_cells_ << ", N_orbs: " << N_orbs_ << ", N_spin: " << N_spin_ << endl;
//cout << "H_size: " << H_size << ", half_size: " << half_size_ << endl;
   
    OPs_.resize(H_size);
    if(Parameters_BHZ_.numberofSoftedges_==1){  ionic_size=lx_-wx_+1;   }
    else{   ionic_size=lx_-2*wx_+2;     }
    ionic.resize(ionic_size);
}

int Hamiltonian_BHZ::makeIndex(int rx, int ry, int orb, int spin){
	int r_val;
	r_val = spin * half_size_ + (orb + 2 * ry + 2 * ly_ * rx);
	return r_val;
}

int Hamiltonian_BHZ::bar_function(int val){
    int return_val;
    if (val == 0){      return_val = 1; }
	else if(val ==1){   return_val = 0; }

	return return_val;
}

void Hamiltonian_BHZ::addOnsiteTerms(int r, int orb, int spin){
	// Adding Onsite Mass Term:-->
	C_mat[r][r] += ( 1.0 * (C_ - 4.0 * D_) + 1.0 * (pow(-1.0, 1.0 * orb)) * (M_ - 4.0 * B_) ) * One_Complex;

    // Edge Confining and Ionic Potential Terms:-->
    if(!Parameters_BHZ_.canonical_){
        int rx, ry;
        ry=((r-spin*half_size_-orb)%(2*ly_))/2;
        rx=(r-spin*half_size_-orb)/(2*ly_);

        if(Parameters_BHZ_.numberofSoftedges_==1){
            if (rx > lx_ - wx_){
                C_mat[r][r] += (Vo_ * 1.0 * (rx + wx_ - lx_) / (1.0 * (wx_ - 1))) * One_Complex;
            }
            if(Parameters_BHZ_.sitewiseIonic==false){
                if(rx <= lx_ - wx_){
                    C_mat[r][r] += -1.0 * no_ * One_Complex;
                }
                if(rx == lx_- wx_){
                    C_mat[r][r] += Parameters_BHZ_.helical_edge_pot * One_Complex;
                    if(spin==0){
                        C_mat[r][r] += Parameters_BHZ_.Hfield * One_Complex;
                    }
                    else{
                        C_mat[r][r] += -1.0 * Parameters_BHZ_.Hfield * One_Complex;
                    }
                }
            }
            else{
                if(rx <= lx_ - wx_){
                    C_mat[r][r] += -1.0 * ionic[rx] * One_Complex;
                }
            }
        }
        else{
            if (rx < wx_ - 1){
                C_mat[r][r] += (Vo_ * 1.0 * (wx_ - 1 - rx) / (1.0 * (wx_ - 1))) * One_Complex;
            }
            if (rx > lx_ - wx_){
                C_mat[r][r] += (Vo_ * 1.0 * (rx + wx_ - lx_) / (1.0 * (wx_ - 1))) * One_Complex;
            }
            if(Parameters_BHZ_.sitewiseIonic==false){
                if (rx >= wx_ - 1 && rx <= lx_ - wx_){
                    C_mat[r][r] += -1.0 * no_ * One_Complex;
                }
                if(rx == lx_- wx_ || rx == wx_-1){
                    C_mat[r][r] += Parameters_BHZ_.helical_edge_pot * One_Complex;
                }
            }
            else{
                if (rx >= wx_ - 1 && rx <= lx_ - wx_){
                    C_mat[r][r] += -1.0 * ionic[rx] * One_Complex;
                }
            }            
        }

    }
    //cout <<"Onsite="<< C_mat[r][r]<<endl;
    //cout<<"Vo,no= "<<Vo_<<","<<no_<<endl;

}

void Hamiltonian_BHZ::addHoppingXTerm(int r, int rpx, int orb){
    // Adding NN Hopping connections:
    C_mat[r][rpx] = ( 1.0 * D_ + 1.0 * (pow(-1.0, 1.0 * orb)) * B_ )* One_Complex;
    C_mat[rpx][r] = conj(C_mat[r][rpx]);
}

void Hamiltonian_BHZ::addHoppingYTerm(int r, int rpy, int orb){
    // Adding NN Hopping connections:
    C_mat[r][rpy] = ( 1.0 * D_ + 1.0 * (pow(-1.0, 1.0 * orb)) * B_ ) * One_Complex;
    C_mat[rpy][r] = conj(C_mat[r][rpy]);
}

void Hamiltonian_BHZ::addOrbMixingXTerm(int r, int rpx, int orb, int spin){
    // Adding NN Orbital-Mixing connections:
    if(spin == 0){
        C_mat[rpx][r] = -(1.0 * A_ / 2.0) * Iota_Complex;
        C_mat[r][rpx] = conj(C_mat[rpx][r]);
    }
    if(spin == 1){
        C_mat[rpx][r] = (1.0 * A_ / 2.0) * Iota_Complex;
        C_mat[r][rpx] = conj(C_mat[rpx][r]);
    }
}

void Hamiltonian_BHZ::addOrbMixingYTerm(int r, int rpy, int orb, int spin){
    // Adding NN Orbital-Mixing connections:
    if(orb == 0){
        C_mat[rpy][r] = (1.0 * A_ / 2.0) * One_Complex;
        C_mat[r][rpy] = conj(C_mat[rpy][r]);
    }
    if(orb ==1){
        C_mat[rpy][r] = -(1.0 * A_ / 2.0) * One_Complex;
        C_mat[r][rpy] = conj(C_mat[rpy][r]);        
    }
}

void Hamiltonian_BHZ::addOnsiteHartreeTerm(int rx, int ry, int orb, int spin){
    //Adding the Hartree approximation of the Onsite interaction term
    int r, r_bar;
    r = makeIndex(rx, ry, orb, spin);
    r_bar = makeIndex(rx, ry, orb, bar_function(spin));

    // Hubbard Hartee Approximation:----->
    C_mat[r][r] += U_ * One_Complex * OPs_[r_bar];

    // Orbital Repulsion Hartree Approximation:----->
    int r_bar1, r_bar2;
    r_bar1 = makeIndex(rx, ry, bar_function(orb), spin);
    r_bar2 = makeIndex(rx, ry, bar_function(orb), bar_function(spin));
    
    C_mat[r][r] += (Up_ - 0.5 * JH_) * One_Complex * (OPs_[r_bar1] + OPs_[r_bar2]);

//    cout<<"U_= "<<U_<<", Up_= "<<Up_<<", JH_= "<<JH_<<endl;

    // Hund's Coupling Hartree Approximation:------> (if required)
    //only Sz*Sz term (so far)
    C_mat[r][r] += -0.5 * JH_ * One_Complex * (OPs_[r_bar1] - OPs_[r_bar2]);
}

void Hamiltonian_BHZ::addNNHartreeXTerm(int rx, int ry, int orb, int spin){
    int r;
    r = makeIndex(rx,ry,orb,spin);

    int rpx1,rpx2,rpx3,rpx4;
    if(rx!=lx_-1){
        rpx1 = makeIndex(rx+1, ry, orb, spin);
        rpx2 = makeIndex(rx+1, ry, orb, bar_function(spin));
        rpx3 = makeIndex(rx+1, ry, bar_function(orb), spin);
        rpx4 = makeIndex(rx+1, ry, bar_function(orb), bar_function(spin));

        C_mat[r][r] += Vnn_*One_Complex*(OPs_[rpx1] + OPs_[rpx2] + OPs_[rpx3] + OPs_[rpx4]);
    }

    int rmx1,rmx2,rmx3,rmx4;
    if(rx!=0){
        rmx1 = makeIndex(rx-1, ry, orb, spin);
        rmx2 = makeIndex(rx-1, ry, orb, bar_function(spin));
        rmx3 = makeIndex(rx-1, ry, bar_function(orb), spin);
        rmx4 = makeIndex(rx-1, ry, bar_function(orb), bar_function(spin));

        C_mat[r][r] += Vnn_*One_Complex*(OPs_[rmx1] + OPs_[rmx2] + OPs_[rmx3] + OPs_[rmx4]);
    }

    if(Parameters_BHZ_.PBC_X==true){
        assert(lx_ > 2 && "lx_ must be greater than 2 for PBC along-x.");

        if(rx==lx_-1){
            rpx1 = makeIndex(0, ry, orb, spin);
            rpx2 = makeIndex(0, ry, orb, bar_function(spin));
            rpx3 = makeIndex(0, ry, bar_function(orb), spin);
            rpx4 = makeIndex(0, ry, bar_function(orb), bar_function(spin));

            C_mat[r][r] += Vnn_*One_Complex*(OPs_[rpx1] + OPs_[rpx2] + OPs_[rpx3] + OPs_[rpx4]);
        }
        if(rx==0){
            rmx1 = makeIndex(lx_-1, ry, orb, spin);
            rmx2 = makeIndex(lx_-1, ry, orb, bar_function(spin));
            rmx3 = makeIndex(lx_-1, ry, bar_function(orb), spin);
            rmx4 = makeIndex(lx_-1, ry, bar_function(orb), bar_function(spin));

            C_mat[r][r] += Vnn_*One_Complex*(OPs_[rmx1] + OPs_[rmx2] + OPs_[rmx3] + OPs_[rmx4]);
        }
    }
}

void Hamiltonian_BHZ::addNNHartreeYTerm(int rx, int ry, int orb, int spin){
    int r;
    r = makeIndex(rx,ry,orb,spin);
    
    int rpy1,rpy2,rpy3,rpy4;
    if(ry!=ly_-1){
        rpy1 = makeIndex(rx, ry+1, orb, spin);
        rpy2 = makeIndex(rx, ry+1, orb, bar_function(spin));
        rpy3 = makeIndex(rx, ry+1, bar_function(orb), spin);
        rpy4 = makeIndex(rx, ry+1, bar_function(orb), bar_function(spin));

        C_mat[r][r] += Vnn_*One_Complex*(OPs_[rpy1] + OPs_[rpy2] + OPs_[rpy3] + OPs_[rpy4]);
    }

    int rmy1,rmy2,rmy3,rmy4;
    if(ry!=0){
        rmy1 = makeIndex(rx, ry-1, orb, spin);
        rmy2 = makeIndex(rx, ry-1, orb, bar_function(spin));
        rmy3 = makeIndex(rx, ry-1, bar_function(orb), spin);
        rmy4 = makeIndex(rx, ry-1, bar_function(orb), bar_function(spin));

        C_mat[r][r] += Vnn_*One_Complex*(OPs_[rmy1] + OPs_[rmy2] + OPs_[rmy3] + OPs_[rmy4]);
    }

    if(Parameters_BHZ_.PBC_Y==true){
        assert(ly_ > 2 && "ly_ must be greater than 2 for PBC along-y.");

        if(ry==ly_-1){
            rpy1 = makeIndex(rx, 0, orb, spin);
            rpy2 = makeIndex(rx, 0, orb, bar_function(spin));
            rpy3 = makeIndex(rx, 0, bar_function(orb), spin);
            rpy4 = makeIndex(rx, 0, bar_function(orb), bar_function(spin));

            C_mat[r][r] += Vnn_*One_Complex*(OPs_[rpy1] + OPs_[rpy2] + OPs_[rpy3] + OPs_[rpy4]);
        }
        if(ry==0){
            rmy1 = makeIndex(rx, ly_-1, orb, spin);
            rmy2 = makeIndex(rx, ly_-1, orb, bar_function(spin));
            rmy3 = makeIndex(rx, ly_-1, bar_function(orb), spin);
            rmy4 = makeIndex(rx, ly_-1, bar_function(orb), bar_function(spin));

            C_mat[r][r] += Vnn_*One_Complex*(OPs_[rmy1] + OPs_[rmy2] + OPs_[rmy3] + OPs_[rmy4]);
        }
    }
}

void Hamiltonian_BHZ::connectionMatrix(){

    for(int i=0;i<H_size;i++){
        for(int j=0;j<H_size;j++){
            C_mat[i][j]=Zero_Complex;
        }
    }

    int r, rpx, rpy, bar_rpx, bar_rpy;
    for (int spin = 0; spin < N_spin_; spin++){
        for(int orb = 0; orb <N_orbs_;orb++){
            for(int rx=0;rx<lx_;rx++){
                for(int ry=0;ry<ly_;ry++){

                    r   = makeIndex(rx, ry, orb, spin);
                    addOnsiteTerms(r, orb, spin);
                    addOnsiteHartreeTerm(rx, ry, orb, spin);
                    addNNHartreeXTerm(rx, ry, orb, spin);
                    addNNHartreeYTerm(rx, ry, orb, spin);

                    if(rx != lx_-1){
                        rpx = makeIndex(rx+1, ry, orb, spin);
                        addHoppingXTerm(r, rpx, orb);

                        bar_rpx = makeIndex(rx+1, ry, bar_function(orb), spin);
                        addOrbMixingXTerm(r, bar_rpx, orb, spin);
                    }
                    if(rx == lx_-1 && Parameters_BHZ_.PBC_X==true){
                        assert(lx_ > 2 && "lx_ must be greater than 2 for PBC along-x.");
                        rpx = makeIndex(0, ry, orb, spin);
                        addHoppingXTerm(r, rpx, orb);

                        bar_rpx = makeIndex(0, ry, bar_function(orb), spin);
                        addOrbMixingXTerm(r, bar_rpx, orb, spin);
                    }
                    if(ry != ly_-1){
                        rpy = makeIndex(rx, ry+1, orb, spin);
                        addHoppingYTerm(r, rpy, orb);

                        bar_rpy = makeIndex(rx, ry+1, bar_function(orb), spin);
                        addOrbMixingYTerm(r, bar_rpy, orb, spin);
                    }
                    if(ry == ly_-1 && Parameters_BHZ_.PBC_Y==true){
                        assert(ly_ > 2 && "ly_ must be greater than 2 for PBC along-y.");
                        rpy = makeIndex(rx, 0, orb, spin);
                        addHoppingYTerm(r, rpy, orb);
                        
                        bar_rpy = makeIndex(rx, 0, bar_function(orb), spin);
                        addOrbMixingYTerm(r, bar_rpy, orb, spin);
                    }
                }
            }
        }
    }

    //cout<<"Hamiltonian matrix constructed succesfully"<<endl;

/*
    for(int i=0;i<H_size;i++){
        for(int j=0;j<H_size;j++){
            cout<<C_mat[i][j]<<" ";
        }
        cout<<endl;
    }
*/
}

void Hamiltonian_BHZ::Diagonalizer(){

//    cout<<"Starting the diagonalizer"<<endl;
    std::vector<std::complex<double>> Ham_(H_size*H_size);
    
    //#pragma omp parallel for default(shared) 
    for (int i = 0; i < H_size; ++i) {
        for (int j = 0; j < H_size; ++j) {
            //Ham_[i*H_size + j] = lapack_make_complex_double(C_mat[i][j].real(), C_mat[i][j].imag());
            Ham_[i*H_size + j] = C_mat[i][j];
        }
    }

    // LAPACK routine variables
    char jobz = 'V'; // Computing both eigenvalues and eigenvectors
    char uplo = 'L'; // Using the lower triangular part of the matrix
    int n = H_size;
    int lda = H_size;
    int info;

    std::vector<double> eigs_(H_size);

    std::vector<std::complex<double>> work(1);
    //std::vector<lapack_complex_double> work(1);
    std::vector<double> rwork(3 * n ); //Should be greater than 3*n - 2;
    int lwork = -1;

    //querying with lwork=-1
    zheev_(&jobz, &uplo, &n, Ham_.data(), &lda, eigs_.data(), work.data(), &lwork, rwork.data(), &info);

//    cout<<"Diagonalization begins"<<endl;
    // Setting workspace size
    lwork = static_cast<int>(work[0].real());
    work.resize(lwork);

    // Perform the eigenvalue decomposition
    zheev_(&jobz, &uplo, &n, Ham_.data(), &lda, eigs_.data(), work.data(), &lwork, rwork.data(), &info);

    // Check for successful execution
    if (info != 0) {
        std::cerr << "LAPACK zheev_ failed with info=" << info << std::endl;
        return;
    }
    else{
    //    cout<<"Ham diagonalized succesfully"<<endl;
    }
//    cout<<"Diagonalization finished"<<endl;
    for(int i=0;i<H_size;i++){
        evals_[i] = eigs_[i];
    }
    eigs_.clear();

    for(int i=0;i<H_size;i++){
        for(int j=0;j<H_size;j++){
            evecs_[i][j] = Ham_[i*H_size+j];
        }
    }
    Ham_.clear();

/*    for (int i = 0; i < H_size; ++i) {
            for (int j = 0; j < H_size; ++j) {
                cout<<Evecs_(i,j)<<" ";
            }
            cout<<endl;
    }
*/

    string Evals_out="Eigenvalues.txt";
    ofstream Evals_file_out(Evals_out.c_str());

    for(int n=0;n<evals_.size();n++){
        Evals_file_out<<n<<"    "<<evals_[n]<<endl;
    }
}


#endif