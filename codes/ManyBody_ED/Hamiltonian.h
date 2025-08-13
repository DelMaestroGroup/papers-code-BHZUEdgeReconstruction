#include "tensors.h"
#include "functions.h"

extern "C" {
//    #include <lapacke.h>
//    #include <cblas.h>
    void zheev_(char *, char *, int *, std::complex<double> *, int *, double *, std::complex<double> *, int *, double *, int *);
}

#ifndef Hamiltonian_class
#define Hamiltonian_class

class Hamiltonian{
    public:
        Hamiltonian(functions &f__): f_(f__){
            initialize();
        }

        functions &f_;

        double A_val,B_val,M_val,Vo_;
        double U_,JbyU_,J_;
        double Up_;

        int Hsize;

        Mat_2_int basis;
        Mat_2_Complex_doub C_mat;

        Mat_2_Complex_doub Sz_mat,S2_mat;

        void initialize();
        int findStateIndex(Mat_1_int &state);
        int getFermionSign(int pos1, int pos2, Mat_1_int &state);
        void addOnsiteTerms();
        void addHoppingTerms();
        void connectionMatrix();
        void createSpinMatrices();
        void getCommutationwithHam(Mat_2_Complex_doub &Mat, string Operator_);
        void Diagonalizer();
    
        double SScorr(Mat_1_int &state, int pos1, int pos2);
        void calculateMagneticObservables();
        void matrixMultiplier(Mat_2_Complex_doub &A, Mat_2_Complex_doub &B, Mat_2_Complex_doub &C);
        void checkHermiticity(Mat_2_Complex_doub &Mat);

        Mat_1_doub evals_;
        Mat_2_Complex_doub evecs_;

};

void Hamiltonian::initialize(){
    basis = f_.createBasis();
    Hsize = basis.size();

    U_=f_.U_;
    J_=f_.J_;
    Up_=f_.Up_;

    //    J_=0.25;
    //    Up_ = 0.0;

    A_val=f_.A_;
    B_val=f_.B_;
    M_val=f_.M_;
    Vo_=f_.Vo_;

    C_mat.resize(Hsize);
    for(int i=0;i<Hsize;i++){
        C_mat[i].resize(Hsize);
    }

    Sz_mat.resize(Hsize);   S2_mat.resize(Hsize);
    for(int i=0;i<Hsize;i++){
        Sz_mat[i].resize(Hsize);    S2_mat[i].resize(Hsize);
    }

    evals_.resize(Hsize);
    evecs_.resize(Hsize);
    for(int i=0;i<Hsize;i++){
        evecs_[i].resize(Hsize);
    }
}

int Hamiltonian::findStateIndex(Mat_1_int &state){
	for(int ind=0;ind<basis.size();ind++){
		if(basis[ind]==state){
			return ind;
		}
	}
	throw runtime_error("State not found in Hamiltonian basis!");
}

int Hamiltonian::getFermionSign(int pos1, int pos2, Mat_1_int &state){
/*
    if(f_.norbs_==1){
        return +1;
    }

    //Remember the format of our states in the basis: [[...,i_{up},...],[...,j_{dn},...]]
    if( (pos1 < f_.ncells_*f_.norbs_ && pos2 < f_.ncells_*f_.norbs_) || ((pos1 >= f_.ncells_*f_.norbs_ && pos2 >= f_.ncells_*f_.norbs_)) ){
        int i,ip;
        i  = max(pos1,pos2);
        ip = min(pos1,pos2);

        int count = 0;
        for(int idx = ip+1; idx < i; idx++){
            if(state[idx] == 1){
                count++;
            }
        }
        if(count % 2 == 0){
            return +1;
        }
        else{
            return -1;
            }       
    }
    else{
        return +1;
    }
*/

    int i = std::max(pos1, pos2);
    int ip= std::min(pos1, pos2);
    int count=0;
    for (int idx = ip+1; idx < i; idx++){
        if (state[idx] == 1) count++;
    }

    double sign = std::pow(-1.0,count);

    return sign;
    
}

void Hamiltonian::addOnsiteTerms(){

    //Adding diagonal part:
    for(int index=0;index<Hsize; index++){
        Mat_1_int &state = basis[index];
        complex<double> onsite_interaction = Zero_Complex;

        for(int cell_=0; cell_<f_.ncells_;cell_++){
            for(int orb_=0;orb_ <f_.norbs_;orb_++){
                //adding Hubbard repulsion at each orb for each cell: (explicitly checked!)
                if( f_.hasUp(state,cell_,orb_) && f_.hasDn(state,cell_,orb_) ){
                    onsite_interaction += U_*One_Complex;
                }
            }

            if(f_.norbs_ > 1){
                //adding orb-orb repulsion for each cell: (explicitly checked!)
                int n0 = f_.occOrbital(state, cell_, 0);
                int n1 = f_.occOrbital(state, cell_, 1);
                onsite_interaction += (Up_-0.5*J_)*(1.0*n0*n1)*One_Complex;

                //adding onsite energy at each orb for each cell: (explicitly checked!)
                onsite_interaction += ((Vo_/3.0) + (M_val-4.0*B_val))*(1.0*n0)*One_Complex;
                onsite_interaction += ((Vo_/3.0) - (M_val-4.0*B_val))*(1.0*n1)*One_Complex;

                //adding sz-sz Hund's coupling term for each cell: (explicitly checked!)
                double sz0 = f_.Sz(state,cell_,0);
                double sz1 = f_.Sz(state,cell_,1);
                onsite_interaction += (-2.0*J_)*(1.0*sz0*sz1)*One_Complex;

                //cout<<"Sz-0 = "<<sz0<<", Sz-1 = "<<sz1<<endl;
            }
        }
        C_mat[index][index] += onsite_interaction;
    }

    if(f_.norbs_ > 1){ 
        //Adding off-diagonal part: -J(S^+_0S^-_1 + S^-_0S^+_1 )
        //S^+_0 S^-_1 flips "0-dn" to "0-up" and "1-up" to "1-dn"
        //S^-_0 S^+_1 flips "0-up" to "0-dn" and "1-dn" to "1-up"

        for(int index=0;index<Hsize; index++){
            Mat_1_int &old_state = basis[index];

            for(int cell_=0;cell_<f_.ncells_;cell_++){
                //S^+_0 S^-_1: c^{dag}_{cell,0,up}c_{cell,0,dn}c^{dag}_{cell,1,dn}c_{cell,1,up}: (explicitly checked!)
                if( f_.hasDn(old_state, cell_, 0) && (!f_.hasUp(old_state, cell_, 0)) &&
                    (!f_.hasDn(old_state, cell_, 1)) && f_.hasUp(old_state, cell_, 1) ){

                    Mat_1_int new_state = old_state;

                    //remove up from (cell,1), add dn to (cell,1)
                    new_state[f_.getDnIndex(cell_,1)] = 1;
                    new_state[f_.getUpIndex(cell_,1)] = 0;
                    //remove dn from (cell,0), add up to (cell,0)
                    new_state[f_.getDnIndex(cell_,0)] = 0;
                    new_state[f_.getUpIndex(cell_,0)] = 1;

                    double fermionSign=1.0;

                    fermionSign *= getFermionSign(f_.getDnIndex(cell_,0),f_.getUpIndex(cell_,0), old_state);
                    fermionSign *= getFermionSign(f_.getUpIndex(cell_,1),f_.getDnIndex(cell_,1), old_state);

                    int new_index = findStateIndex(new_state);

                //    cout<<"annihilating-state= "<<index<<", creating-state= "<<new_index<<endl;

                    C_mat[new_index][index] += (-J_*1.0*fermionSign)*One_Complex;
                }

                //S^-_0 S^+_1: c^{dag}_{cell,0,dn}c_{cell,0,up}c^{dag}_{cell,1,up}c_{cell,1,dn}: (explicitly checked!)
                if( f_.hasUp(old_state, cell_, 0) && (!f_.hasDn(old_state, cell_, 0)) &&
                    (!f_.hasUp(old_state, cell_, 1)) && f_.hasDn(old_state, cell_, 1) ){

                    Mat_1_int new_state = old_state;

                    //remove up from (cell,1), add dn to (cell,1)
                    new_state[f_.getDnIndex(cell_,1)] = 0;
                    new_state[f_.getUpIndex(cell_,1)] = 1;
                    //remove dn from (cell,0), add up to (cell,0)
                    new_state[f_.getDnIndex(cell_,0)] = 1;
                    new_state[f_.getUpIndex(cell_,0)] = 0;

                    double fermionSign = 1.0;

                    fermionSign *= getFermionSign(f_.getUpIndex(cell_,0),f_.getDnIndex(cell_,0), old_state);
                    fermionSign *= getFermionSign(f_.getDnIndex(cell_,1),f_.getUpIndex(cell_,1), old_state);

                    int new_index = findStateIndex(new_state);

                //    cout<<"annihilating-state= "<<index<<", creating-state= "<<new_index<<endl;

                    C_mat[new_index][index] += (-J_*1.0*fermionSign)*One_Complex;
                }
            }
        }

        //Adding Pair-hopping Term:  +J(c^{dag}_{cell,1,dn}c_{cell,1,up}c^{dag}_{cell,0,up}c_{cell,0,dn} + h.c.)
        for(int index=0;index<Hsize; index++){
            Mat_1_int &old_state = basis[index];

            for(int cell_=0;cell_<f_.ncells_;cell_++){
                //S^+_0 S^-_1: c^{dag}_{cell,0,up}c_{cell,0,dn}c^{dag}_{cell,1,dn}c_{cell,1,up}: (explicitly checked!)
                if( f_.hasDn(old_state, cell_, 0) && (f_.hasUp(old_state, cell_, 0)) &&
                    (!f_.hasDn(old_state, cell_, 1)) && !f_.hasUp(old_state, cell_, 1) ){

                    Mat_1_int new_state = old_state;

                    //remove up and dn from (cell,0)
                    new_state[f_.getDnIndex(cell_,0)] = 0;
                    new_state[f_.getUpIndex(cell_,0)] = 0;
                    //remove up and dn to (cell,1)
                    new_state[f_.getDnIndex(cell_,1)] = 1;
                    new_state[f_.getUpIndex(cell_,1)] = 1;

                    int new_index = findStateIndex(new_state);

                    double fermionSign = 1.0;

                    fermionSign *= getFermionSign(f_.getUpIndex(cell_,0),f_.getUpIndex(cell_,1), old_state);
                    fermionSign *= getFermionSign(f_.getDnIndex(cell_,0),f_.getDnIndex(cell_,1), old_state);

                    C_mat[new_index][index] += (J_*1.0*fermionSign)*One_Complex;
                    C_mat[index][new_index] += conj(C_mat[new_index][index]);
                }
            }
        }
    }
}

void Hamiltonian::addHoppingTerms(){

    int pos_from, pos_to;
    for (int index=0; index<Hsize; index++){
        Mat_1_int old_state = basis[index];

        for (int cell= 0; cell<f_.ncells_; cell++){
            int next_cell = (cell + 1) % f_.ncells_; //to include PBC in this loop

            //IntraOrbital Hopping: (explicitly checked!)
            for (int orb=0; orb<f_.norbs_; orb++){
                for (int spin=0; spin<f_.nspin_; spin++){
                    
                    if(spin == 0){
                        pos_from = f_.getUpIndex(cell, orb);
                        pos_to   = f_.getUpIndex(next_cell, orb);
                    }
                    else{
                        pos_from = f_.getDnIndex(cell, orb);
                        pos_to   = f_.getDnIndex(next_cell, orb);
                    }

                    //if there's an electron at pos_from and an empty site at pos_to, then:
                    if(old_state[pos_from] == 1 && old_state[pos_to] == 0){
                        
                        Mat_1_int new_state = old_state;
                        new_state[pos_from] = 0;
                        new_state[pos_to]   = 1;

                        int new_index = findStateIndex(new_state);

                        double fermionSign = getFermionSign(pos_from, pos_to, old_state);                        
                        
                        double hop_intra;
                        if(orb == 0){   hop_intra =  B_val;    }
                        else{           hop_intra = -B_val;   }
                        
                        C_mat[new_index][index] += hop_intra * fermionSign *  One_Complex;
                        C_mat[index][new_index] += conj(C_mat[new_index][index]);

                    }
                }
            }

            if(f_.norbs_ > 1){
            //InterOrbital Hopping: (explicitly checked!)
            for(int spin= 0; spin<f_.nspin_; spin++){
                //Term-1: A/2 ( c^{\dagger}_{i+1,0,\sigma} c_{i,1,\sigma} + h.c.)
                if(spin == 0){
                    pos_from = f_.getUpIndex(cell, 1); 
                    pos_to   = f_.getUpIndex(next_cell, 0); 
                }
                else{
                    pos_from = f_.getDnIndex(cell, 1); 
                    pos_to   = f_.getDnIndex(next_cell,0);
                }
                
                if(old_state[pos_from] == 1 && old_state[pos_to] == 0){
                    Mat_1_int new_state = old_state;
                    new_state[pos_from] = 0;
                    new_state[pos_to]   = 1;

                    int new_index = findStateIndex(new_state);
                    double fermionSign = getFermionSign(pos_from, pos_to, old_state);
                    
                    C_mat[new_index][index] += (A_val/2.0) * fermionSign * One_Complex;
                    C_mat[index][new_index] += conj(C_mat[new_index][index]);
                }

                //Term-1: -A/2 ( c^{\dagger}_{i+1,1,\sigma} c_{i,0,\sigma} + h.c.)
                if(spin == 0){
                    pos_from = f_.getUpIndex(cell, 0);
                    pos_to   = f_.getUpIndex(next_cell,1);
                }
                else{
                    pos_from = f_.getDnIndex(cell, 0);
                    pos_to   = f_.getDnIndex(next_cell,1);
                }

                if(old_state[pos_from] == 1 && old_state[pos_to] == 0){
                    Mat_1_int new_state = old_state;
                    new_state[pos_from] = 0;
                    new_state[pos_to]   = 1;
                    
                    int new_index = findStateIndex(new_state);
                    double fermionSign = getFermionSign(pos_from, pos_to, old_state);
                    
                    C_mat[new_index][index] += -(A_val/2.0) * fermionSign * One_Complex;
                    C_mat[index][new_index] += conj(C_mat[new_index][index]);
                }
            }
        
            }
        }
    }

}

void Hamiltonian::connectionMatrix(){
    for(int i=0;i<Hsize;i++){
        for(int j=0;j<Hsize;j++){
            C_mat[i][j]=Zero_Complex;
        }
    }

    addOnsiteTerms();
    addHoppingTerms();

/*    for(int i=0;i<Hsize;i++){
        for(int j=0;j<Hsize;j++){
            cout<<C_mat[i][j]<<" ";
        }
        cout<<endl;
    }*/

    checkHermiticity(C_mat);
}

void Hamiltonian::createSpinMatrices(){
    for(int i=0;i<Hsize;i++){
        for(int j=0;j<Hsize;j++){
            Sz_mat[i][j] = Zero_Complex;
            S2_mat[i][j] = Zero_Complex;
        }
    }

    //Sz matrix:
    for(int index=0;index<Hsize;index++){
        double valSz = 0.0;
        for(int cell=0;cell<f_.ncells_; cell++){
            for(int orb=0;orb<f_.norbs_;orb++){

                valSz += f_.Sz(basis[index], cell, orb);
            }
        }
        Sz_mat[index][index] = valSz*One_Complex;
    }

    //S+, S- matrices:
    Mat_2_Complex_doub Sp_mat(Hsize, Mat_1_Complex_doub (Hsize, Zero_Complex));
    Mat_2_Complex_doub Sm_mat(Hsize, Mat_1_Complex_doub (Hsize, Zero_Complex));

    for(int index=0;index<Hsize;index++){
        Mat_1_int &oldState = basis[index];

        for(int cell=0;cell<f_.ncells_;cell++){
            for(int orb=0;orb<f_.norbs_;orb++){
                
                if( f_.hasDn(oldState, cell, orb) && !f_.hasUp(oldState, cell, orb) ){

                    Mat_1_int newState = oldState;
                    newState[f_.getDnIndex(cell,orb)] = 0;
                    newState[f_.getUpIndex(cell,orb)] = 1;
                    
                    int newIndex = findStateIndex(newState);
                    double fermionSign = getFermionSign(f_.getDnIndex(cell,orb),f_.getUpIndex(cell,orb), oldState);
                //    double fermionSign=1.0;
                                        
                    Sp_mat[newIndex][index] += 1.0*fermionSign*One_Complex;
                }
                if( f_.hasUp(oldState, cell, orb) && !f_.hasDn(oldState, cell, orb) ){

                    Mat_1_int newState = oldState;
                    newState[f_.getUpIndex(cell,orb)] = 0;
                    newState[f_.getDnIndex(cell,orb)] = 1;
                    
                    int newIndex = findStateIndex(newState);
                    double fermionSign = getFermionSign(f_.getUpIndex(cell,orb),f_.getDnIndex(cell,orb), oldState);
                //    double fermionSign=1.0;

                    Sm_mat[newIndex][index] += 1.0*fermionSign*One_Complex;
                }
            }
        }
    }

    Mat_2_Complex_doub Sx_mat(Hsize, Mat_1_Complex_doub (Hsize, Zero_Complex));
    Mat_2_Complex_doub Sy_mat(Hsize, Mat_1_Complex_doub (Hsize, Zero_Complex));

    for(int ind1=0;ind1<Hsize;ind1++){
        for(int ind2=0;ind2<Hsize;ind2++){
            Sx_mat[ind1][ind2] = (One_Complex/2.0)*(Sp_mat[ind1][ind2] + Sm_mat[ind1][ind2]);
            Sy_mat[ind1][ind2] = (1.0/(2.0*Iota_Complex))*(Sp_mat[ind1][ind2] - Sm_mat[ind1][ind2]);
        }
    }

    Mat_2_Complex_doub Sx2_mat(Hsize, Mat_1_Complex_doub (Hsize, Zero_Complex));
    Mat_2_Complex_doub Sy2_mat(Hsize, Mat_1_Complex_doub (Hsize, Zero_Complex));
    Mat_2_Complex_doub Sz2_mat(Hsize, Mat_1_Complex_doub (Hsize, Zero_Complex));

    matrixMultiplier(Sx_mat, Sx_mat, Sx2_mat);
    matrixMultiplier(Sy_mat, Sy_mat, Sy2_mat);
    matrixMultiplier(Sz_mat, Sz_mat, Sz2_mat);

    for(int ind1=0;ind1<Hsize;ind1++){
        for(int ind2=0;ind2<Hsize;ind2++){
            S2_mat[ind1][ind2] = Sz2_mat[ind1][ind2] + Sx2_mat[ind1][ind2] + Sy2_mat[ind1][ind2];
        }
    }

    Sp_mat.clear();     Sm_mat.clear();
    Sx_mat.clear();     Sy_mat.clear();
    Sx2_mat.clear();    Sy2_mat.clear();    Sz2_mat.clear();

    string file_S2_="S2_matrix.txt";
    ofstream file_S2(file_S2_.c_str());

    for(int i=0;i<Hsize;i++){
        for(int j=0;j<Hsize;j++){
            file_S2<<S2_mat[i][j]<<" ";
        }
        file_S2<<endl;
    }

    getCommutationwithHam(Sz_mat, "S^z");
    getCommutationwithHam(S2_mat, "S^2");
}

void Hamiltonian::matrixMultiplier(Mat_2_Complex_doub &A, Mat_2_Complex_doub &B, Mat_2_Complex_doub &C){
    int dim = A.size();
    //Assuming A,B,C are dxd matrices:
    complex<double> sumVal;
    for(int i=0;i<dim;i++){
        for(int j=0;j<dim;j++){
            sumVal=Zero_Complex;

            for(int k=0;k<dim;k++){
                sumVal += A[i][k]*B[k][j];
            }
            C[i][j] = sumVal;
        }
    }    
}

void Hamiltonian::getCommutationwithHam(Mat_2_Complex_doub &Mat, string Operator_){
    Mat_2_Complex_doub HO_mat(Hsize, Mat_1_Complex_doub (Hsize, Zero_Complex));
    Mat_2_Complex_doub OH_mat(Hsize, Mat_1_Complex_doub (Hsize, Zero_Complex));
    Mat_2_Complex_doub Comm_mat(Hsize, Mat_1_Complex_doub (Hsize, Zero_Complex));

    matrixMultiplier(C_mat, Mat, HO_mat);
    matrixMultiplier(Mat, C_mat, OH_mat);

/*    for(int i=0; i<Hsize; i++){
        for(int j=0; j<Hsize; j++){
            Comm_mat[i][j] = HO_mat[i][j] - OH_mat[i][j];
            cout<<Comm_mat[i][j]<<" ";
        }
        cout<<endl;
    }
*/
//    checkHermiticity(Comm_mat);

    double normComm=0.0;
    for(int i=0; i<Hsize; i++){
        for(int j=0; j<Hsize; j++){
            normComm +=norm(Comm_mat[i][j]);
        }
    }
    normComm= sqrt(normComm);
    cout<<"Norm of commutator [H, " <<Operator_<< " ] = "<<normComm<<endl;
}

void Hamiltonian::Diagonalizer(){

//    cout<<"Starting the diagonalizer"<<endl;
    std::vector<std::complex<double>> Ham_(Hsize*Hsize);
//    cout<<Hsize<<endl;
//    cout<<"1"<<endl;

    //#pragma omp parallel for default(shared) 
    for (int i = 0; i < Hsize; ++i) {
        for (int j = 0; j < Hsize; ++j) {
            //Ham_[i*H_size + j] = lapack_make_complex_double(C_mat[i][j].real(), C_mat[i][j].imag());
            Ham_[i*Hsize + j] = C_mat[i][j];
        }
    }

    // LAPACK routine variables
    char jobz = 'V'; // Computing both eigenvalues and eigenvectors
    char uplo = 'L'; // Using the lower triangular part of the matrix
    int n = Hsize;
    int lda = Hsize;
    int info;

    std::vector<double> eigs_(Hsize);

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
    for(int i=0;i<Hsize;i++){
        evals_[i] = eigs_[i];
    }
    eigs_.clear();

    for(int i=0;i<Hsize;i++){
        for(int j=0;j<Hsize;j++){
            evecs_[i][j] = Ham_[i*Hsize+j];
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

double Hamiltonian::SScorr(Mat_1_int &state, int pos1, int pos2){

    double corr=0.0;

    int cell1, cell2, orb1, orb2;
    cell1 = pos1%f_.ncells_;    orb1 = pos1/f_.ncells_;
    cell2 = pos2%f_.ncells_;    orb2 = pos2/f_.ncells_;

    double Sz1 = f_.Sz(state, cell1, orb1);
    double Sz2 = f_.Sz(state, cell2, orb2);

    corr += Sz1 * Sz2;

    //S^+_1 S^-_2 => 1 if pos-1 has down and pos-2 has up:
    bool pos1_hasDn = f_.hasDn(state, cell1, orb1);
    bool pos2_hasUp = f_.hasUp(state, cell2, orb2);
    int SpSm = (pos1_hasDn && pos2_hasUp) ? 1 : 0;

    //S^-_1 S^+_2 => 1 if pos-1 has up and pos-2 has dn:
    bool pos1_hasUp = f_.hasUp(state, cell1, orb1);
    bool pos2_hasDn = f_.hasDn(state, cell2, orb2);
    int SmSp = (pos1_hasUp && pos2_hasDn) ? 1 : 0;
    
    corr += 0.5*(SpSm + SmSp);
    return corr;
}

void Hamiltonian::calculateMagneticObservables(){

    string Magnetic_out_="Magnetic_observables_for_each_state.txt";
    ofstream Magnetic_out(Magnetic_out_.c_str());
    Magnetic_out<<"#index    Sz_tot    S2_local    S2_tot"<<endl;

    int total_sites = f_.ncells_*f_.norbs_;

    Mat_1_doub Sz_total(Hsize,0.0),S2_local(Hsize,0.0),S2_total(Hsize,0.0);
    
    //evaluating total Sz for each eigenvector |m>: sum_{i} <m|Sz_i|m>=sum_{b}|Psi_m(b)|^2 sum_{i}Sz_i(b)
    //evaluating local S2 for each eigenvector |m>: sum_{i} <m|S2_i|m>=sum_{i}<m|(Sz_i)^2|m> + 0.5<m|(Sp_i*Sm_i + Sm_i*Sp_i)|m>
    //(a) sum_{i}<m|(Sz_i)^2|m> = sum_{b}|Psi_m(b)|^2 sum_{i}(Sz_i(b))^2
    //(b) sum_{i}<m|Sp_i*Sm_i|m> = sum_{b}|Psi_m(b)|^2 sum_{i}(Sp_i*Sm_i)(b)
    //(c) similar to (b)
    for(int n=0; n<Hsize; n++){
        double Total_Sz_n=0.0;
        double Local_S2_n=0.0;
        double Corr_SiSj_n=0.0;

        for(int m=0;m<Hsize;m++){
            double prob = norm(evecs_[n][m]);

            double stateSz = 0.0;
            double stateS2 = 0.0;
            
            for(int cell=0; cell<f_.ncells_; cell++){
                for (int orb = 0; orb < f_.norbs_; orb++){
                    
                    stateSz += f_.Sz(basis[m], cell, orb);

                    //SzSz term:
                    stateS2 += f_.Sz(basis[m], cell, orb)*f_.Sz(basis[m], cell, orb);

                    //0.5*(SpSm + SmSp) term:
                    stateS2 += 0.5*(f_.local_SpSm(basis[m], cell, orb) + f_.local_SmSp(basis[m], cell, orb));

                    //0.5(SpSm + SmSp) term: will be implemented
//                    stateSpSm += f_.local_SpSm(basis[m], cell, orb);
//                    stateSmSp += f_.local_SmSp(basis[m], cell, orb);
                }
            }
            Total_Sz_n += prob * stateSz;
            Local_S2_n += prob * stateS2;


            //evaluating total S^2 for each eigenvector: <m|S^2|m> = sum_{i}<m|S^2_{i}|m> + sum_{i,j}<m|S_{i}.S_{j}|m>
            double stateSiSj = 0.0;
            for(int site1=0; site1<total_sites; site1++){
                for (int site2=0; site2<total_sites; site2++){

                    if(site1 != site2){
                        stateSiSj +=  SScorr(basis[m], site1, site2);
                    }
                }
            }
            Corr_SiSj_n += prob * stateSiSj;
        }
        Sz_total[n] = Total_Sz_n;
        S2_local[n] = Local_S2_n;
        S2_total[n] = Local_S2_n + Corr_SiSj_n;

        Magnetic_out<<n<<"      "<<Sz_total[n]<<"           "<<S2_local[n]<<"       "<<S2_total[n]<<endl;
    }

}

void Hamiltonian::checkHermiticity(Mat_2_Complex_doub &Mat){
    int dim = Mat.size();

    for(int i=0;i<dim;i++){
        for(int j=0;j<dim;j++){
            if( abs(Mat[i][j]-conj(Mat[j][i]))>1e-6 ){
                cout<<"Error: Matrix is not Hermitian. Check Mat["<<i <<"]["<<j <<"]"<<endl;
            }
            assert(abs(Mat[i][j]-conj(Mat[j][i]))<1e-6 );
        }
    }
}


#endif
