#ifndef function_class
#define function_class

class functions{

	public:
	int ncells_,norbs_,nspin_=2,nparticles_;
	double A_,B_,M_,U_,Up_,J_,Vo_;

	Mat_2_int createBasis();

	int getUpIndex(int cell_, int orb_){	return cell_ + ncells_*orb_;					}
	int getDnIndex(int cell_, int orb_){	return cell_ + ncells_*orb_ + ncells_*norbs_;	}

	bool hasUp(const Mat_1_int &state, int cell_, int orb_){
		// spin index s=0 => up
		int ind = getUpIndex(cell_,orb_);
		return (state[ind] == 1);
	}
	bool hasDn(const Mat_1_int &state, int cell_, int orb_){
		// spin index s=1 => dn
		int ind = getDnIndex(cell_,orb_);
		return (state[ind] == 1);
	}

	//counting total occupancy in orbital "orb_" on cell "cell_"
	int occOrbital(const Mat_1_int &state, int cell_, int orb_){
		return (int) ( hasUp(state,cell_,orb_) + hasDn(state,cell_,orb_) );
	}

	//evaluating (local) Sz for orbital "orb_" on cell "cell_" = (n_up - n_down)/2
	double Sz(const Mat_1_int &state, int cell_, int orb_){
		return 0.5 * ( (hasUp(state,cell_,orb_) ? 1 : 0) - (hasDn(state,cell_,orb_) ? 1 : 0) );
		//return 0.5 * ( hasUp(state,cell_, orb_) - hasDn(state,cell_,orb_));
	}

	//evaluating (local) SpSm for orbital "orb_" on cell "cell_"
	double local_SpSm(const Mat_1_int &state, int cell_, int orb_){
		return (hasUp(state, cell_, orb_) && !hasDn(state, cell_, orb_)) ? 1.0 : 0.0;
	}

	//evaluating (local) SmSp for orbital "orb_" on cell "cell_"
	double local_SmSp(const Mat_1_int &state, int cell_, int orb_){
		return (hasDn(state, cell_, orb_) && !hasUp(state, cell_, orb_)) ? 1.0 : 0.0;
	}

	void printBasis();
};

Mat_2_int functions::createBasis(){
	int Size_ = ncells_*norbs_*nspin_;
	int Np_ = nparticles_;

	if(Np_>Size_){
		cout<<"Number of particles cannot exceed the total size of the system (ncells_*norb_*nspin_)"<<endl;
		assert(Size_>=Np_);
	}
	
	Mat_2_int basis;
	Mat_1_int combination_(Np_);
	//Initializing the first combination [0, 1, 2, ..., Np_-1]
	for(int i=0; i<Np_; i++){
		combination_[i] = i;
	}

	while(true){
		//building a state vector of length Size_ with 1's in the positions
		Mat_1_int state;
		state.resize(Size_,0);

		for(int ind : combination_){
			state[ind] = 1;
		}
		basis.push_back(state);

		//generating it's next combination in lexicographic order
		int i;
		for(i=Np_-1; i>=0; i--){
			//The key formula: combination_[i] < Size_ - (Np_ - i)
			//ensures there's room to increment.
			if(combination_[i] < Size_ - (Np_-i) ){
				combination_[i]++;
				//Resetting all subsequent elements
				for(int j=i+1; j< Np_; j++){
					combination_[j] = combination_[j-1] + 1;
				}
				break;
			}
		}
		if(i < 0){
			break;
		}

	}
	return basis;
}

void functions::printBasis(){

	Mat_2_int basis = createBasis();
	string basis_out_="Basis_Schema.txt";
	ofstream basis_out(basis_out_.c_str());

	basis_out<<"Total number of state = "<< basis.size() << endl;
	//Disclaimer: Simple beautification is performed
	for(int i=0; i<basis.size(); i++){
		basis_out << i << ": [ ";
		for(int b : basis[i]){
			basis_out << b << " ";
		}
		basis_out << "]"<< endl;
	}
}


#endif