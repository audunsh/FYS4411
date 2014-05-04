#include <basis.h>


basis::basis(int N)
{
    Nstates = N;
    Nstates2 = 2*Nstates;
    field<mat> vn(Nstates, Nstates); //no-spin basis
    for (int i = 0; i < Nstates; ++i) {
        for (int j = 0; j < Nstates; ++j) {
            vn(i,j) = zeros<mat>(Nstates,Nstates); // fill V with 3x3 mx elements
        }
    }

    field<mat> Vn(Nstates2, Nstates2); //spin-basis
    for (int i = 0; i < Nstates2; ++i) {
        for (int j = 0; j < Nstates2; ++j) {
            Vn(i,j) = zeros<mat>(Nstates2,Nstates2); // fill V with 3x3 mx elements
        }
    }
    v = vn;  // global variables for class basis.cpp
    V = Vn;
    h.set_size(Nstates,Nstates);
    h.zeros();
    H.set_size(Nstates2,Nstates2);
    H.zeros();

    //creating the overlap matrix
    S.set_size(Nstates,Nstates);
    S.zeros();
}

void basis::set_orthonormal(){
    //sets the overlap matrix equal to the identity matrix
    for(int i=0;i < Nstates2;i++){
        for(int j=0;j < Nstates2;j++){
            if(i==j){
                S(i,j) = 1.0;
            }
        }
    }
}

void basis::init_overlap(){
    //calculates the overlap matrix for the generated basis
    for(int i=0;i < Nstates2;i++){
        for(int j=0;j < Nstates2;j++){
            if(i==j){
                S(i,j) = 1.0;
            }
        }
    }
}

void basis::read(string filename, int Zn){
    Z = Zn;
    //read basis from file
    //this also means you have to provide the one-body energy and the overlap matrix (use h0 and ...)
    cout << "Loading predefined basis from file: " << filename << endl;
    ifstream myfile;
    myfile.open(filename.c_str());
    if(!myfile.good()){
        cout << "file is not open" << endl;
    }
    if (myfile.is_open()){
        int p,q,r,s;
        double value;
        while (!myfile.eof()){
            myfile >> p;
            myfile >> q;
            myfile >> r;
            myfile >> s;
            myfile >> value;
            //cout << p << q << r << s << endl;
            try{
                v(p,q)(r,s) = Z*value;
                if(p==q){
                    h(p,q) = h0(p*2,q*2);
                }
            }
            catch(int e){
                cout << "Failed to load basis." << endl;
            }

        }
    }
    else
        cout << "Did not manage to open file in HFSolve::init()"<< endl;
    expand();
}

void basis::expand(){
    h.set_size(Nstates2);
    S.set_size(Nstates2,Nstates2);
    S.zeros();
    //expanding basis to include spin
    double D = 0;
    double Ex = 0;
    for (int p = 0; p < Nstates2; ++p) {
        for (int q = 0; q < Nstates2; ++q) {
            for (int r = 0; r < Nstates2; ++r) {
                for (int s= 0; s < Nstates2; ++s) {
                    try{
                        D = v(p/2,q/2)(r/2,s/2);  // Direct term
                        Ex =v(p/2,q/2)(s/2,r/2); // Exchange term
                        //cout << D << Ex << " ";
                        //cout << p << " " << q << " " << r << " " << s << endl;

                        V(p,q)(r,s) = state(p,q,r,s,D,Ex);
                    }
                    catch(int e){
                        cout << "Failed to read basis." << endl;
                    }
                }
            }
            H(p,q) = h0(p,q);
        }
    }
}

double basis::h0(int i, int j){
    // the one-body interaction
    double h = 0;
    if (i == j){
        double n = i/2 + 1.0;
        h = -(Z*Z)/(2*n*n);
    }
    return h;
}

double basis::state(int p, int q, int r, int s, double D, double Ex){
    //Evaluating spin configuration, returning direct and/or exchange term or 0
    double S = 0;
    int s1 = 0; int s2 = 0;int s3 = 0;int s4 = 0;
    s1 = p%2;
    s2 = q%2;
    s3 = r%2;
    s4 = s%2;
    if (s1 == s2){
        if (s3 == s4){
            if ( s1 == s3){
                S = D-Ex;
            }
            else{
                S = 0;
            }
        }
    }
    else if (s1 != s2){
        if (s3 != s4){
            if (s1 == s3){
                S = D;
            }
            else{
                S = -Ex;
            }
        }
        else{
            S = 0;
        }
    }
    return S;
}


double basis::get(int p, int q, int r, int s){
    //returns matrix element pqrs
}

void basis::generate(){
    //function to generate the basis using STO/GTOs
}

double basis::eval(int p, int q, int r, int s){
    //function to calculate the matrix elements, possibly on the fly
    //This function also needs a variable to specify type of orbitals
}

