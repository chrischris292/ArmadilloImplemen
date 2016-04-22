//
//  Sparse_Representation.h
//  Research
//
//  Created by Poop on 1/13/16.
//  Copyright © 2016 Christopher Chan. All rights reserved.
//

#ifndef Sparse_Representation_h
#define Sparse_Representation_h

#include <iostream>
#include <fstream>
#include <string>
#include "Dependencies/cnpy/cnpy.h"
#include <cmath>
#include <complex>
#include <iomanip>
#include "BigD.hpp"
#include "Dependencies/Voigt.hpp"
//#include "Header.h"

using namespace std;

class Sparse_Representation {
public:
    Sparse_Representation();
    void load_file();
    int area() {return 0;}//return width*height;}
    arma::mat prepareDictionary(double dz,int numGalaxy,int pdfSize,vector<double>& z);
    arma::mat create_voigt_dict(vector<double> &zfine, tuple<double,double> mu, int Nmu, tuple<double,double> sigma, int Nsigma, int Nv,double cut = 1.e-5);
    tuple<arma::vec,arma::vec> sparse_basis(arma::mat& dictionary,arma::vec query_vec,int n_basis, int tolerance = 0);
    double normUnsquared(arma::vec& input);

private:
    string fname = "Data/CFHTLens_sample.P.npy";
    vector<arma::vec >pdfs;
    arma::vec linspace( double a,  double b, int n);
    double norm(arma::vec & input);
    int argMax(const arma::vec & input);
    void swapVectorVar(arma::vec &input, int one, int two);
    arma::vec roundEigenVec(arma::vec &input);
    arma::vec combine_int(arma::vec index, arma::vec dind);
    int get_Ncoef(int longn);
    int getNbase(int longn);

};
Sparse_Representation::Sparse_Representation(){
    //load_file();
}

#endif /* Sparse_Representation_h */
