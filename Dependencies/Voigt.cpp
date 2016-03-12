//
//  Voigt.cpp
//  Research
//
//  Created by Christopher Chan on 1/8/16.
//  Copyright © 2016 Christopher Chan. All rights reserved.
//

#include "Voigt.hpp"
#include "Faddeeva.hpp"
using namespace std;


vector<double> Voigt::voigtprofile(vector<double> x, double x_mean,double sigma, double gamma){
    //x=x-x_mean
    for(int i = 0;i<x.size();i++){
        x[i] = x[i]-x_mean;
    }
    vector<complex<double>> z;
    vector<double> returnVal;
    for(int i = 0;i<x.size();i++){
        complex<double> temp(x[i],gamma);
        temp.operator/=(sqrt(2.)*sigma);
        z.push_back(temp);
        returnVal.push_back(Faddeeva::w(z[i]).real());
    }
    return returnVal;
}
