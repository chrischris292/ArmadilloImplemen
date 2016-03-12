#include <iostream>
#include <armadillo>
#include "SparseRepresentation.h"
using namespace std;
using namespace arma;

// Armadillo documentation is available at:
// http://arma.sourceforge.net/docs.html




int
main(int argc, char** argv)
  {    
    arma::mat A = arma::mat(3,3);
    arma::vec B = arma::vec(3);
    
    B(0)=0.26165705285168211480240074706671293824911117553711;
    B(1)=0.12205090710037808099386325011437293142080307006836;
    B(2) = 0.07786977;
    A(0,0) = 1.;
    A(0,1)=  0.;
    A(0,2) = 0.;
    A(1,0) = 0.42448897901098936458197385945823043584823608398438;
    A(1,1) = 0.90543310448547653646045318964752368628978729248047;
    A(1,2) = 0.;
    A(2,0) =  0.32514374;
    A(2,1) = -0.15169085;
    A(2,2) = 0.93341922;

    arma::vec y = arma::solve(A,B);
    cout << arma::solve(A.t(),y);
  
  return 0;

  }

