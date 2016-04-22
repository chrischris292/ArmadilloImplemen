#include <iostream>
#define ARMA_DONT_USE_WRAPPER
#include <armadillo>
#include "SparseRepresentation.h"
#include <time.h>
#include <chrono>

using namespace std;
using namespace arma;
using namespace std::chrono;
// Armadillo documentation is available at:
// http://arma.sourceforge.net/docs.html

void Sparse_Representation::load_file () {
    


    cnpy::NpyArray arr; 
    arr =cnpy::npy_load("Data/CFHTLens_sample.P.npy");// = cnpy::npy_load("Data/CFHTLens_sample.P.npy");
    int numGalaxy = arr.shape[0]-1; //for redshift
    int redshiftpos = arr.shape[0];
    int pdfsize = arr.shape[1];
    double* loaded_data = reinterpret_cast<double*>(arr.data);
    int curr = 0;
    for(int i = 0;i<numGalaxy;i++){
        arma::vec row(pdfsize);
        for(int i =0;i<pdfsize;i++){
            row[i] = loaded_data[curr];
            curr++;
        }
        pdfs.push_back(row);
    }
    //retrieve the last row of data for z... equivalent of z[-1]
    vector<double> z;
    for(int i =0;i<pdfsize;i++){
        z.push_back(loaded_data[curr]);
        curr++;
    }
    double dz = z[1]-z[0];
    arr.destruct();
    //int numGalaxy = 100;
    cnpy::NpyArray Dict;
    Dict= cnpy::npy_load("Data/dictionary.npy");
    double * loaded_Dict = reinterpret_cast<double*>(Dict.data);
    arma::mat D(200,60000);
    int dictI= 0;
    for(int i = 0;i<Dict.shape[0];i++){
        for(int j = 0; j<Dict.shape[1];j++){
            D(i,j) = loaded_Dict[dictI];
            dictI++;
        }
    }
    //arma::mat Dmine = prepareDictionary(dz, numGalaxy, pdfsize,z);
    
    //arma::mat Dictionary();
    //Define some tolerance and or Max number of Basis to be used, when tolerance is reached rest of basis is 0.
    double toler = 1.e-10;
    int Nsparse = 20;
    int Ncoef = 32001;
    arma::vec AA;
    AA = arma::linspace<vec>(0,1,Ncoef);
    double Da = AA[1]-AA[0];
    int s0= 0;

    BigD bigD;
    arma::vec Dvalm;
    arma::vec index;
    //Python code does this in parallel. first doing this single threaded to ensure correctness'
    high_resolution_clock::time_point t1 = high_resolution_clock::now();

    for(int ik = 0;ik<numGalaxy;ik++){
        int k = s0 + ik;
        arma::vec  pdf0 = pdfs[ik];
        int np = Nsparse;
        high_resolution_clock::time_point t2 = high_resolution_clock::now();
        tuple<arma::vec,arma::vec> Dtuple =  sparse_basis(D,pdf0,np); 
        high_resolution_clock::time_point t3 = high_resolution_clock::now();

        auto duration = std::chrono::duration_cast<microseconds>( t3 - t2 ).count();

        cout <<duration<<","<<setprecision(10);

        arma::vec Dind = get<0>(Dtuple);
        arma::vec Dval = get<1>(Dtuple);
        bigDTuple temp;

        bigD.array.push_back(temp);
        bigD.array[ik].sparse =Dtuple;
        if(Dval.max()>0){
            double dval0 = Dval[0];
            // cout << Dval;
            //cout << Dval.max();
            auto dividend =  Dval.max()*Da;
            Dvalm = (Dval/dividend);
            index = arma::round(Dvalm);
            int index0 = (int)round(dval0/Da);
            index[0] = index0;

        }
        else{
            index = arma::zeros<vec>(Dind.size());
        }
          
        bigD.array[k].sparse_ind = combine_int(index,Dind);
        /*
        string fileName = "sparseInd"+to_string(ik)+".txt";
        arma::vec sparseResult = combine_int(index,Dind);
        std::ofstream file(fileName);
         if (file.is_open())
         {
            for(int i = 0;i<sparseResult.n_elem;i++){
                file << setprecision(50)<<sparseResult[i]<<endl;
            }
         }
         file.close();
        */
        //D[:, [Dind]] = D[:, [arange(len(Dind))]]
        for(int i = 0;i<Dind.size();i++){
            D.col(Dind[i]) = D.col(i);
        }
    }
        high_resolution_clock::time_point t4 = high_resolution_clock::now();

        auto duration = chrono::duration_cast<microseconds>( t4 - t1 ).count();

        cout << "Total: "<<duration<<endl;

}
tuple<arma::vec,arma::vec> Sparse_Representation::sparse_basis(arma::mat &dictionary,arma::vec query_vec,int n_basis, int tolerance){

    arma::vec a_n(dictionary.n_cols);   //arma::zeros(dictionary.n_cols);

    arma::vec alpha = dictionary.t()*query_vec;
    arma::vec res = query_vec;
    arma::vec idxs;
    arma::vec gamma;
    int n_active;

    idxs = arma::linspace(0,dictionary.n_cols-1,dictionary.n_cols);
    //cout << idxs;
    arma::mat L = arma::zeros<arma::mat>(n_basis,n_basis);

    L(0,0) = 1;
    for( n_active = 0;n_active<n_basis;n_active++){ //9.7ms,
        //abs(dot(dictionary.T, res))
        //cout << "n_active" << n_active<<endl;
        uword lam;/*
        std::clock_t start;
        start = std::clock();
        std::clock_t startTemp;
        startTemp = std::clock();
        dictionary.t();
        std::cout << "Time for "<<n_active<<" tranpose: " << (std::clock() - startTemp) / (double)(CLOCKS_PER_SEC / 1000) << " ms" << std::endl;

        startTemp = std::clock();
        auto temp = dictionary.t()*res;
        std::cout << "Time for "<<n_active<<" multiplication: " << (std::clock() - startTemp) / (double)(CLOCKS_PER_SEC / 1000) << " ms" << std::endl;


        startTemp = std::clock();
        auto result111 = arma::abs(dictionary.t()*res);
        std::cout << "Time for "<<n_active<<" absolute value: " << (std::clock() - startTemp) / (double)(CLOCKS_PER_SEC / 1000) << " ms" << std::endl;

        //start = std::clock();
        startTemp = std::clock();

        arma::abs(dictionary.t()*res).max(lam); //12.097 ms
        std::cout << "Time for "<<n_active<<" totalFunction: " << (std::clock() - startTemp) / (double)(CLOCKS_PER_SEC / 1000) << " ms" << std::endl;


        std::clock_t start1;
        start1 = std::clock();
*/
       // high_resolution_clock::time_point t2 = high_resolution_clock::now();
        arma::abs(dictionary.t()*res).max(lam); 
       // high_resolution_clock::time_point t3 = high_resolution_clock::now();

        //auto duration = std::chrono::duration_cast<microseconds>( t3 - t2 ).count();

        //cout << "abs: "<<duration<<endl;

        //lam = argMax(dictionary.t()*res); //12.097 ms
        //std::cout << "Time for "<<n_active<<" multiplication: " << (std::clock() - start1) / (double)(CLOCKS_PER_SEC / 1000) << " ms" << std::endl;

        //cout<< "5796: "<<absVectorParam[5796]<<endl;
        //cout<< "6036: "<<absVectorParam[6036]<<endl;

        if(n_active>0){
            //L[n_active, :n_active] = dot(dictionary[:, :n_active].T, dictionary[:, lam])
            //cout << "preresultA: "<<dictionary.cols(0,n_active-1).t();
            //cout << "preresultB: "<<dictionary.col(lam);
            //cout << "Preresult: "<<result<<endl;
            //cout << "n_active: "<< n_active<<endl;
            //cout << "result shape:"<<result.n_rows << ","<<result.n_cols<<endl;
            //cout << "other shape:"<<L.row(n_active).head(n_active).n_rows << ","<<L.row(n_active).head(n_active).n_cols<<endl;
            L.row(n_active).head(n_active) = (dictionary.cols(0,n_active-1).t()*dictionary.col(lam)).t();
            //cout << "finished: "<<L.row(n_active).head(n_active);
            //cout <<"ASSIGNED VARIABLE"<<endl;
            //sla.solve_triangular(L[:n_active, :n_active], L[n_active, :n_active], lower=True, overwrite_b=True)
            //L.row(n_active).leftCols(n_active) = L.topLeftCorner(n_active,n_active).triangularView<Eigen::Lower>().solve(L.row(n_active).leftCols(n_active));
            //
            //cout << L.row(n_active).cols(0,n_active-1).n_rows<<","<<L.row(n_active).cols(0,n_active-1).n_cols<<endl;
            arma::vec b= (L.row(n_active).cols(0,n_active-1)).t();
            arma::mat a = L(span(0,n_active-1),span(0,n_active-1));

            /*cout << "TriangularA:"<<endl;
            cout<< L.topLeftCorner(n_active,n_active)<<endl;
            cout<<"TriangularB"<<endl;
            cout << b<<endl;
            cout<<"result:"<<endl;*/
            //L.row(n_active).leftCols(n_active) =L.topLeftCorner(n_active,n_active).triangularView<Eigen::Lower>().solve(L.row(n_active).leftCols(n_active).transpose()).transpose();
            //cout << "triangular A:"<<a;
            //cout << "triangular B:"<<b;
            //b =a.triangularView<Eigen::Lower>().solve(b);
            L.row(n_active).cols(0,n_active-1) = (arma::solve(arma::trimatl(a),b)).t();
            //b = arma::solve(arma::trimatl(a),y);
            //cout<<"triangular result:"<< L.row(n_active).cols(0,n_active-1)<<endl;
            //cout<<L.row(n_active).leftCols(n_active)<<endl;
            //L.row(n_active).leftCols(n_active).norm()
            //arma::vec nextRow = L.row(n_active).leftCols(n_active);
            arma::vec nextRow = (L.row(n_active).cols(0,n_active-1)).t();

            double v = pow(arma::norm(nextRow),2); //.003ms
            //cout << "V: "<<v<<endl;
            //double v = pow(L.row(n_active).leftCols(n_active).norm(),2);
            if(v<1){
                L(n_active,n_active) = sqrt(1-v);
            }
            else{
                //matrix already solved
                break;
            }
            //cout << L <<endl;
            //cout<< L.row(n_active).leftCols(n_active).cols()<<endl;
            //cout << L.row(n_active).leftCols(n_active).rows()<<endl;
            //cout <<L.row(n_active).leftCols(n_active).cols();
            //x = A.triangularView<Lower>().solve(b);
        }
        //dictionary[:, [n_active, lam]] = dictionary[:, [lam, n_active]]
        //cout<<"lamCol:"<<dictionary.col(lam);
        //cout<<"n_active_col:"<<dictionary.col(n_active);

        //dictionary.col(lam).swap(dictionary.col(n_active));
        dictionary.swap_cols(lam,n_active);
        //cout<<"n_active_col:"<<dictionary.col(n_active);
        //alpha[[n_active, lam]] = alpha[[lam, n_active]]
        swapVectorVar(alpha,lam,n_active);
        //idxs[[n_active, lam]] = idxs[[lam, n_active]]
        swapVectorVar(idxs,lam,n_active);
        //gamma = sla.cho_solve((L[:n_active + 1, :n_active + 1], True), alpha[:n_active + 1], overwrite_b=False)x = A.llt() .solve(b));  // A sym. p.d.      #include <Eigen/Cholesky>
        
        //gamma = L.topLeftCorner(n_active+1, n_active+1).llt().solve(alpha.head(n_active+1));
        arma::vec y = arma::solve(L(span(0,n_active),span(0,n_active)),alpha(span(0,n_active)));
        gamma = arma::solve(L(span(0,n_active),span(0,n_active)).t(),y);

        //res = query_vec - dot(dictionary[:, :n_active + 1], gamma)
        res = query_vec - (dictionary.cols(0,n_active)*gamma);
        //std::cout << "Time for "<<n_active<<" shorter block: " << (std::clock() - start1) / (double)(CLOCKS_PER_SEC / 1000) << " ms" << std::endl;
        //std::cout << "Time for "<<n_active<<" Total: " << (std::clock() - start) / (double)(CLOCKS_PER_SEC / 1000) << " ms" << std::endl;
    }

    tuple<arma::vec,arma::vec> result(idxs.head(n_active),gamma);

    //cout <<"idxs"<< idxs.head(n_active)<<endl; 
    //cout <<"gamma: " <<gamma<<endl;
    return result;

}
void Sparse_Representation::swapVectorVar(arma::vec &input, int one, int two){
    double temp = input[one];
    input[one] = input[two];
    input[two] = temp;
}

int Sparse_Representation::argMax(const arma::vec & input){
    int maxIndice = 0;
    for(int i = 0;i<input.size();i++){
        if(abs(input[maxIndice])<abs(input[i])){
            maxIndice = i;
        }
    }
    return maxIndice;
}

arma::vec Sparse_Representation::combine_int(arma::vec index, arma::vec dind){
    arma::vec result(index.size());
    for(int i = 0;i<index.size();i++){
        int ncoef = index[i];
        int nbase = dind[i];
        result[i] =ncoef<<16|nbase;
    }
    return result;
}
arma::vec Sparse_Representation::roundEigenVec(arma::vec & input){
    arma::vec returnVal(input.size());
    for(int i = 0; i<input.size();i++){
        returnVal[i]= round(input[i]);
    }
    return returnVal;
}
/*
arma::vec Sparse_Representation::reconstruct_pdf_int(long index, Header header, cut = 1.e-5){
    Vals =  arma::linspace<vec>(0,1,Ncoef);
    dVals = Vals[1]-Vals[0];
    //sp_ind = 
}

int Sparse_Representation::get_Ncoef(int longn){
    return (longn>>16);
}
int Sparse_Representation::getNbase(int longn){
    return (longn & (pow(2,16)-1));
}*/

int
main(int argc, char** argv)
  {    
        cout << "HI";
    Sparse_Representation init;
    init.load_file();  
  return 0;

  }

