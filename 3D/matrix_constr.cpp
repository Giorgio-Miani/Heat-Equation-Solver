#include <chrono>
#include <cstdlib>
#include <iostream>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <unsupported/Eigen/SparseExtra>

    using SpMat = Eigen::MatrixXd;
    using SpVec = Eigen::VectorXd;


int main(int argc, char* argv[])
{
    int const N = 300;
    int const dim = N*N*N;
    // SpMat M(dim, dim);
    // Eigen::SparseMatrix<double> dense_to_sparceM(dim,  dim);    
    Eigen::SparseMatrix<double> sparseM(dim,  dim);
    std::vector<Eigen::Triplet<double>> triplets;
    double deltaT = 0.1;
    double deltaX = 0.1;
    double coef = 1/(2*deltaX*deltaX);


    //Matrix construction using dense matrix and trasforming it into a sparse one
    //Comment for N > 30
    auto start_time
            = std::chrono::system_clock::now();

    // for (int k = 0; k < N; k++){
    //     for(int j = 0; j < N; j++){
    //         for(int i = 0; i < N; i ++){
    //             M(i + (j)*N  + (k)*N*N, i + (j)*N  + (k)*N*N) += 1/deltaT + 3/(deltaX*deltaX);
    //         if(i > 0)
    //             M(i + (j)*N  + (k)*N*N, i - 1 + (j)*N  + (k)*N*N) += coef;
    //         if(i < N - 1)
    //             M(i + (j)*N  + (k)*N*N, i + 1 + (j)*N  + (k)*N*N) += coef;
    //         if(j > 0)
    //             M(i + (j)*N  + (k)*N*N, i + (j - 1)*N  + (k)*N*N) += coef;
    //         if(j < N - 1)
    //             M(i + (j)*N  + (k)*N*N, i + (j + 1)*N  + (k)*N*N) += coef;
    //         if(k > 0)
    //             M(i + (j)*N  + (k)*N*N, i + (j)*N  + (k - 1)*N*N) += coef;
    //         if(k < N - 1)
    //             M(i + (j)*N  + (k)*N*N, i + (j)*N  + (k + 1)*N*N) += coef;
    //         }
    //     }
    // }

    // dense_to_sparceM = M.sparseView();

    auto end_time
        = std::chrono::system_clock::now();

    std::chrono::duration<double> dense_to_sparse 
        = end_time - start_time;

    //Matrix construction using triplets and sparse matrix
    
    start_time
        = std::chrono::system_clock::now();

    for (int k = 0; k < N; k++){
        for(int j = 0; j < N; j++){
            for(int i = 0; i < N; i ++){
            triplets.push_back(Eigen::Triplet<double>(i + (j)*N  + (k)*N*N, i + (j)*N  + (k)*N*N, 1/deltaT + 3/(deltaX*deltaX))); 
            if(i > 0)
                triplets.push_back(Eigen::Triplet<double>(i + (j)*N  + (k)*N*N, i - 1 + (j)*N  + (k)*N*N, coef));
            if(i < N - 1)
                triplets.push_back(Eigen::Triplet<double>(i + (j)*N  + (k)*N*N, i + 1 + (j)*N  + (k)*N*N, coef)); 
            if(j > 0)
                triplets.push_back(Eigen::Triplet<double>(i + (j)*N  + (k)*N*N, i + (j - 1)*N  + (k)*N*N, coef) );
            if(j < N - 1)
                triplets.push_back(Eigen::Triplet<double>(i + (j)*N  + (k)*N*N, i + (j + 1)*N  + (k)*N*N, coef) );
            if(k > 0)
                triplets.push_back(Eigen::Triplet<double>(i + (j)*N  + (k)*N*N, i + (j)*N  + (k - 1)*N*N, coef) );
            if(k < N - 1)
                triplets.push_back(Eigen::Triplet<double>(i + (j)*N  + (k)*N*N, i + (j)*N  + (k + 1)*N*N, coef) );
            }
        }
    }

    sparseM.setFromTriplets(triplets.begin(), triplets.end());

    end_time
        = std::chrono::system_clock::now();

    std::chrono::duration<double> sparse_to_sparse 
        = end_time - start_time;

    std::cout << "Total time dense to sparse: "
                    << dense_to_sparse.count() << " seconds " << std::endl << std::endl;
    std::cout << "Total time sparse to sparse: "
                << sparse_to_sparse.count() << " seconds " << std::endl << std::endl;

    

    //Risolvere con 1milione di punti

    // std::cout << sparseM << std::endl;

}
