#include <chrono>
#include <cstdlib>
#include <iostream>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <unsupported/Eigen/SparseExtra>
#include <Eigen/IterativeLinearSolvers>

using SpMat = Eigen::SparseMatrix<double>;
using SpVec = Eigen::VectorXd;
using triplet = Eigen::Triplet<double>;

void system_assemble(SpMat &A, SpVec &rhs, const int N, const double deltaX)
{

    std::vector<Eigen::Triplet<double>> triplets;
    double const diag_coef = 6/(deltaX*deltaX);
    double const coef = -1/(deltaX*deltaX);

    for (int k = 0; k < N; k++){
        for(int j = 0; j < N; j++){
            for(int i = 0; i < N; i ++){
                triplets.push_back(Eigen::Triplet<double>(i + (j)*N  + (k)*N*N, i + (j)*N  + (k)*N*N, diag_coef)); 
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


                //Righten side
                rhs(i + (j)*N  + (k)*N*N) = 3 * std::sin((i + 1) * deltaX) * std::sin((j + 1) * deltaX)*std::sin((k + 1) * deltaX);

                // Boundary conditions

                if(i == N - 1){
                    rhs(i + (j)*N  + (k)*N*N) +=  (std::sin(1.0) * std::sin((j + 1) * deltaX) * std::sin((k + 1) * deltaX))/(deltaX*deltaX);// i == N-1
                }
                if(j == N - 1){
                    rhs(i + (j)*N  + (k)*N*N) +=  (std::sin((i + 1) * deltaX) * std::sin((k + 1) * deltaX) * std::sin(1.0))/(deltaX*deltaX); // k == N-1
                }
                if(k == N - 1){
                    rhs(i + (j)*N  + (k)*N*N) +=  (std::sin(1.0) * std::sin((i + 1) * deltaX) * std::sin((j + 1) * deltaX))/(deltaX*deltaX); // j == N-1
                }
            }                        
        }
    }

    A.setFromTriplets(triplets.begin(), triplets.end());

    // std::cout << A << std::endl;

}

void system_solver(SpMat &A, SpVec &rhs, SpVec &solution){

    //Create preconditioner
    // Eigen::IncompleteLUT<double> preconditioner(A);
    // preconditioner.compute(A); 
    // Use Conjugate Gradient solver to solve the linear system

    Eigen::BiCGSTAB<SpMat> solver;

    solver.setTolerance(1e-10);

    solver.setMaxIterations(1000);

    solver.compute(A);

    if (solver.info() != Eigen::Success) {
        // decomposition failed
        std::cout << "Failed computing factorization" << std::endl;
    }

    solution = solver.solve(rhs); // Solve Ax = b

    if (solver.info() != Eigen::Success) {
        // solving failed
        std::cout << "Failed solving" << std::endl;
    }
}

double error_computation(SpVec &solution, SpVec &real_solution, const int N, const double deltaX){
    double error;

    for (int k = 0; k < N; k++){
        for(int j = 0; j < N; j++){
            for(int i = 0; i < N; i ++){
                real_solution(i + (j)*N  + (k)*N*N) = std::sin((i + 1) * deltaX) * std::sin((j + 1) * deltaX) * std::sin((k + 1) * deltaX);
            }
        }
    }    

    return error = (solution - real_solution).norm();
}

int main(int argc, char* argv[]){

    std::vector<double> errors;
    std::vector<double> deltaXs;

    int const  N_vec[5] = {20, 40, 60, 80, 100};
    for (int N:N_vec) {
 
        int const dim = N*N*N;
        double const deltaX = 1.0/(double(N) + 1);

        SpMat A(dim,  dim);
        SpVec rhs(dim); 
        SpVec solution(dim);
        SpVec real_solution(dim);

        std::cout << N << " points" << std::endl;

        auto start_time
            = std::chrono::system_clock::now();

        system_assemble(A, rhs, N, deltaX);

        auto end_time
            = std::chrono::system_clock::now();

        std::chrono::duration<double> matrix_assemble 
            = end_time - start_time;


        std::cout << "Matrix assemble time: "
                    << matrix_assemble.count() << " seconds " << std::endl << std::endl;    

            start_time
            = std::chrono::system_clock::now();

        system_solver(A, rhs, solution);

        end_time
            = std::chrono::system_clock::now();

        std::chrono::duration<double> solving_time 
            = end_time - start_time;


        std::cout << "System solving time: "
                    << solving_time.count() << " seconds " << std::endl << std::endl;    

        double error = error_computation(solution, real_solution, N, deltaX)/dim;

        std::cout << "Error with " << N << " points: " << error << " delta x: " << deltaX << std::endl << std::endl; 

        errors.push_back(error);
        deltaXs.push_back(deltaX);

    }

    for(int j=1; j < 5; j++){
        std::cout << (log(errors.at(j)) - log(errors.at(j - 1)))/(3 * (log(deltaXs.at(j)) - log(deltaXs.at(j - 1)))) << std::endl;
    }
    
}