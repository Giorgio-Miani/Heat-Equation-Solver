#include <chrono>
#include <cstdlib>
#include <iostream>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <unsupported/Eigen/SparseExtra>

    using SpMat = Eigen::SparseMatrix<double>;
    using SpVec = Eigen::VectorXd;


int main(int argc, char* argv[])
{

    
    int N[7] = {100,1000,10000,100000, 1000000, 10000000};
    for (int n:N) {
    const double deltaX = 1.0/(n+1);

    SpVec f(n);
    SpVec solution(n);
    SpVec b(n);
    double a;

    auto start_time_assembly
        = std::chrono::system_clock::now();

    a = -1/(deltaX*deltaX);

    for (int i = 0; i <= n - 1; i++){
        b(i) = 2/(deltaX*deltaX);
        f(i) = std::sin((i + 1) * deltaX);
    }
    f(n-1) += std::sin(1.0)/(deltaX*deltaX);

    auto start_time_solving
        = std::chrono::system_clock::now();

    double w;
    for(int i = 1; i < n; i++){
        w = a/b(i - 1);
        b(i) -= w*a;
        f(i) -= w*f(i-1);
    }

    solution(n - 1) = f(n - 1)/b(n-1); 
    for(int i = n - 2; i >= 0; i-- ){
        solution(i) = (f(i) - a*solution(i + 1))/b(i); 
    }

    auto end_time
        = std::chrono::system_clock::now();
    
    std::chrono::duration<double> assembly_duration 
        = start_time_solving - start_time_assembly;

    std::chrono::duration<double> solving_duration 
        = end_time - start_time_solving;
    //std::cout << solver.error();

   SpVec real_solution(n);
    for(int i=0; i < n; i++){
        real_solution(i) = std::sin((i + 1) * deltaX); 
    }
    std::cout << "Error with " << n << " points: " << ((solution - real_solution)).norm()/n << " Delta x: " << deltaX << std::endl; 

    std::cout << "Assembly duration: "
              << assembly_duration.count()*1000 << " milliseconds " << std::endl;

    std::cout << "Solving duration: "
              << solving_duration.count()*1000 << " milliseconds " << std::endl;


    std::cout << "Total time: "
              << (assembly_duration.count() + solving_duration.count())/n << " seconds/points " << std::endl << std::endl;
    }


}