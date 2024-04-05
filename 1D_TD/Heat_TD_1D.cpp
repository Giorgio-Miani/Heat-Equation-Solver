#include <chrono>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <unsupported/Eigen/SparseExtra>

    using SpMat = Eigen::SparseMatrix<double>;
    using SpVec = Eigen::VectorXd;


int main(int argc, char* argv[])
{

    std::vector<double> errors;
    std::vector<double> deltaXs;

    std::ofstream fout("data.csv");
    fout << "deltaX,"  << "err" << "\n";

    int N[5] = {11, 21, 41, 81, 161};
    for (int n:N) {
        const double deltaX = 1.0/(n-1);
        const double deltaT = 0.25*deltaX*deltaX;
        double Nt = 1/deltaT;
        double t = 0.0;
        double Nx = 1/deltaX;
        SpVec solution(n);
        SpVec new_solution(n);
        SpVec real_solution(n);
        SpVec force(n);
        double max_error = 0.0;
        double error = 0.0;
        double c = (deltaT/(deltaX*deltaX));


        //Initialization
        for(int i=0; i < n; i++){
            solution(i) = std::sin(i * deltaX)*std::sin(0); 
            force(i) = deltaT*std::sin(i * deltaX);
        }

        auto start_time_solving
            = std::chrono::system_clock::now();

        t += deltaT;
        while(t <= 1.0){
            double f = (std::cos(t) + std::sin(t));
            double old_value = solution(0);
            for (int i = 1; i < n - 1; i++){
                double new_value = solution(i) * (1 - 2*c) + c*(old_value + solution(i+1)) + force(i)*f;
                old_value = solution(i);
                solution(i) = new_value;
            }

            solution(0) = std::sin(0)*std::sin(t);
            solution(n - 1) = std::sin(1) * std::sin(t);
            t += deltaT;

            // for(int i=0; i < n; i++){
            //     real_solution(i) = std::sin(i * deltaX)*std::sin(t); 
            // }

            //Calculate Error
            //error = 0.0;
            // for(int i=1; i < n - 1; i++)
            //         error += (solution(i) - real_solution(i)) * (solution(i) - real_solution(i));
            // error = std::sqrt(error);
            // error = error / n;

            //Save only the max error for every temporal step
            // if(error > max_error){
            //     max_error = error;
            // }
            
        }

        auto end_time
            = std::chrono::system_clock::now();

        std::chrono::duration<double> solving_duration 
            = end_time - start_time_solving;

        for(int i=0; i < n; i++){
            real_solution(i) = std::sin(i * deltaX)*std::sin(t); 
        }

        //Calculate Error
        error = 0.0;
        for(int i=1; i < n - 1; i++)
                error += (solution(i) - real_solution(i)) * (solution(i) - real_solution(i));
        error = std::sqrt(error);
        max_error = error / n;


        std::cout << "Nx: " << Nx << " Nt: " << Nt << std::endl;
        std::cout << "Error with " << n << " points: " << max_error << " delta x: " << deltaX << std::endl; 
        std::cout << "Total time: "
                    << (solving_duration.count())/(Nx*Nt) << " seconds/points " << std::endl << std::endl;

        fout << deltaX << "," << max_error << "\n";

        errors.push_back(max_error);
        deltaXs.push_back(deltaX);
    }

    fout.close();

    for(int j=1; j < 5; j++){
        std::cout << (log(errors.at(j)) - log(errors.at(j - 1)))/(log(deltaXs.at(j)) - log(deltaXs.at(j - 1))) << std::endl;
    }
}