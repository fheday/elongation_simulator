 #include <iostream>
 #include "gillespie.cpp"
 #include "reactionsset.h"
 #include <boost/numeric/ublas/matrix.hpp>
 #include <boost/numeric/ublas/vector.hpp>
 #include <boost/numeric/ublas/io.hpp>

 
int main(int argc, char **argv) {
    std::cout << "Hello, world!" << std::endl;
    Eigen::MatrixXi initial_population(1,1);
    initial_population(0) = 20;
    Eigen::MatrixXi reactions(1,1);
    reactions(0) = -1;
    Eigen::ArrayXf ks(1);
    ks(0)=0.1;
//     std::shared_ptr<boost::numeric::ublas::matrix<int>> reactions (new boost::numeric::ublas::matrix<int>(1, 1));
//     (*reactions)(0, 0) = -1;
//     std::shared_ptr<boost::numeric::ublas::vector<float>> ks (new boost::numeric::ublas::vector<float>(1));
//     (*ks)(0) = 0.1;
//     Eigen::ArrayXXi population;
//     (*population)(0, 0) = 20; // initial population
//     std::unique_ptr<Simulations::Gillespie> gillespie_simulator(new Simulations::Gillespie(300, population, reactions, ks));
    Simulations::Gillespie simulation(300, initial_population, reactions, ks);
    simulation.run();

    return 0;
}
