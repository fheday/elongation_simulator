 #include <iostream>
 #include "gillespie.cpp"
 #include "reactionsset.cpp"
 #include <boost/numeric/ublas/matrix.hpp>
 #include <boost/numeric/ublas/vector.hpp>
 #include <boost/numeric/ublas/io.hpp>

 
int main(int argc, char **argv) {
    std::cout << "Hello, world!" << std::endl;
//     Eigen::MatrixXi initial_population(1,1);
//     initial_population(0) = 20;
//     Eigen::MatrixXi reactions(1,1);
//     reactions(0) = -1;
//     Eigen::ArrayXf ks(1);
//     ks(0)=0.1;
//     Simulations::Gillespie simulation(300, initial_population, reactions, ks);
//     simulation.run();
    Eigen::MatrixXi reactions(4,10);
    reactions(1,1) = -1;
    reactions(2,1) = -1;
    Simulations::ReactionsSet reactions_set;
    reactions_set.addReaction(reactions, 0.1);

    return 0;
}
