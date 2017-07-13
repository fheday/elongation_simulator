 #include <iostream>
 #include "gillespie.cpp"
 #include "reactionsset.cpp"
 #include "concentrations_reader.cpp"
 #include <vector>
 #include "ribosomesimulator.cpp"


 void testGillespie()
 {
    // create one reaction matrix
    Eigen::MatrixXi reactions(4,10);
    reactions.fill(0);
    reactions(1,1) = -1;
    reactions(2,1) = 1;
    // create a ReactionSet object and add the reaction there.
    Simulations::ReactionsSet reactions_set;
    // the reaction propensity is the second argument.
    reactions_set.addReaction(reactions, .1);
    // create a matrix with the initial species population.
    Eigen::MatrixXi population(4,10);
    population.fill(0);
    population(1,1) = 800;
    std::cout<<"initial population = "<<population<<"\n------\n";
    // create the Gillespie simulator object with the population and reactions
    Simulations::Gillespie simulation(300, population, reactions_set);
    //run the simulation.
    simulation.run();
 }
 
int main(int argc, char **argv) {
    std::cout << "Hello, world!" << std::endl;
    csv_utils::concentrations_reader cr;
    std::vector<csv_utils::concentration_entry> concentrations_vector;
    Simulations::RibosomeSimulator rs(cr);
//         for (std::pair<std::string, ReactionsSet> element:rs.reactions_map){
//         std::cout<<"Starting new reactions set: "<< element.first <<"\n ++++++++++++\n";
//         for (Eigen::MatrixXi m: element.second.reactions_vector) {
//             std::cout<<"Reaction Matrix = \n"<< m<<"\n---------\n";
//         }
//         std::cout<<"Finished reactions set\n ++++++++++++\n";
//     }
    rs.setIterationLimit(1000);
    // create a matrix with the initial species population.
    Eigen::MatrixXi population(32, 1);
    population.fill(0);
    population(0,0) = 1;
    rs.setInitialPopulation(population);
    std::string codon = "AAA";
    rs.setCodonForSimulation(codon);
    float decoding, translocating;
    rs.run_and_get_times(decoding, translocating);
    std::cout<<" decoding time = "<< decoding << ", translocating time = " << translocating<<", Total time = "<<rs.total_time<<"\n";
//     for (float dt:rs.dt_history) {
//         std::cout<<"dt = "<<dt<<"\n";
//     }
    return 0;
}
