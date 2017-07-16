 #include <iostream>
 #include <iomanip>
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
    cr.load_concentrations("../../../Projects/RSim/data_with_times/concentrations.csv");
    Simulations::RibosomeSimulator rs("../../../Projects/RSim/data_with_times/concentrations.csv");
    rs.setIterationLimit(1000);
    // create a matrix with the initial species population.
    rs.setNumberOfRibosomes(1);
    
    float decoding, translocating;
    std::vector<std::string> codons;
    cr.get_codons_vector(codons);
     std::cout<<std::setprecision(10);
    for (std::string codon:codons){
        rs.setCodonForSimulation(codon);
        std::cout<<"+++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n"<<"                      "<<codon<<"                        \n"<<"+++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
        for (int i = 0 ; i < 10000; i++){
            rs.run_and_get_times(decoding, translocating);
            std::cout<<" decoding time = "<< decoding << ", translocating time = " << translocating<<", Total time = "<<(decoding + translocating)<<"\n";
        }
    }

    return 0;
}
