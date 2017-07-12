 #include <iostream>
 #include "gillespie.cpp"
 #include "reactionsset.cpp"
 #include "concentrations_reader.cpp"
//  #include <boost/numeric/ublas/matrix.hpp>
 #include <vector>
//  #include <boost/numeric/ublas/io.hpp>

 
int main(int argc, char **argv) {
    std::cout << "Hello, world!" << std::endl;
    csv_utils::concentrations_reader cr;
    std::vector<csv_utils::concentration_entry> concentrations_vector;
    cr.get_contents(concentrations_vector);
    //print
    for (csv_utils::concentration_entry item: concentrations_vector) {
        std::cout<<"codon: " << item.codon<<", "<<item.three_letter<< ", WC cognate: " << item.wc_cognate_conc << ", wobble = "<< item.wobblecognate_conc<<" nearcognate = "<<item.nearcognate_conc<<"\n";
    }
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
    return 0;
}
