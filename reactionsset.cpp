#include "reactionsset.h"
using namespace Simulations;

ReactionsSet::ReactionsSet()
{
    reactions_vector.clear(); //remove all reactions from vector.
    ks.clear(); // remove all propensities from vector.
    k_pop_index.clear(); // also remove the pre-calculated values.
}

void ReactionsSet::addReaction(Eigen::MatrixXi reaction,  double k, std::string reaction_id)
{
    addReaction(reaction, k); // add the reaction.
    decrptions.push_back(reaction_id); //add reaction's description. e.g.: in a mRNA it would be the codon.
}

void ReactionsSet::addReaction(Eigen::MatrixXi reaction,  double k)
{
    reactions_vector.push_back(reaction); // add reaction to the list
    ks.push_back(k); // add propensity coefficient to the list.
    // now we get the index of the consumed species from the reaction.
    std::vector<std::tuple<int,int>> reactants;
    for (int i = 0; i < reaction.cols(); i++){
        for (int j = 0; j < reaction.rows();  j++) {
            if (reaction(j, i) < 0 ){
                reactants.push_back(std::make_tuple(j, i));
            }
        }
    }
    k_pop_index.push_back(std::make_tuple(k, reactants));
}

void ReactionsSet::getAlphas(const Eigen::MatrixXi& species, Eigen::VectorXd& as_vector, Eigen::VectorXi& reaction_number_vector)
{
    std::vector<double> as;
    std::vector<int> reaction_number;
    int i = 0;
    for (std::tuple<double, std::vector<std::tuple<int, int>>> k_and_index : k_pop_index){
        double k = std::get<0>(k_and_index);
        if (std::get<1>(k_and_index).empty() && k > 0) {
            //zero order reaction. no need for species.
            as.push_back(std::get<0>(k_and_index)); // a = propensity.
            reaction_number.push_back(i);
        } else {
            // first order reaction. dependent on species.
            int specie_population = 1;
            for (std::tuple<int, int> element:std::get<1>(k_and_index)) specie_population *= species(std::get<0>(element), std::get<1>(element));
            // store the reaction's propensity and the locations of the reactants in the species matrix.
            if (k!=0 && specie_population != 0){
                as.push_back(k * specie_population);
                reaction_number.push_back(i); //save index number.
            }
        }
        i++;
    }
    as_vector = Eigen::VectorXd::Map(as.data(), as.size());
    reaction_number_vector = Eigen::VectorXi::Map(reaction_number.data(), reaction_number.size());
}


Eigen::MatrixXi ReactionsSet::getReaction(int index)
{
    return reactions_vector.at(index);
}


ReactionsSet::~ReactionsSet()
{
}
