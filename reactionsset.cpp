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
    addReaction(reaction, k);
    decrptions.push_back(reaction_id); //reaction's description. e.g.: in a mRNA it would be the codon.
}

void ReactionsSet::addReaction(Eigen::MatrixXi reaction,  double k)
{
    reactions_vector.push_back(reaction); // add reaction to the list
    ks.push_back(k); // add propensity coefficient to the list.
    // now we get the index of the consumed species from the reaction.
    Eigen::MatrixXi indexes = (reaction.array() < 0).cast<int>();
    //get location of maximum (we expect only one value=1, others zero.
    Eigen::MatrixXd::Index maxRow, maxCol;
    double max = indexes.maxCoeff(&maxRow, &maxCol);
    //add the tuple k, row, col to the k_pop_index vector.
    //this will help us to calculate the alphas later.
    if (max == 0){ 
        //zero order reactions will have a negative col, row as marker
        maxRow = -1;
        maxCol = -1;
    }
    k_pop_index.push_back(std::make_tuple(k, maxRow, maxCol));
}

void ReactionsSet::getAlphas(const Eigen::MatrixXi& species, Eigen::VectorXd& as_vector, Eigen::VectorXi& reaction_number_vector)
{
    std::vector<double> as;
    std::vector<int> reaction_number;
    for (unsigned int i = 0; i < k_pop_index.size(); i++){
        auto k_and_index = k_pop_index.at(i); //get the first element.
        if (std::get<1>(k_and_index) < 0){
            //zero order reaction. no need for species.
            as.push_back(std::get<0>(k_and_index));
            reaction_number.push_back(i); //save index number.
        } else {
            // first order reaction. dependent on species.
            int specie_population = species(std::get<1>(k_and_index), std::get<2>(k_and_index));
            if (specie_population != 0){
                as.push_back(std::get<0>(k_and_index) * specie_population);
                reaction_number.push_back(i); //save index number.
            }
        }
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
