#include "reactionsset.h"
using namespace Simulations;

ReactionsSet::ReactionsSet()
{
    reactions_vector.clear(); //remove all reactions from vector.
    ks.clear(); // remove all propensities from vector.
    k_pop_index.clear(); // also remove the pre-calculated values.
}

void ReactionsSet::addReaction(Eigen::MatrixXi reaction,  float k)
{
    reactions_vector.push_back(reaction); // add reaction to the list
    ks.push_back(k); // add propensity coefficient to the list.
    // now we get the index of the consumed species from the reaction.
    Eigen::MatrixXi indexes = (reaction.array() < 0).cast<int>();
    //get location of maximum (we expect only one value=1, others zero.
    Eigen::MatrixXf::Index maxRow, maxCol;
    float max = indexes.maxCoeff(&maxRow, &maxCol);
    //add the tuple k, row, col to the k_pop_index vector.
    //this will help us to calculate the alphas later.
    if (max == 0){ 
        //zero order reactions will have a negative col, row as marker
        maxRow = -1;
        maxCol = -1;
    }
    k_pop_index.push_back(std::make_tuple(k, maxRow, maxCol));
//     std::cout << "Max: " << max <<  ", at: " << maxRow << "," << maxCol;
}

void ReactionsSet::getAlphas(const Eigen::MatrixXi& species, Eigen::VectorXf& as_vector, Eigen::VectorXi& reaction_number_vector)
{
    std::vector<float> as;
    std::vector<int> reaction_number;
    for (int i = 0; i < k_pop_index.size(); i++){
        auto k_and_index = k_pop_index[i]; //get the first element.
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
    as_vector = Eigen::VectorXf::Map(as.data(), as.size());
    reaction_number_vector = Eigen::VectorXi::Map(reaction_number.data(), reaction_number.size());    
    std::cout<< " as = "<< as_vector << "\nreactions number = "<< reaction_number_vector;
//     for ( int i =0; i < as.size(); i++){
//         std::cout<< " as = "<< as[i] << "\nreactions number = "<<reaction_number[i];
//     }
}


Eigen::MatrixXi ReactionsSet::getReaction(int index)
{
    return reactions_vector[index];
}


ReactionsSet::~ReactionsSet()
{
}
