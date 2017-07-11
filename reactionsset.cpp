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
    k_pop_index.push_back(std::make_tuple(k, maxRow, maxCol));
//     std::cout << "Max: " << max <<  ", at: " << maxRow << "," << maxCol;
      
}

Eigen::MatrixXi ReactionsSet::getReaction(int index)
{
    return reactions_vector[index];
}


ReactionsSet::~ReactionsSet()
{
}
