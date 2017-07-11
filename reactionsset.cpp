#include "reactionsset.h"
using namespace Simulations;

ReactionsSet::ReactionsSet()
{
    reactions_vector.clear();
    ks.clear();
}

void ReactionsSet::addReaction(Eigen::MatrixXi reaction,  float k)
{
    reactions_vector.push_back(reaction);
    ks.push_back(k);
}

Eigen::MatrixXi ReactionsSet::getReaction(int index)
{
    return reactions_vector[index];
}

ReactionsSet::~ReactionsSet()
{
}
