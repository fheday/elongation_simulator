#ifndef REACTIONSSET_H
#define REACTIONSSET_H

#include <eigen3/Eigen/Dense>

namespace Simulations {

class ReactionsSet
{
public:
ReactionsSet();
void addReaction(Eigen::MatrixXi reaction,  float k);
Eigen::MatrixXi getReaction(int);
~ReactionsSet();
    std::vector<Eigen::MatrixXi> reactions_vector;
    Eigen::ArrayXf ks;

private:
    Eigen::ArrayXf k_times_reaction;
};
}
#endif // REACTIONSSET_H
