#ifndef REACTIONSSET_H
#define REACTIONSSET_H

#include <eigen3/Eigen/Dense>
#include <vector>

namespace Simulations {

class ReactionsSet
{
public:
ReactionsSet();
void addReaction(Eigen::MatrixXi reaction,  float k);
Eigen::MatrixXi getReaction(int);
~ReactionsSet();
    std::vector<Eigen::MatrixXi> reactions_vector;
    std::vector<float> ks;

private:
    std::vector<std::tuple<float, int, int>> k_pop_index;
};
}
#endif // REACTIONSSET_H
