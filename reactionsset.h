#ifndef REACTIONSSET_H
# define REACTIONSSET_H

# include <eigen3/Eigen/Dense>
# include <vector>

namespace Simulations {

class ReactionsSet
{
public:
    ReactionsSet();
    void addReaction(Eigen::MatrixXi,  double);
    void addReaction(Eigen::MatrixXi,  double, std::string);
    Eigen::MatrixXi getReaction(int);
    void getAlphas(const Eigen::MatrixXi&, Eigen::VectorXd&,  Eigen::VectorXi&);
    ~ReactionsSet();
    std::vector<Eigen::MatrixXi> reactions_vector;
    std::vector<double> ks;
    std::vector<std::string> decrptions;

private:
    std::vector<std::tuple<double, int, int>> k_pop_index;
};
}
#endif                                                      // REACTIONSSET_H
