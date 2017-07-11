#ifndef SIMULATIONS_GILLESPIE_H
#define SIMULATIONS_GILLESPIE_H

#include <memory>
#include <list>
#include <eigen3/Eigen/Dense>

namespace Simulations {

class Gillespie
{
public:
Gillespie(int, const Eigen::MatrixXi&, const Eigen::MatrixXi&, const Eigen::ArrayXf&);
void run();
~Gillespie();
    int iteration_limit;
    Eigen::MatrixXi initial_populations;
    Eigen::MatrixXi reactions;
    Eigen::ArrayXf ks;
    std::list<float> dt_history;
    std::list<Eigen::MatrixXi> population_history;

};
    void get_an(const Eigen::MatrixXi&, const Eigen::MatrixXi&, const Eigen::ArrayXf&, Eigen::ArrayXf&);
}

#endif // SIMULATIONS_GILLESPIE_H
