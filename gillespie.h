#ifndef SIMULATIONS_GILLESPIE_H
#define SIMULATIONS_GILLESPIE_H

#include <memory>
#include <list>
#include <eigen3/Eigen/Dense>
#include "reactionsset.h"

namespace Simulations {

class Gillespie
{
public:
Gillespie(int, const Eigen::MatrixXi&, const ReactionsSet&);
void run();
~Gillespie();
    int iteration_limit;
    Eigen::MatrixXi initial_populations;
    ReactionsSet reactions;
    std::list<float> dt_history;
    std::list<Eigen::MatrixXi> population_history;

};
}

#endif // SIMULATIONS_GILLESPIE_H
