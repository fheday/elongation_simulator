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
        Gillespie();
        Gillespie(int, const Eigen::MatrixXi&, const ReactionsSet&);
        void run();
        void setReactionsSet(const ReactionsSet&);
        void setInitialPopulation(const Eigen::MatrixXi&);
        void setIterationLimit(int);
        ~Gillespie();
        int iteration_limit;
        Eigen::MatrixXi initial_populations;
        ReactionsSet reactions;
        std::vector<float> dt_history;
        std::vector<Eigen::MatrixXi> population_history;
        float total_time;
        
    };
}

#endif // SIMULATIONS_GILLESPIE_H
