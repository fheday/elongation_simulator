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
        virtual double getReactionTime(double, double, std::string);
        ~Gillespie();
        int iteration_limit;
        Eigen::MatrixXi initial_populations;
        ReactionsSet reactions;
        std::vector<double> dt_history;
        std::vector<Eigen::MatrixXi> population_history;
        double total_time;
        
    };
}

#endif // SIMULATIONS_GILLESPIE_H
