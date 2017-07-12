#include "gillespie.h"
#include <vector>
#include <random>
#include <math.h>
#include <boost/algorithm/cxx11/any_of.hpp>
#include <boost/multi_array.hpp>
#include <eigen3/Eigen/Dense>
using namespace Simulations;

Gillespie::Gillespie()
{
    iteration_limit = 0;
    dt_history.clear();
    population_history.clear();
}

Gillespie::Gillespie(int itera, const Eigen::MatrixXi& popul, const ReactionsSet& reac)
{
    iteration_limit = itera;
    dt_history.clear();
    population_history.clear();
    initial_populations = popul;
    reactions = reac;
}

void Gillespie::setInitialPopulation(const Eigen::MatrixXi& popul)
{
    initial_populations = popul;
}

void Gillespie::setReactionsSet(const ReactionsSet& reac)
{
    reactions = reac;
}

void Gillespie::run()
{
    std::cout<<"Starting simulation... iteration_limit = "<<iteration_limit<<"\n";
    Eigen::MatrixXi populations = initial_populations;
    Eigen::MatrixXi updated_populations;
    // initialize the random generator
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(0, 1);
    float r1, r2, a0, t = 0, tau = 0;
    for (int i =0; i < iteration_limit; i++)
    {
//         std::cout<<"iteration = "<<i<<"\n";
        population_history.push_back(populations);
        dt_history.push_back(tau);
        // randomly generate parameter for calculating dt
        r1 = dis(gen);
        // randomly generate parameter for selecting reaction
        r2 = dis(gen);
        // calculate an
        Eigen::VectorXf as;
        Eigen::VectorXi reactions_index;
        reactions.getAlphas(populations, as, reactions_index);
        
        if (!as.size())
        {
            // no available reactions, quit loop prematurely.
            break;
        }
        float a0 = as.sum();
        tau = 1.0/a0 * log(1.0/r1);  // calculate time of next reaction
        // select next reaction to execute
        float cumsum = 0;
        int selected_index = -1;
        do {
//             std::cout<<"selecting reaction. cumsum = "<<cumsum <<", selected_index = "<<selected_index<<", a0*r2 = "<<(a0*r2)<<", r2 = "<<r2<<", a0= "<<a0<<"\n";
            selected_index++;
            cumsum += as[selected_index]; 
        } while (cumsum < a0 * r2);
        selected_index = reactions_index[selected_index]; //index of selected reaction.
        // put the new population on a temporary variable: this is to help
        // detecting negative populations.
        Eigen::MatrixXi reac = reactions.getReaction(selected_index);
        updated_populations = populations + reac;
        if ((updated_populations.array() < 0).any())
        {
            break;
        }  
        else
        {
            // update time / clock
            t += tau;
            // update population
            populations = updated_populations;
        }

    }
    // finished. Plot population snapshots
    std::cout<<"Finished simulation.\n";
//     for (int i = 0; i < population_history.size(); i++){
//         std::cout<<"population = "<<population_history[i]<<"    delta_time = "<< dt_history[i] <<"\n";
//     }
//     std::cout<<"Total time: "<<t;
}


Gillespie::~Gillespie()
{
    
}
