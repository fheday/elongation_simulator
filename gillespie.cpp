#include "gillespie.h"
#include <vector>
#include <random>
#include <math.h>
#include <boost/algorithm/cxx11/any_of.hpp>
#include <boost/multi_array.hpp>
#include <eigen3/Eigen/Dense>
 using namespace Simulations;

Gillespie::Gillespie(int itera, const Eigen::MatrixXi& popul, const Eigen::MatrixXi& reac, const Eigen::ArrayXf& _ks)
{
    iteration_limit = itera;
    dt_history.clear();
    population_history.clear();
    initial_populations = popul;
    reactions = reac;
    ks = _ks;
}

void Gillespie::run()
{
    float dt = 0, t = 0;
    Eigen::ArrayXi populations = initial_populations;

    // initialize the random generator
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(1, 2);
    float r1, r2, a0, tau;
//     bool lessThan0 ( int i ) { return i < 0; };
    for (int i =0; i < iteration_limit; i++)
    {
        //population_history->push_back(populations);
        dt_history.push_back(dt);
        // randomly generate parameter for calculating dt
        r1 = dis(gen);
        // randomly generate parameter for selecting reaction
        r2 = dis(gen);
        
        // calculate an
        Eigen::ArrayXf as = Eigen::ArrayXf();
        Simulations::get_an(reactions, populations, ks, as);
        if (a0 == 0)
        {
            // no available reactions, quit loop prematurely.
            break;
        }
        tau = 1/a0 * log(1/r1);  // calculate time of next reaction
        // select next reaction to execute
//         select_next_reaction(r2);
        // put the new population on a temporary variable: this is to help
        // detecting negative populations.
//         updated_populations = populations + reactions.get_selected_reaction_population_change();
        if ((populations == 0).any())
        {
            break;
        }  
        else
        {
            // update time / clock
            t += tau;
            // update population
//             populations = updated_populations;
            
        }
        
    // finished. Plot population snapshots
//     plt.step(clock_times, population_size)
//     plt.show()
    // return history of population size and clock times
//     return population_size, clock_times

    }

}


void Simulations::get_an(const Eigen::MatrixXi& reac, const Eigen::MatrixXi& pop, const Eigen::ArrayXf& ks, Eigen::ArrayXf& as)
{
        as = 0;
        as.resize(1, pop.cols());
        Eigen::MatrixXi indices = (reac.array() < 0).cast<int>();
        std::cout<<"indices : "<< indices<<" Done!";
        std::cout<<"population : "<< pop<<" Done!";
        std::cout<<"reactions : "<< reac<<" Done!";
        std::cout<<"ks : "<< ks<<" Done!\n";
//         for ()
        
//         for (int i = 0; i < reac.size(); i++){
//             
//             if (reac(i) < 0 ) {
//                 // compute an.
//                 for (int j = 0 
//                 
//             }
//             else {
//                 as(i) = 0;
//             }
//         }
        
        
//         alphas = pop * ks;
//         return alphas.sum();

}

Gillespie::~Gillespie()
{
    
}
