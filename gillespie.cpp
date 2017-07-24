#include "gillespie.h"
#include <vector>
#include <random>
#include <math.h>
#include <eigen3/Eigen/Dense>
#include <limits>
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

void Gillespie::setIterationLimit(int i)
{
    iteration_limit = i;
}


void Gillespie::setReactionsSet(const ReactionsSet& reac)
{
    reactions = reac;
}

void Gillespie::run()
{
    dt_history.clear();
    population_history.clear();
    Eigen::MatrixXi populations = initial_populations;
    Eigen::MatrixXi updated_populations;
    Eigen::VectorXd as;
    Eigen::VectorXi reactions_index;
    // initialize the random generator
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(0, 1);
     
    double r1 = 0, r2 = 0;
    double tau = 0, clock = 0.0;
    for (int i = 0; i < iteration_limit; i++)
    {
        population_history.push_back(populations);
        dt_history.push_back(tau);
        // randomly generate parameter for calculating dt
        r1 = dis(gen);//dis(gen);
        // randomly generate parameter for selecting reaction
        r2 = dis(gen);//dis(gen);
        // calculate an
        reactions.getAlphas(populations, as, reactions_index);
        if (as.size() == 0)
        {
            // no available reactions, quit loop prematurely.
            break;
        }
        double a0 = as.sum();
        // select next reaction to execute
        double cumsum = 0;
        int selected_index = -1;
        do {
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
            std::string codon_name = reactions.decrptions[selected_index];
            tau = getReactionTime(a0, r1, codon_name);// calculate time of next reaction
            // update time / clock
            clock += tau;
            // update population
            populations = updated_populations;
        }

    }
    // finished. Plot population snapshots
    total_time = clock;
}

double Gillespie::getReactionTime(double a0, double r1, std::string codon)
{
    return (1.0/a0) * log(1.0/r1);  // calculate time of next reaction
}


Gillespie::~Gillespie()
{
    
}
