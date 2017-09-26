#ifndef SIMULATIONS_RIBOSOMESIMULATOR_H
#define SIMULATIONS_RIBOSOMESIMULATOR_H

#include <vector>
#include <map>
#include <tuple>
#include <random>
#include "concentrationsreader.h"
#include "reactionsset.h"
#include "gillespie.h"

namespace Simulations {
    
    class RibosomeSimulator : public Gillespie
    {
    public:
        RibosomeSimulator();
        int getState();
        void setState(int);
        void getAlphas(std::vector<double>&, std::vector<int>&);
        void loadConcentrations(std::string);
        void setCodonForSimulation(const std::string&);
        void run_and_get_times(double&, double&);
        void getDecodingAndTranslocationTimes(double&, double&);
        csv_utils::ConcentrationsReader concentrations_reader;
        std::mt19937 gen;
        std::uniform_real_distribution<> dis;
    private:
        std::vector<std::vector<std::tuple<double, int>>> createReactionsGraph(const csv_utils::concentration_entry&);
        std::map<std::string, std::vector<std::vector<std::tuple<double, int>>>> reactions_map;
        std::vector<std::vector<std::tuple<double, int>>> reactions_graph; //vector where the index is the ribosome's current state and the content is a vector of tuples containing the propensity and next state of each possible reaction.
    };
}

#endif // SIMULATIONS_RIBOSOMESIMULATOR_H
