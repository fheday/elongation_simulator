#ifndef SIMULATIONS_RIBOSOMESIMULATOR_H
#define SIMULATIONS_RIBOSOMESIMULATOR_H

#include <vector>
#include <map>
#include <tuple>
#include "concentrationsreader.h"
#include "reactionsset.h"
#include "gillespie.h"

namespace Simulations {
    
    class RibosomeSimulator : public Gillespie
    {
    public:
        RibosomeSimulator();
        void loadConcentrations(std::string);
        void setNumberOfRibosomes(int);
        void setCodonForSimulation(const std::string&);
        void run_and_get_times(double&, double&);
        void getDecodingAndTranslocationTimes(double&, double&);
        csv_utils::ConcentrationsReader concentrations_reader;
        std::map<std::string, ReactionsSet> reactions_map;
    private:
        ReactionsSet createReactionSet(const csv_utils::concentration_entry&);
    };
}

#endif // SIMULATIONS_RIBOSOMESIMULATOR_H