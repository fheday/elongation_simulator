#ifndef SIMULATIONS_RIBOSOMESIMULATOR_H
#define SIMULATIONS_RIBOSOMESIMULATOR_H

#include <vector>
#include <map>
#include "concentrations_reader.h"
#include "reactionsset.h"
#include "gillespie.h"

namespace Simulations {
    
    class RibosomeSimulator : public Gillespie
    {
    public:
        RibosomeSimulator(csv_utils::concentrations_reader&);
        void setCodonForSimulation(const std::string&);
        void run_and_get_times(float&, float&);
        void getDecodingAndTranslocationTimes(float&, float&);
        csv_utils::concentrations_reader concentrations_reader;
        std::map<std::string, ReactionsSet> reactions_map;
    private:
        ReactionsSet createReactionSet(const csv_utils::concentration_entry&);
    };
}

#endif // SIMULATIONS_RIBOSOMESIMULATOR_H
