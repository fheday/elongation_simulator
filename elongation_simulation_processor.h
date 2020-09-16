#ifndef ELONGATION_SIMULATION_READER_H
#define ELONGATION_SIMULATION_READER_H

#include <string>
#include <vector>

namespace Simulations {
    class SimulationProcessor {
        public:
        SimulationProcessor() = delete;
        SimulationProcessor(std::string);
        std::vector<float>& getClock();
        std::vector<std::vector<int>>& getRibosomes();
        private:
        std::string fasta_file;
        float initiation_rate;
        float termination_rate;
        std::vector<float> clock;
        std::vector<std::vector<int>> elongating_ribosomes;
    };
}
#endif