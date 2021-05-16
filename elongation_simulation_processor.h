#ifndef ELONGATION_SIMULATION_READER_H
#define ELONGATION_SIMULATION_READER_H

/*
 * @file  elongation_simulation_processor.h
 * 
 * @brief manager class to process parallel simulations output and save its results
 *
 * @author Fabio Hedayioglu
 * Contact: fheday@gmail.com
 *
 */

#include <string>
#include <vector>
#include "json/json.h"

namespace Simulations {
    class SimulationProcessor {
        public:
        SimulationProcessor() = delete;
        SimulationProcessor(std::string);
        SimulationProcessor(Json::Value, std::string);
        std::vector<float>& getClock();
        std::vector<std::vector<int>>& getRibosomes();
        void removeRibosomePositions();
        void calculateRibosomeCollisions();
        std::vector<std::vector<int>>& getCollidingRibosomes();
        std::vector<std::vector<int>>& getStalledRibosomes();
        void packData();
        void save();

        private:
        std::string configuration_file_name;
        std::string fasta_file;
        float initiation_rate;
        float termination_rate;
        std::vector<float> clock;
        std::vector<std::vector<int>> elongating_ribosomes;
        std::vector<std::vector<int>> colliding_ribosomes;
        std::vector<std::vector<int>> stalled_ribosomes;
        void parseJson(Json::Value&);
    };
}
#endif