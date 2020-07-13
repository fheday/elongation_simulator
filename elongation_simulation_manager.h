#ifndef ELONGATION_MANAGER_ELONGATIONSIMULATIONMANAGER_H
#define ELONGATION_MANAGER_ELONGATIONSIMULATIONMANAGER_H

#include <string>
#include <vector>
#include <tuple>

namespace Elongation_manager {
    enum stop_condition_enum {ITERATION, TIME, RIBOSOMES};
    class SimulationManager{
        public:
        SimulationManager() = delete; // no default constructor
        SimulationManager(std::string); //constructor with configuration file name.
        private:
        std::string configuration_file_path;
        std::string concentration_file_path;
        bool pre_populate = false;
        std::vector<std::tuple<std::string, std::string, float, float, float>> simulations_configurations;
        stop_condition_enum stop_condition_type;
        float stop_condition_value;
        std::size_t history_size;
        
    };
}
#endif  // ELONGATION_MANAGER_ELONGATIONSIMULATIONMANAGER_H