#ifndef ELONGATION_MANAGER_ELONGATIONSIMULATIONMANAGER_H
#define ELONGATION_MANAGER_ELONGATIONSIMULATIONMANAGER_H

#include "translation.h"
#include <future>
#include <string>
#include <tuple>
#include <vector>

namespace Elongation_manager {
enum stop_condition_enum { ITERATION, TIME, RIBOSOMES };
class SimulationManager {
public:
  SimulationManager() = delete;   // no default constructor
  SimulationManager(std::string); // constructor with configuration file name.
  std::string get_concentration_file_path();
  std::string get_configuration_file_path();
  bool get_pre_populate();
  std::vector<std::tuple<std::string, std::string, float, float, float>> &
  get_simulations_configurations();
  stop_condition_enum get_stop_condition_type();
  float get_stop_condition_value();
  std::size_t get_history_size();
  bool start();

private:
  std::string configuration_file_path;
  std::string concentration_file_path;
  std::string directory;
  bool pre_populate = false;
  std::vector<std::tuple<std::string, std::string, float, float, float>>
      simulations_configurations;
  stop_condition_enum stop_condition_type;
  float stop_condition_value;
  std::size_t history_size;
  std::vector<std::future<Simulations::Translation>> simulations;

  bool is_simulation_valid();
  bool save_sim(Simulations::Translation&);
};
} // namespace Elongation_manager
#endif // ELONGATION_MANAGER_ELONGATIONSIMULATIONMANAGER_H