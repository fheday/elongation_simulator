#ifndef ELONGATION_SIMULATION_MANAGER_H
#define ELONGATION_SIMULATION_MANAGER_H

/*
 * @file  elongation_simulation_manager.h
 * 
 * @brief manager class to load and execute parallel simulations
 *
 * @author Fabio Hedayioglu
 * Contact: fheday@gmail.com
 *
 */

#include "translation.h"
#include <future>
#include <string>
#include <tuple>
#include <vector>

namespace Elongation_manager {
enum stop_condition_enum { ITERATION, TIME, RIBOSOMES, STEADY_STATE_TIME, STEADY_STATE_RIBOSOMES };
class SimulationManager {
public:
  SimulationManager() = delete;   // no default constructor
  SimulationManager(std::string); // constructor with configuration file name.
  std::string get_concentration_file_path();
  std::string get_configuration_file_path();
  std::map<std::string, float> get_reactions_modifiers();
  bool get_pre_populate();
  std::vector<std::tuple<std::string, std::string, float, float, float>> &
  get_simulations_configurations();
  stop_condition_enum get_stop_condition_type();
  float get_stop_condition_value();
  std::size_t get_history_size();
  bool start(bool verbose=false, unsigned int n_threads = std::thread::hardware_concurrency());
  void set_save_collisions(bool);
  void set_remove_ribosome_positions(bool);

private:
  std::string configuration_file_path;
  std::string concentration_file_path;
  std::string directory;
  bool pre_populate = false;
  bool save_collisions = false;
  bool remove_ribosome_positions = false;
  std::vector<std::tuple<std::string, std::string, float, float, float>>
      simulations_configurations;
  stop_condition_enum stop_condition_type;
  float stop_condition_value;
  std::size_t history_size;
  std::vector<std::future<Simulations::Translation>> simulations;
  std::map<std::string, float> reactions_modifiers;
    std::array<std::string, 44> reactions_identifiers = {
      {"non1f",    "near1f",     "wobble1f", "WC1f",     "non1r",    "near1r",
       "near2f",   "near2r",     "near3f",   "near4f",   "near5f",   "neardiss",
       "near6f",   "wobble1r",   "wobble2f", "wobble2r", "wobble3f", "wobble4f",
       "wobble5f", "wobblediss", "wobble6f", "WC1r",     "WC2f",     "WC2r",
       "WC3f",     "WC4f",       "WC5f",     "WCdiss",   "WC6f",     "dec7f",
       "trans1f",  "trans1r",    "trans2",   "trans3",   "trans4",   "trans5",
       "trans6",   "trans7",     "trans8",   "trans9"}};

  bool is_simulation_valid();
  bool save_sim(Simulations::Translation&);
};
} // namespace Elongation_manager
#endif // ELONGATION_MANAGER_ELONGATIONSIMULATIONMANAGER_H