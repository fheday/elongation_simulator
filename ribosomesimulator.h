#ifndef SIMULATIONS_RIBOSOMESIMULATOR_H
#define SIMULATIONS_RIBOSOMESIMULATOR_H

#include <map>
#include <random>
#include <tuple>
#include <vector>
#include "concentrationsreader.h"
#include "gillespie.h"
#include "reactionsset.h"

namespace Simulations {

class RibosomeSimulator {
 public:
  RibosomeSimulator();
  int getState();
  void setState(int);
  void getAlphas(std::vector<double>&, std::vector<int>&);
  void getDecodingAlphas(std::vector<double>&, std::vector<int>&);
  void loadConcentrations(std::string);
  void setCodonForSimulation(const std::string&);
  void run_and_get_times(double&, double&);
  csv_utils::ConcentrationsReader concentrations_reader;

 private:
  std::vector<std::vector<std::tuple<double, int>>> createReactionsGraph(
      const csv_utils::concentration_entry&);
  std::map<std::string, std::vector<std::vector<std::tuple<double, int>>>>
      reactions_map;
  std::vector<std::vector<std::tuple<double, int>>>
      reactions_graph;  // vector where the index is the ribosome's current
                        // state and the content is a vector of tuples
                        // containing the propensity and next state of each
                        // possible reaction.
  int current_state = 0;
};
}  // namespace Simulations

#endif  // SIMULATIONS_RIBOSOMESIMULATOR_H
