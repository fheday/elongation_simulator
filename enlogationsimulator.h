#ifndef SIMULATIONS_ENLOGATIONSIMULATOR_H
#define SIMULATIONS_ENLOGATIONSIMULATOR_H

#include "enlogationsimulator.h"
#include "mrna_reader.h"
#include "ribosomesimulator.h"

namespace Simulations {

class EnlogationSimulator {
 public:
  EnlogationSimulator();
  void setTerminationRate(double);
  void setInitiationRate(double);
  void setMRnaFileName(const std::string&);
  void setConcentrationsFileName(const std::string&);
  void setAverageTimesFileName(const std::string&);
  double getReactionTime(double&, double&, std::string&);
  std::tuple<std::vector<double>, std::vector<int>, std::vector<int>>
  getEnlogationDuration();
  void calculateAverageTimes();
  double termination_rate = -1;
  double initiation_rate = -1;
  std::string mRNA_file_name;
  std::string concentrations_file_name;
  std::string average_times_file_name;
  mRNA_utils::mRNAReader mrna_reader;
  RibosomeSimulator ribosome_simulator;
  std::vector<std::vector<int>> ribosome_positions_history;
  std::vector<double> dt_history;
  // array with the total times the ribosomes spent in the codons
  std::vector<double> total_time;
  // number of times a codon was occupied
  std::vector<int> n_times_occupied;
  // average occupation time
  std::vector<double> codons_average_occupation_time;

 private:
  void intializeMRNAReader();
  std::vector<double> translocation_times;
  std::vector<std::string> stop_codons = {"UAG", "UAA",
                                          "UGA"};  // list of stop codons.
};
}  // namespace Simulations

#endif  // SIMULATIONS_ENLOGATIONSIMULATOR_H
