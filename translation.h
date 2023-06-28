#ifndef TRANSLATION_H
#define TRANSLATION_H
/*
 * @file  translation.h
 * 
 * @brief class implementing the translation process
 *
 * @author Fabio Hedayioglu
 * Contact: fheday@gmail.com
 *
 */

#include <memory>
#include <vector>
#include "circularbuffer.h"
#include "elongation_codon.h"

namespace Simulations {

class Translation {
 public:
  void loadMRNA(const std::string&);
  void loadMRNA(const std::string&, const std::string&);

  void inputMRNA(std::string); // get mRNA informed by the user
  void loadConcentrations(const std::string&);
  void loadConcentrationsFromString(const std::string&);

  void setInitiationRate(double);
  void setTerminationRate(double);

  void setIterationLimit(int);
  void setTimeLimit(double);
  void setFinishedRibosomes(int);
  void setSimulateToSteadyState(bool);
  void setSteadyStateTime(float);
  void setSteadyStateTerminations(int);

  void setPrepopulate(bool);

  void setHistorySize(std::size_t);

  void run();

  void getAverageTimes();
  std::tuple<std::vector<double>, std::vector<int>> getElongationDuration();
  void getInitiationElongationTermination();

  std::vector<int> getRibosomesPositions();
  void setRibosomePositions(std::vector<int>);

  void setLogCodonStates(bool log);
  std::vector<std::tuple<std::vector<int>, std::vector<double>>>
  getLogCodonStates();

  void getRibosomeCollisions();

  std::vector<double> initiations_durations, elongations_durations,
      terminations_durations;
  std::vector<int> initiation_iteration, termination_iteration;

  double termination_rate = -1;
  double initiation_rate = -1;
  int iteration_limit = -1;

  int finished_ribosomes_limit = -1;
  double time_limit = -1;
  bool no_noCognate = false;

  std::vector<double> alphas;  // reactions alphas - all available ones.
  std::vector<std::size_t> codon_index;  // indexes of the codon where the alpha belongs to.
  std::vector<std::size_t> reaction_index;  // in the codon, the index of the reaction.
  std::size_t global_size = 0; // written size of alphas, codon_index, reaction_index.
  std::vector<std::unique_ptr<Simulations::mRNAElement>> codons_vector;
  std::string mrna_file_name;
  std::string gene_name;
  std::string mrna_input = "";
  std::string concentrations_file_name;
  std::string concentrations_string;
  std::string concentrations_source = "None";

  std::string saccharomyces_cerevisiae_concentrations = 
    "concentrations/Saccharomyces_cerevisiae.csv";


  std::vector<double> dt_history;
  std::vector<std::vector<int>> ribosome_positions_history;

  void setPropensities(std::vector<std::map<std::string, double>> prop);
  void setNoNonCognate(bool noNonCog);
  std::vector<std::map<std::string, double>> getPropensities();

  // array with the total times the ribosomes spent in the codons
  std::vector<double> total_time;
  // number of times a codon was occupied
  std::vector<int> n_times_occupied;
  // average occupation time
  std::vector<double> codons_average_occupation_time;
  // ribosomes colliding
  std::vector<std::vector<int>> colliding_ribosomes;
  //stalled ribosomes (blocking others to elongate)
  std::vector<std::vector<int>> stalled_ribosomes;

 private:
  void initializeMRNAReader();
  void insertRibosome(std::size_t, bool);
  bool pre_populate = false;
  bool changed_propensities = false;
  bool is_logging_codon_state = false;
  bool is_initiation_set = false;
  bool is_termination_set = false;
  bool is_mRNA_valid = true;
  bool is_collisions_calculated = false;
  std::size_t history_size = 100000;
  bool simulate_to_steady_state = false;
  int steady_state_terminations = -1;
  float steady_state_time = -1.0f;
  void getAlphas(utils::circular_buffer<std::vector<int>>&);
};
}  // namespace Simulations
#endif  // TRANSLATION_H
