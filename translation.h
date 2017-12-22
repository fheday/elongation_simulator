#ifndef TRANSLATION_H
#define TRANSLATION_H

#include <memory>
#include <vector>
#include "enlongation_codon.h"

namespace Simulations {

class Translation {
 public:
  void loadMRNA(const std::string&);
  void loadConcentrations(const std::string&);

  void setInitiationRate(double);
  void setTerminationRate(double);

  void setIterationLimit(int);
  void setTimeLimit(double);
  void setFinishedRibosomes(int);

  void setPrepopulate(bool);

  void getAlphas();
  void run();

  void calculateAverageTimes();
  std::tuple<std::vector<double>, std::vector<int>> getEnlogationDuration();
  void getInitiationEnlongationTermination();

  std::vector<double> initiations_durations, enlongations_durations,
      terminations_durations;
  std::vector<int> initiation_iteration, termination_iteration;

  double termination_rate = -1;
  double initiation_rate = -1;
  int iteration_limit = -1;
  double time_limit = -1;
  int finished_ribosomes_limit = -1;

  std::vector<double> alphas =
      std::vector<double>(10);  // reactions alphas - all available ones.
  std::vector<int> codon_index =
      std::vector<int>(10);  // indexes of the codon where the alpha belongs to.
  std::vector<int> reaction_index =
      std::vector<int>(10);  // in the codon, the index of the reaction.
  std::vector<std::unique_ptr<Simulations::mRNAElement>> codons_vector;
  std::string mrna_file_name;
  std::string concentrations_file_name;

  std::vector<double> dt_history;
  std::vector<std::vector<int>> ribosome_positions_history;

  // array with the total times the ribosomes spent in the codons
  std::vector<double> total_time;
  // number of times a codon was occupied
  std::vector<int> n_times_occupied;
  // average occupation time
  std::vector<double> codons_average_occupation_time;

 private:
  void initializeMRNAReader();
  bool pre_populate = false;
};
}  // namespace Simulations
#endif  // TRANSLATION_H
