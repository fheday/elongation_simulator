#ifndef MRNA_UTILS_MRNA_READER_H
#define MRNA_UTILS_MRNA_READER_H

#include <eigen3/Eigen/Dense>
#include <string>
#include "ratecalculator.h"
#include "reactionsset.h"

namespace mRNA_utils {

class mRNAReader {
 public:
  mRNAReader();
  void loadmRNAFile(std::string);
  void generateInitialPopulation();
  void generateReactions();
  void setInitiationRate(double);
  void setTerminationRate(double);
  std::string getCodon(int);
  std::string mRNA_sequence;
  double termination_rate;
  double initiation_rate;
  std::string mRNA_file_name;
  Eigen::MatrixXi initial_population;
  Simulations::ReactionsSet reactions_set;
};
}  // namespace mRNA_utils

#endif  // MRNA_UTILS_MRNA_READER_H
