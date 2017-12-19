#ifndef MRNA_UTILS_MRNA_READER_H
#define MRNA_UTILS_MRNA_READER_H

#include <string>

namespace mRNA_utils {

class mRNAReader {
 public:
  mRNAReader();
  void loadmRNAFile(const std::string&);
  void setInitiationRate(double);
  void setTerminationRate(double);
  std::string getCodon(int);
  int sizeInCodons();

 private:
  std::string mRNA_sequence;
  double termination_rate;
  double initiation_rate;
  std::string mRNA_file_name;
};
}  // namespace mRNA_utils

#endif  // MRNA_UTILS_MRNA_READER_H
