#ifndef MRNA_UTILS_MRNA_READER_H
#define MRNA_UTILS_MRNA_READER_H

#include <string>
#include <vector>

namespace mRNA_utils {

class mRNAReader {
 public:
  mRNAReader();
  void loadmRNAFile(const std::string&); // load a file containing only one gene.
  void loadGene(const std::string&, const std::string&); // load a gene inside a file.
  std::vector<std::string> get_names_in_file(const std::string&); // get all gene names in passed file.
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
