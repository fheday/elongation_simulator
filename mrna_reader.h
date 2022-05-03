#ifndef MRNA_UTILS_MRNA_READER_H
#define MRNA_UTILS_MRNA_READER_H

/*
 * @file  mrna_reader.h
 * 
 * @brief class to read mRNA FASTA files mRNA string sequences.
 *
 * @author Fabio Hedayioglu
 * Contact: fheday@gmail.com
 *
 */

#include <string>
#include <vector>

namespace mRNA_utils {

class mRNAReader {
 public:
  mRNAReader();
  void loadmRNAFile(const std::string&); // load a file containing only one gene.
  void loadGene(const std::string&, const std::string&); // load a gene inside a file.
  std::vector<std::string> get_names_in_file(const std::string&); // get all gene names in passed file.
  void readMRNA(std::string); // load mRNA informed by the user.
  std::string getCodon(int);
  int sizeInCodons();

 private:
  std::string mRNA_sequence;
  double termination_rate;
  double initiation_rate;
  std::string mRNA_file_name;
  void post_process_sequence(std::string);
};
}  // namespace mRNA_utils

#endif  // MRNA_UTILS_MRNA_READER_H
