#include "mrna_reader.h"
#include <fstream>

using namespace mRNA_utils;

mRNAReader::mRNAReader() {}

void mRNAReader::loadmRNAFile(std::string mRNA_file_name) {
  std::ifstream ist{mRNA_file_name};

  if (!ist) {
    throw std::runtime_error("can't open input file: " + mRNA_file_name);
  }
  mRNA_sequence.clear();
  std::string line;
  while (ist.good()) {
    std::getline(ist, line);
    // some file formats start with a '>' symbol on the first line.
    // we need to skip that line.
    if (line[0] == '>') continue;
    mRNA_sequence.append(line);
  }
  // replace all T's for U's.
  std::size_t found = mRNA_sequence.find('T');
  while (found != std::string::npos) {
    mRNA_sequence[found] = 'U';
    found = mRNA_sequence.find('T', found + 1);
  }
}

void mRNAReader::setInitiationRate(double val) {
  if (val > 0) {
    initiation_rate = val;
  } else {
    throw std::runtime_error("invalid initiation rate.");
  }
}

void mRNAReader::setTerminationRate(double val) {
  if (val > 0) {
    termination_rate = val;
  } else {
    throw std::runtime_error("invalid termination rate.");
  }
}

std::string mRNA_utils::mRNAReader::getCodon(int codon_number) {
  return mRNA_sequence.substr(static_cast<std::size_t>(codon_number * 3), 3);
}
