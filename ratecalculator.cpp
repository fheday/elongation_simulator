#include "ratecalculator.h"
#include <error.h>
#include <algorithm>
#include <fstream>
#include <vector>

csv_utils::RateCalculator::RateCalculator() {}

/**
 * @brief Load a file with the average times and calculates the rates of the
 * codons.
 *
 * @param file_name string with the location of the file to be open.
 */
void csv_utils::RateCalculator::loadRates(const std::string& file_name) {
  std::ifstream ist{file_name};

  if (!ist) {
    throw std::runtime_error("can't open file: " + file_name);
  }
  codon_rates.clear();
  std::string codon;
  double decoding_time;
  std::string tmp_str;
  std::vector<std::string> stop_codons = {"UAG", "UAA",
                                          "UGA"};  // list of stop codons.
  bool header = true;
  while (ist.good()) {
    if (!std::getline(ist, codon, ',')) {
      break;
    }

    std::getline(ist, tmp_str, '\n');
    decoding_time = std::atof(tmp_str.c_str());
    codon.erase(std::remove(codon.begin(), codon.end(), '\"'),
                codon.end());  // remove \" from string.
    if (!header) {
      auto result = std::find(stop_codons.begin(), stop_codons.end(), codon);
      // only add if not a stop codon.
      if (result == end(stop_codons)) {
        codon_rates[codon] = 1 / decoding_time;
      } else {
        // it is a stop codon. It has a fixed decoding rate.
        codon_rates[codon] = 1.0;
      }
    } else {
      header = false;
    }
  }
  // check if stop codons where added or not: This will depend if the file has
  // them or not: if it does, they where added in the previous loop.
  for (const std::string& codon : stop_codons) {
    auto result = std::find(stop_codons.begin(), stop_codons.end(), codon);
    if (result != end(stop_codons)) {
      // This stop codon was not previously added. Add now. It has a fixed
      // decoding rate.
      codon_rates[codon] = 1.0;
    }
  }
}
