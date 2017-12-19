#include "concentrationsreader.h"
#include <error.h>
#include <algorithm>
#include <fstream>
#include <string>

csv_utils::ConcentrationsReader::ConcentrationsReader() { contents.clear(); }

void csv_utils::ConcentrationsReader::loadConcentrations(
    const std::string& file_name) {
  std::ifstream ist{file_name};

  if (!ist) {
    throw std::runtime_error("can't open input file: " + file_name);
  }
  contents.clear();
  std::string codon;
  std::string three_letter;
  double wc_cognate_conc;
  double wobblecognate_conc;
  double nearcognate_conc;
  std::string tmp_str;
  std::vector<std::string> stop_codons = {"UAG", "UAA",
                                          "UGA"};  // list of stop codons.
  bool header = true;
  while (ist.good()) {
    if (!std::getline(ist, codon, ',')) {
      break;
    }
    std::getline(ist, three_letter, ',');
    std::getline(ist, tmp_str, ',');
    wc_cognate_conc = std::atof(tmp_str.c_str());
    std::getline(ist, tmp_str, ',');
    wobblecognate_conc = std::atof(tmp_str.c_str());
    std::getline(ist, tmp_str, '\n');
    nearcognate_conc = std::atof(tmp_str.c_str());
    codon.erase(std::remove(codon.begin(), codon.end(), '\"'),
                codon.end());  // remove \" from string.
    three_letter.erase(
        std::remove(three_letter.begin(), three_letter.end(), '\"'),
        three_letter.end());  // remove \" from string.
    if (!header) {
      auto result = std::find(stop_codons.begin(), stop_codons.end(), codon);
      // only add if not a stop codon.
      if (result == end(stop_codons)) {
        contents.push_back(
            concentration_entry{codon, three_letter, wc_cognate_conc,
                                wobblecognate_conc, nearcognate_conc});
      }
    } else {
      header = false;
    }
  }
}

void csv_utils::ConcentrationsReader::getContents(
    std::vector<concentration_entry>& result) {
  result = contents;
}

void csv_utils::ConcentrationsReader::getCodonsVector(
    std::vector<std::string>& codons_vector) {
  codons_vector.clear();
  for (concentration_entry entry : contents) {
    codons_vector.push_back(entry.codon);
  }
}
