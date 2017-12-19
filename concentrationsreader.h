#ifndef CONCENTRATIONS_READER_H
#define CONCENTRATIONS_READER_H

#include <string>
#include <vector>

namespace csv_utils {
struct concentration_entry {
  std::string codon;
  std::string three_letter;
  double wc_cognate_conc;
  double wobblecognate_conc;
  double nearcognate_conc;
};

class ConcentrationsReader {
 public:
  ConcentrationsReader();
  void loadConcentrations(const std::string&);
  void getContents(std::vector<concentration_entry>&);
  void getCodonsVector(std::vector<std::string>&);

 private:
  std::vector<concentration_entry> contents;
};
}  // namespace csv_utils
#endif  // CONCENTRATIONS_READER_H
