#ifndef CONCENTRATIONS_READER_H
#define CONCENTRATIONS_READER_H

/*
 * @file  concentrationsreader.h
 * 
 * @brief interface for reading concentration file
 *
 * @author Fabio Hedayioglu
 * Contact: fheday@gmail.com
 *
 */

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
  void loadConcentrationsFromString(const std::string&);
  void getContents(std::vector<concentration_entry>&);
  void getCodonsVector(std::vector<std::string>&);

 private:
  std::vector<concentration_entry> contents;
  void readConcentratonsStream(std::istream&);
};
}  // namespace csv_utils
#endif  // CONCENTRATIONS_READER_H
