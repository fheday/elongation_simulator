#include "concentrationsreader.h"
#ifndef _MSC_VER
#ifndef __APPLE__
#include <error.h>
#endif
#endif

/*
 * @file  concentrationsreader.cpp
 * 
 * @brief Methods to read concentration file
 *
 * @author Fabio Hedayioglu
 * Contact: fheday@gmail.com
 *
 */

#include <algorithm>
#include <fstream>
#include <string>
#include <sstream>


csv_utils::ConcentrationsReader::ConcentrationsReader() { contents.clear(); }

void csv_utils::ConcentrationsReader::loadConcentrations(
    const std::string& file_name) {
  std::filebuf fb;
	if (!fb.open (file_name,std::ios::in)) {
		throw std::runtime_error("can't open input file: " + file_name);
	}
	std::istream ist(&fb);
	readConcentratonsStream(ist);
	fb.close();
}

void csv_utils::ConcentrationsReader::loadConcentrationsFromString(const std::string& data) {
  std::istringstream ifs (data);
	readConcentratonsStream(ifs);
}

void csv_utils::ConcentrationsReader::readConcentratonsStream(std::istream& inputStream){
  contents.clear();
  std::string codon;
  std::string three_letter;
  double wc_cognate_conc;
  double wobblecognate_conc;
  double nearcognate_conc;
  int wc_cognate_conc_index = -1;
  int wobblecognate_conc_index = -1;
  int nearcognate_conc_index = -1;
  int codon_index = -1;
  int three_letter_index = 0;
  std::string tmp_str;
  std::vector<std::string> stop_codons = {"UAG", "UAA",
                                          "UGA"};  // list of stop codons.
  bool header = true;
  while (inputStream.good()) {
    if (header){
      // get the indexes of the columns of interest.
      std::getline(inputStream,tmp_str);
      // make header all lowercase.
      transform(tmp_str.begin(), tmp_str.end(), tmp_str.begin(), ::tolower);
      //remove \r, \n and non-printable characters from the header.
      tmp_str.erase(std::remove_if(tmp_str.begin(), tmp_str.end(), 
        [](unsigned char x){return std::isspace(x);}), tmp_str.end());
      // remove double quotes.
      tmp_str.erase(std::remove(tmp_str.begin(), tmp_str.end(), '\"'), tmp_str.end());
      std::stringstream ss(tmp_str);
      int index = 0;
      std::string column_name;
      while (std::getline(ss, column_name, ',')) {
        if (column_name == "codon"){
          codon_index = index;
        } else if (column_name == "three.letter") {
          three_letter_index = index;
        } else if (column_name == "wccognate.conc") {
          wc_cognate_conc_index = index;
        } else if (column_name == "wobblecognate.conc"){
          wobblecognate_conc_index = index;
        } else if (column_name == "nearcognate.conc") {
          nearcognate_conc_index = index;
        }
        index++;
      }
      if (codon_index < 0) throw std::runtime_error("no codon column in csv file.");
      if (three_letter_index < 0) throw std::runtime_error("no three.letter column in csv file.");
      if (wc_cognate_conc_index < 0) throw std::runtime_error("no WCcognate.conc column in csv file.");
      if (wobblecognate_conc_index < 0) throw std::runtime_error("no wobblecognate.conc column in csv file.");
      if (nearcognate_conc_index < 0) throw std::runtime_error("no nearcognate.conc column in csv file.");

      header = false;
    } else {
      std::getline(inputStream,tmp_str);
      std::stringstream ss(tmp_str);
      int curr_index = 0;
      while (std::getline(ss, tmp_str, ',')) {
        if (curr_index == codon_index) {
          tmp_str.erase(std::remove(tmp_str.begin(), tmp_str.end(), '\"'), tmp_str.end());
          codon = tmp_str;
        } else if (curr_index == three_letter_index) {
          three_letter = tmp_str;
        } else if (curr_index == wc_cognate_conc_index){
            wc_cognate_conc = std::atof(tmp_str.c_str());
        } else if (curr_index == wobblecognate_conc_index) {
            wobblecognate_conc = std::atof(tmp_str.c_str());
        } else if (curr_index == nearcognate_conc_index) {
            nearcognate_conc = std::atof(tmp_str.c_str());
        }
        curr_index++;
      }
      auto result = std::find(stop_codons.begin(), stop_codons.end(), codon);
      // only add if not a stop codon.
      if (result == end(stop_codons)) {
        contents.push_back(
            concentration_entry{codon, three_letter, wc_cognate_conc,
                                wobblecognate_conc, nearcognate_conc});
      }
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
