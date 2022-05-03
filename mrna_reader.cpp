/*
 * @file  mrna_reader.cpp
 * 
 * @brief class to read mRNA FASTA files mRNA string sequences.
 *
 * @author Fabio Hedayioglu
 * Contact: fheday@gmail.com
 *
 */

#include "mrna_reader.h"
#include <fstream>
#include <regex>

mRNA_utils::mRNAReader::mRNAReader() {
  termination_rate = 0;
  initiation_rate = 0;
}

void mRNA_utils::mRNAReader::post_process_sequence(std::string raw_sequence) {
  // replace all T's for U's.
    std::size_t found = raw_sequence.find('T');
    while (found != std::string::npos) {
      raw_sequence[found] = 'U';
      found = raw_sequence.find('T', found + 1);
    }
    // remove all non-printable characgters, if any.
    raw_sequence.erase(std::remove_if(raw_sequence.begin(), raw_sequence.end(), 
        [](unsigned char x){return std::isspace(x);}), raw_sequence.end());
    mRNA_sequence = raw_sequence;
}

void mRNA_utils::mRNAReader::loadmRNAFile(const std::string& mRNA_file_name) {
  std::ifstream ist{mRNA_file_name};
  if (!ist) {
    throw std::runtime_error("can't open input file: " + mRNA_file_name);
  }
  std::string raw_sequence;
  std::string line;
  bool found_first_gene = false;
  while (ist.good()) {
    std::getline(ist, line);
    // some file formats start with a '>' symbol on the first line.
    // we need to skip that line.
    if (line[0] == '>') {
      if (!found_first_gene){
        found_first_gene = true;
        continue;
      } else {
        //we are entering a new gene. No need.
        break;
      }
    }
    raw_sequence.append(line);
  }
  ist.close(); // close file.
  post_process_sequence(raw_sequence);
}

void mRNA_utils::mRNAReader::readMRNA(std::string user_mRNA) {
  post_process_sequence(user_mRNA);
}

void mRNA_utils::mRNAReader::loadGene(const std::string& mRNA_file_name, const std::string& gene_name) {
  std::ifstream ist{mRNA_file_name};
  if (!ist) {
    throw std::runtime_error("can't open input file: " + mRNA_file_name);
  }
  std::string raw_sequence;
  std::string line;
  bool found_gene = false;
  while (ist.good()) {
    std::getline(ist, line);
    // some file formats start with a '>' symbol on the first line.
    // we need to skip that line.
    if (!found_gene && line[0] == '>' && line.find(gene_name) !=  std::string::npos)  {
      found_gene = true;
      continue;
    }
    if (found_gene){
      if (line[0] != '>'){
        raw_sequence.append(line);
      } else {
        // we reached a new gene. no need to progress further.
        break;
      }
      
    }
  }
  ist.close(); // close file.
  post_process_sequence(raw_sequence);
}

std::vector<std::string> mRNA_utils::mRNAReader::get_names_in_file(const std::string& mRNA_file_name) {
  std::ifstream ist{mRNA_file_name};
  if (!ist) {
    throw std::runtime_error("can't open input file: " + mRNA_file_name);
  }
  std::string line;
  std::vector<std::string> result;
  std::regex word_regex("(\\w+)");
  
  while (ist.good()) {
    std::getline(ist, line);
    if (line[0] == '>')  {
      auto words_first = std::sregex_iterator(line.begin(), line.end(), word_regex);
      result.push_back(words_first->str());
      continue;
    }
  }
  ist.close(); // close file.
  return result;
}

std::string mRNA_utils::mRNAReader::getCodon(int codon_number) {
  return mRNA_sequence.substr(static_cast<std::size_t>(codon_number * 3), 3);
}

int mRNA_utils::mRNAReader::sizeInCodons() {
  return static_cast<int>(mRNA_sequence.size() / 3);
}
