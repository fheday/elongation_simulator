#include <getopt.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include "concentrationsreader.h"
#include <map>
#include "ribosomesimulator.h"


/**
 * @brief use the concentrations informed in concentrations_file_name, execute
 * the number of informed iterations and then calculates the average decoding
 * and translocating times, writing the output as a csv file into
 * average_times_file_name. This function is usually called from the function
 * calculate_codons_propensities and used by it.
 *
 * This procedure should usually be used only for initializing values for the
 * ElogationSimulator class.
 *
 * @param concentrations_file_name string containing the path to the csv file
 * containing the concentrations in the cell.
 * @param iterations number of iterations to run per codon base.
 * @param average_times_file_name string containing the path to write the
 * average times calculated by the algorithm.
 * @param translocating_times boolean if true, all codons have only decoding
 * times and the translocating time is represented by the codon 'tra'. If false,
 * all the codons times are decoding + translocating.
 * @param states_file_name string containing the path to write the csv file with
 * the times spent in the ribosome's states
 * @return std::map< std::__cxx11::string, double > a map with codons and
 * average decoding times. Average Translocating time is given by entry 'tra'
 * 
 * @author Fabio Hedayioglu
 * Contact: fheday@gmail.com
 */
std::map<std::string, double> calculate_codons_times(
    const std::string& concentrations_file_name, int iterations,
    const std::string& average_times_file_name,
    const std::string& times_vector_file_name,
    const std::string& sates_file_name, bool translocating_times) {
  Simulations::RibosomeSimulator ribosome;
  ribosome.loadConcentrations(concentrations_file_name);
  csv_utils::ConcentrationsReader cr;
  cr.loadConcentrations(concentrations_file_name);
  double decoding = 0, translocating = 0;
  std::vector<std::string> codons;
  std::map<std::string, double> codons_times;
  cr.getCodonsVector(codons);

  double total_translocating = 0, codon_total_decoding = 0, n = 0;
  std::ofstream average_times_file;
  std::ofstream times_vector_file;
  // set numbers precision in the files.
  average_times_file << std::setprecision(15);
  times_vector_file << std::setprecision(15);
  // open the files for writing.
  average_times_file.open(average_times_file_name);
  times_vector_file.open(times_vector_file_name);
  // create header line.
  average_times_file << "codon, time\n";
  times_vector_file << "codon";
  for (int i = 0; i < 2 * iterations; i++) {
    times_vector_file << ", V" << i;
  }
  times_vector_file << "\n";
  // calculate times and generate the vectors.
  std::vector<double> vector(static_cast<std::size_t>(2 * iterations), 0);
  double codon_total_translocating = 0;
  for (const std::string& codon : codons) {
    codon_total_decoding = 0;
    codon_total_translocating = 0;
    ribosome.setCodonForSimulation(codon);
    ribosome.setNonCognate(0.0);
    average_times_file << "\"" << codon << "\"";
    std::cout << "Starting codon: " << codon;
    if (!translocating_times) {
      codon_total_decoding = ribosome.run_repeatedly_get_average_time(iterations) * iterations;
      codon_total_translocating = 0;
    } else {
      for (int i = 0; i < iterations; i++) {
        ribosome.setState(0);
        ribosome.run_and_get_times(decoding, translocating);
      if (decoding * translocating <= std::numeric_limits<double>::epsilon()) {
        throw std::runtime_error("decoding nor translocation cannot be zero.");
      }
        codon_total_decoding += decoding;
        codon_total_translocating += translocating;
        n++;
        // save vector.
        vector[static_cast<std::size_t>(i)] = decoding;
        vector[static_cast<std::size_t>(iterations + i)] = translocating;
      }
    }
    total_translocating += codon_total_translocating;
    // write times and vector to files.
    average_times_file << ", ";
    if (translocating_times) {
      average_times_file << (codon_total_decoding) / iterations;
    } else {
      average_times_file << (codon_total_decoding + codon_total_translocating) /
                                iterations;
    }
    average_times_file << "\n";
    for (int j = 0; j < (2 * iterations) - 1; j++) {
      times_vector_file << vector[static_cast<std::size_t>(j)] << ",";
    }
    times_vector_file << vector[static_cast<std::size_t>((2 * iterations)) - 1]
                      << vector[static_cast<std::size_t>((2 * iterations)) - 1]
                      << "\n";
    codons_times[codon] = decoding;
    std::cout << ". Finished. Average time: " ;
    if (translocating_times) {
      std::cout << codon_total_decoding / iterations;
    } else {
      std::cout << (codon_total_decoding + codon_total_translocating) /
                    iterations;
    }
    std::cout << "\n";
  }
  // save translocation times.
  if (translocating_times) {
    average_times_file << "tra, " << (total_translocating / n) << "\n";
    codons_times["tra"] = (total_translocating / n);
    std::cout << "Average translocating time: " << (total_translocating / n)
              << "\n";
  }

  if (!sates_file_name.empty()) {
    // save states file.
  }
  // close files.
  average_times_file.close();
  times_vector_file.close();
  return codons_times;
}

void printHelp() {
  std::cout << "Wrong number of parameters informed.\n";
  std::cout << "Usage: calculateCodonsTimes  -c concentrations_file_name "
               "-l iterations -a average_times_file_name -v "
               "times_vector_file_name -t\n";
  std::cout << "concentrations_file_name  - The path to the csv file "
               "containing the concentrations in the cell.\n";
  std::cout << "iterations - Number of iterations to run per codon base.\n";
  std::cout << "average_times_file_name - The path to write the average "
               "times calculated by the algorithm.\n";
  std::cout << "times_vector_file_name - The path to write the vector with "
               "the calculated times by the algorithm.\n";
  std::cout << "-t - if passed, all codons have only "
               "decoding times and the translocating time is represented by "
               "the codon 'tra'. If ommited, all the codons times are "
               "decoding + translocating.\n";
}

int main(int argc, char** argv) {
  const char* const short_opts = "c:l:ta:v:s:h";
  const option long_opts[] = {
      {"concentration", 1, nullptr, 'c'}, {"iterations", 1, nullptr, 'l'},
      {"translocation", 0, nullptr, 't'}, {"averageFile", 1, nullptr, 'a'},
      {"vector", 1, nullptr, 'v'},        {"statesFile", 1, nullptr, 's'},
      {"help", 0, nullptr, 'h'},          {nullptr, 0, nullptr, 0}};

  std::string concentration_file, average_file, vector_file, states_file;
  int iterations = 0;
  bool translocation = false;

  if (argc == 1) {
    // no options given. Print help and exit program.
    printHelp();
    return 0;
  }
  while (optind < argc) {
    const auto opt = getopt_long(argc, argv, short_opts, long_opts, nullptr);
    if (opt != -1) {
      // Option argument
      switch (opt) {
        case 'c': {
          concentration_file = std::string(optarg);
          break;
        }
        case 'l':
          iterations = std::stoi(optarg);
          break;
        case 't':
          translocation = true;
          break;
        case 'a':
          average_file = std::string(optarg);
          break;
        case 'v':
          vector_file = std::string(optarg);
          break;
        case 's':
          states_file = std::string(optarg);
          break;
        case 'h':  // -h or --help
        case '?':  // Unrecognized option
        default:
          printHelp();
          return 0;
      }
    } else {
      break;
    }
  }
  // check if we have all we need:
  if (concentration_file.empty() || iterations == 0) {
    printHelp();
    return 0;
  }
  calculate_codons_times(concentration_file, iterations, average_file,
                         vector_file, states_file, translocation);

  return 0;
}
