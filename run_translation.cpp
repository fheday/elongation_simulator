/*
 * @file  run_translation.cpp
 * 
 * @brief command line mRNA simulator
 *
 * @author Fabio Hedayioglu
 * Contact: fheday@gmail.com
 *
 */

#include "concentrationsreader.h"
#include "translation.h"
#include <execinfo.h>
#include <fstream>
#include <getopt.h>
#include <iomanip>
#include <iostream>
#include <signal.h>
#include <tuple>
#include <unistd.h>
#include <vector>

#define EIGEN_NO_DEBUG // disables Eigen's assertions.

[[noreturn]] void handler(int sig) {
  void *array[10];
  int size;

  // get void*'s for all entries on the stack
  size = backtrace(array, 10);

  // print out all the frames to stderr
  fprintf(stderr, "Error: signal %d:\n", sig);
  backtrace_symbols_fd(array, size, STDERR_FILENO);
  exit(1);
}

void execute_translation(const std::string &concentrations_file,
                         const std::string &mrna_file, 
                         const std::string &gene_name, double initiation_rate,
                         double termination_rate, int time_limit,
                         int number_iterations, int number_ribosomes,
                         bool pre_fill_mRNA, bool save_ribosomes_positions,
                         const std::string &output_file_name) {
  // separate the path from the file name.
  std::size_t found = output_file_name.find_last_of("/\\");
  std::string path = "./"; // current path.
  std::string file_name = output_file_name;
  if (found != std::string::npos) {
    // there is a path.
    path = output_file_name.substr(0, found + 1);
    file_name = output_file_name.substr(found + 1);
  }

  // prepare and run the simulation.
  Simulations::Translation ts;
  ts.loadConcentrations(concentrations_file);
  if (gene_name.empty()){
    ts.loadMRNA(mrna_file);
  } else {
    ts.loadMRNA(mrna_file, gene_name);
  }
  ts.setInitiationRate(initiation_rate);
  ts.setTerminationRate(termination_rate);
  if (time_limit > 0) {
    ts.setTimeLimit(time_limit);
  } else if (number_iterations > 0) {
    ts.setIterationLimit(number_iterations);
  } else if (number_ribosomes > 0) {
    ts.setFinishedRibosomes(number_ribosomes);
  }
  ts.setPrepopulate(pre_fill_mRNA); // simulations pre-populate the mRNA by
                                    // default. This can be changed in the
                                    // future.
  ts.run();
  if (save_ribosomes_positions == false) {
    ts.getAverageTimes();
    std::vector<double> elongation_duration;
    std::vector<int> iteration_initiation;
    std::tie(elongation_duration, iteration_initiation) =
        ts.getElongationDuration();
    // save elongation data into csv file.

    std::vector<double> clock(ts.dt_history.size()), clock_at_initiation(iteration_initiation.size());

    std::partial_sum(ts.dt_history.begin(), ts.dt_history.end(), clock.begin(), std::plus<double>());

    // get the clock at the initiation of each terminating ribosome.
    std::transform(iteration_initiation.begin(), iteration_initiation.end(), clock_at_initiation.begin(), [&] (auto iteration){return clock[static_cast<std::size_t>(iteration)];});

    // now we save the clock_at_initiation and elongation_duration.
    std::ofstream clock_and_elongation_csv_file;
    clock_and_elongation_csv_file.open(output_file_name);
    // header
    clock_and_elongation_csv_file
        << "Clock at initiation, Ribosome elongation duration\n";
    // data
    for (std::size_t i = 0, total = clock_at_initiation.size(); i < total;
    ++i) {
      clock_and_elongation_csv_file << std::fixed << std::setprecision(10)
                                    << clock_at_initiation[i] << ", "
                                    << elongation_duration[i] << "\n";
    }
    clock_and_elongation_csv_file.close();

    std::ofstream codon_average_time_file;
    codon_average_time_file.open(path + "codon_average_time_" + file_name);
    // header
    codon_average_time_file << "codon average time\n";
    // data
    for (auto average_occupation_time : ts.codons_average_occupation_time) {
      codon_average_time_file << std::fixed << std::setprecision(10)
                              << average_occupation_time << "\n";
    }
    codon_average_time_file.close();
  } else {
    // save file with dt,  ribosome positions.
    std::ofstream ribosome_positions_csv_file;
    ribosome_positions_csv_file.open(output_file_name);
    // header
    ribosome_positions_csv_file << "dt,ribosomes\n";
    // data
    for (std::size_t i = 0, total = ts.dt_history.size(); i < total; ++i) {
      ribosome_positions_csv_file << std::fixed << std::setprecision(10)
                                    << ts.dt_history[i] << ",\"[";
                                    bool first = true;
                                    for (auto elem:ts.ribosome_positions_history[i]) {
                                      if (!first) {
                                        ribosome_positions_csv_file<<", ";
                                      };
                                      first = false;
                                      ribosome_positions_csv_file << elem;
                                    }
                                    ribosome_positions_csv_file << "]\"\n";
    }
    ribosome_positions_csv_file.close();

  }
}

void printHelp() {
  std::cout << "Wrong number of parameters informed.\n";
  std::cout << "run_translation concentrations mRNA Initiation Termination "
               "Time output\n";
  std::cout << "Concentrations = path to the file containing the "
               "concentrations to be used in the simulation.\n";
  std::cout << "mRNA = path to the file with the mRNA to be used.\n";
  std::cout << "Initiation = value to be used as the initiation factor.\n";
  std::cout << "Termination = value to be used as the termination factor.\n";
  std::cout << "Time limit = time limit for when the simulation should stop. "
               "This is in yeast time, not in real life time\n";
  std::cout << "output = file to be created with the simulation results.\n";
}

int main(int argc, char **argv) {
  signal(SIGSEGV, handler); // install our handler
  const char *const short_opts = "c:m:g:i:t:y:r:l:peo:h";
  const option long_opts[] = {
      {"concentration", 1, nullptr, 'c'}, {"mrna", 1, nullptr, 'm'}, {"gene", 1, nullptr, 'g'},
      {"initiation", 1, nullptr, 'i'},    {"termination", 1, nullptr, 't'},
      {"yeasttime", 1, nullptr, 'y'},     {"ribosomes", 1, nullptr, 'r'},
      {"iterations", 1, nullptr, 'l'},    {"positions", 0, nullptr, 'p'},
      {"init_empty", 0, nullptr, 'e'},    {"output", 1, nullptr, 'o'},
      {"help", 0, nullptr, 'h'},          {nullptr, 0, nullptr, 0}};

  std::string concentration_file, mrna_file, gene_name, output_file;
  double initiation = 0.0, termination = 0.0;
  int yeast_time, ribosomes, iterations;
  bool stop_condition_passed = false, pre_fill_mRNA = true,
       ribosomes_positions = false;
  yeast_time = ribosomes = iterations = -1;

  std::string halting_condition_error =
      "only one of the following halting options can be used: yeast time, "
      "terminating ribosomes, or iteration limit\n";
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
      case 'c': 
        concentration_file = std::string(optarg);
        break;
      case 'm': 
        mrna_file = std::string(optarg);
        break;
      case 'g': 
        gene_name = std::string(optarg);
        break;
      case 'i':
        initiation = std::stod(optarg);
        break;
      case 't':
        termination = std::stod(optarg);
        break;
      case 'y':
        if (!stop_condition_passed) {
          yeast_time = std::stoi(optarg);
          stop_condition_passed = true;
        } else {
          std::cout << halting_condition_error;
          return -1;
        }
        break;
      case 'r':
        if (!stop_condition_passed) {
          ribosomes = std::stoi(optarg);
          stop_condition_passed = true;
        } else {
          std::cout << halting_condition_error;
          return -1;
        }
        break;
      case 'l':
        if (!stop_condition_passed) {
          iterations = std::stoi(optarg);
          stop_condition_passed = true;
        } else {
          std::cout << halting_condition_error;
          return -1;
        }
        break;
      case 'e':
        pre_fill_mRNA = false;
        break;
      case 'p':
        ribosomes_positions = true;
        break;
      case 'o':
        output_file = std::string(optarg);
        break;
      case 'h': // -h or --help
      case '?': // Unrecognized option
      default:
        printHelp();
        return 0;
      }
    } else {
      break;
    }
  }
  // check if we have all we need:
  if (concentration_file.empty() || mrna_file.empty() || initiation == 0.0 ||
      termination == 0.0 || yeast_time == 0.0 || ribosomes == 0.0 ||
      iterations == 0.0) {
    printHelp();
    return 0;
  }
  execute_translation(concentration_file, mrna_file, gene_name, initiation, termination,
                      yeast_time, iterations, ribosomes, pre_fill_mRNA,
                      ribosomes_positions, output_file);
  return 0;
}
