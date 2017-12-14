#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>
#include "concentrationsreader.h"
#include "enlogationsimulator.h"
#include "gillespie.h"
#include "mrna_reader.h"
#include "ratecalculator.h"
#include "reactionsset.h"
#include "ribosomesimulator.h"

#include "translation.h"

void testGillespie() {
  // create one reaction matrix
  Eigen::MatrixXi reactions(4, 10);
  reactions.fill(0);
  reactions(1, 1) = -1;
  reactions(2, 1) = 1;
  // create a ReactionSet object and add the reaction there.
  Simulations::ReactionsSet reactions_set;
  // the reaction propensity is the second argument.
  reactions_set.addReaction(reactions, .1);
  // create a matrix with the initial species population.
  Eigen::MatrixXi population(4, 10);
  population.fill(0);
  population(1, 1) = 800;
  std::cout << "initial population = " << population << "\n------\n";
  // create the Gillespie simulator object with the population and reactions
  Simulations::Gillespie simulation(300, population, reactions_set);
  // run the simulation.
  simulation.run();
}

std::map<std::string, double> load_average_times_file(
    std::string average_times_file_name) {
  std::ifstream ist{average_times_file_name};
  std::map<std::string, double> codons_times;

  if (!ist) {
    throw std::runtime_error("can't open input file: " +
                             average_times_file_name);
  }
  bool header = true;
  std::string codon, tmp_str;
  while (ist.good()) {
    if (!std::getline(ist, codon, ',')) {
      break;
    }
    std::getline(ist, codon, ',');
    std::getline(ist, tmp_str, ',');
    if (!header) {
      codons_times[codon] = std::atof(tmp_str.c_str());
    } else {
      header = false;
    }
  }
  return codons_times;
}

// void test_mRNA_reader()
// {
//     mRNA_utils::mRNAReader mrr;
//     mrr.loadRateCalculatorFile("../data/codons/average_time.csv");
//     mrr.loadmRNAFile("../../PolisomeSimulator/mRNA/S288C_ _FSF1_coding.fsa");
//     std::cout<<"mRNA = "<<"\n"<<mrr.mRNA_sequence;
//     mrr.generateInitialPopulation();
//     std::cout<<"initial population = \n";
//     std::cout<<mrr.initial_population<<"\n";
//     mrr.generateReactions();
// }

void test_enlogation_simulator(
    std::string average_times_file_name = "../data/codons/average_time.csv",
    std::string concentrations_file_name =
        "../../RSim/data_with_times/concentrations.csv",
    double init_rate = 0.0001, int term_rate = 10, int iterations = 1000000) {
  Simulations::Translation t;
  t.loadConcentrations(concentrations_file_name);
  //      t.loadMRNA("../../PolisomeSimulator/mRNA/mRNA_sample_MinCFLuc.txt");
  t.loadMRNA("../../PolisomeSimulator/mRNA/test_small_1.txt");
  t.setInitiationRate(init_rate);
  t.setTerminationRate(term_rate);
  t.setIterationLimit(iterations);
  // run enlogation simulator.
  t.run();
  int i = 0;
  for (auto rib_pos_vec : t.ribosome_positions_history) {
    std::cout << "iteration = " << i << "positions = ";
    for (int pos : rib_pos_vec) std::cout << " " << pos;
    std::cout << "\n";
    i++;
  }
  //      es.calculateAverageTimes();
  //      i = 0;
  //      for (double time:es.codons_average_occupation_time) {
  //          std::cout<<"codon = "<< i <<", average time spent = "<<time<<"\n";
  //          i++;
  //      }
}

int main(int argc, char **argv) {
  std::cout << "Hello, world!" << std::endl;
  test_enlogation_simulator();
  //      std::map<std::string, double> decoding_times_map;
  //      calculate_codons_times(argv[1], std::stoi(argv[2]), argv[3], argv[4]);
  //       calculate_codons_times("../../RSim/vito.tRNASeq.B4/concentrations.csv",
  //       10000, "../data/codons/average_time_vito.tRNASeq.B4.csv",
  //       "../data/codons/sample_vectors_vito.tRNASeq.B4.csv");

  return 0;
}
