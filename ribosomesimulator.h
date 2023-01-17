#ifndef SIMULATIONS_RIBOSOMESIMULATOR_H
#define SIMULATIONS_RIBOSOMESIMULATOR_H

/*
 * @file  ribosomesimulator.h
 * 
 * @brief class where a codon is represented and could be individually simulated
 *
 * @author Fabio Hedayioglu
 * Contact: fheday@gmail.com
 *
 */

#include <array>
#include <functional>
#include <map>
#include <random>
#include <tuple>
#include <vector>
#include "concentrationsreader.h"

namespace Simulations {

class RibosomeSimulator {
 public:
  RibosomeSimulator();
  int getState();
  void setState(int);
  void getAlphas(std::vector<double>&, std::vector<int>&);
  void getDecodingAlphas(std::vector<double>&, std::vector<int>&);

  void setPropensities(std::map<std::string, double> prop);
  void setPropensity(std::string&, const double);
  double getPropensity(std::string);
  void setNonCognate(double noNonCog);

  std::map<std::string, double> getPropensities();
  void loadConcentrations(const std::string&);
  void loadConcentrationsFromString(const std::string&);
  void setCodonForSimulation(const std::string&);
  void run_and_get_times(double&, double&);
  double run_repeatedly_get_average_time(const int);
  std::vector<double> dt_history;
  std::vector<int> ribosome_state_history;
  std::string saccharomyces_cerevisiae_concentrations = 
       "concentrations/Saccharomyces_cerevisiae.csv";
  // propensity identifyers
  std::array<std::string, 44> reactions_identifiers = {
      {"non1f",    "near1f",     "wobble1f", "WC1f",     "non1r",    "near1r",
       "near2f",   "near2r",     "near3f",   "near4f",   "near5f",   "neardiss",
       "near6f",   "wobble1r",   "wobble2f", "wobble2r", "wobble3f", "wobble4f",
       "wobble5f", "wobblediss", "wobble6f", "WC1r",     "WC2f",     "WC2r",
       "WC3f",     "WC4f",       "WC5f",     "WCdiss",   "WC6f",     "dec7f",
       "trans1f",  "trans1r",    "trans2",   "trans3",   "trans4",   "trans5",
       "trans6",   "trans7",     "trans8",   "trans9"}};

 private:
  std::random_device rd;
  std::mt19937 gen;
  std::uniform_real_distribution<> dis;
  void buildReactionsMap();
  std::string simulation_codon_3_letters = "";
  csv_utils::ConcentrationsReader concentrations_reader;
  std::vector<std::vector<std::tuple<std::reference_wrapper<double>, int>>>
  createReactionsGraph(const csv_utils::concentration_entry&);
  std::map<
      std::string,
      std::vector<std::vector<std::tuple<std::reference_wrapper<double>, int>>>>
      reactions_map;
  std::vector<std::vector<std::tuple<std::reference_wrapper<double>, int>>>
      reactions_graph;  // vector where the index is the ribosome's current
                        // state and the content is a vector of tuples
                        // containing the propensity and next state of each
                        // possible reaction.
  int current_state = 0;

  std::vector<std::string> stop_codons = {"UAG", "UAA", "UGA"};
  // constants for WCcognate interaction in 1/sec
  std::map<std::string, double> WC1f;
  double WC1r = 85;
  double WC2f = 190;
  double WC2r = 0.23;
  double WC3f = 260;
  double WC4f = 1000;
  double WC5f = 1000;
  double WCdiss = 60;
  double WC6f = 1000;
  double dec7f = 200;

  // constants for wobblecognate interaction in 1/sec
  std::map<std::string, double> wobble1f;
  double wobble1r = 85;
  double wobble2f = 190;
  double wobble2r = 1;
  double wobble3f = 25;
  double wobble4f = 1000;
  double wobble5f = 1000;
  double wobblediss = 1.1;
  double wobble6f = 6.4;

  // constants for nearcognate interaction in 1/sec
  std::map<std::string, double> near1f;
  double near1r = 85;
  double near2f = 190;
  double near2r = 80;
  double near3f = 0.4;
  double near4f = 1000;
  double near5f = 1000;
  double neardiss = 1000;
  double near6f = 60;

  double totalconc = 1.9e-4;

  // constants for noncognate interaction in 1/sec.
  // Non-cognates are assumed to not undergo any significant
  // interaction but to simply dissociate quickly.
  std::map<std::string, double> non1f;
  double non1r = 1e5;

  // based on yeast value of 226000 molecules per cell as determined
  // in von der Haar 2008 (PMID 18925958)
  double eEF2conc = 1.36e-5;
  // constants for translocation in 1/sec
  // 150 uM-1 s-1 = is from Fluitt et al 2007 (PMID 17897886)
  double trans1f = eEF2conc * 1.5e8;
  double trans1r = 140;
  double trans2 = 250;
  double trans3 = 350;
  double trans4 = 1000;
  double trans5 = 1000;
  double trans6 = 1000;
  double trans7 = 1000;
  double trans8 = 1000;
  double trans9 = 1000;

  // TEMPORARY SOLUTION: in order to provide a mechanism easily change only one
  // propensity I've created a map of references to the propensities variables,
  // but this must be updated manually if in the future we add or remove
  // reactions.
  std::map<std::string, double*> propensities_map;
};

}  // namespace Simulations

#endif  // SIMULATIONS_RIBOSOMESIMULATOR_H
