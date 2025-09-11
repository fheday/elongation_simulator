#ifndef SIMULATIONS_CODONSIMULATOR_H
#define SIMULATIONS_CODONSIMULATOR_H

/*
 * @file  codonsimulator.h
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

class CodonSimulator {
 public:
  CodonSimulator();
  [[nodiscard]] int getState() const;
  void setState(int);
  void getAlphas(std::vector<float>&, std::vector<int>&);
  void getDecodingAlphas(std::vector<float>&, std::vector<int>&);

  void setPropensities(std::map<std::string, float>& prop);
  void setPropensity(std::string&, const float&);
  float getPropensity(std::string&);
  void setNonCognate(float noNonCog);

  std::map<std::string, float> getPropensities();
  void loadConcentrations(const std::string&);
  void loadConcentrationsFromString(const std::string&);
  void setCodonForSimulation(const std::string&);
  void run_and_get_times(float&, float&);
  float run_repeatedly_get_average_time(const int&);
  std::vector<float> dt_history;
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
  std::string simulation_codon_3_letters;
  csv_utils::ConcentrationsReader concentrations_reader;
  std::vector<std::vector<std::tuple<std::reference_wrapper<float>, int>>>
  createReactionsGraph(const csv_utils::concentration_entry&);
  std::map<
      std::string,
      std::vector<std::vector<std::tuple<std::reference_wrapper<float>, int>>>>
      reactions_map;
  std::vector<std::vector<std::tuple<std::reference_wrapper<float>, int>>>
      reactions_graph;  // vector where the index is the ribosome's current
                        // state and the content is a vector of tuples
                        // containing the propensity and next state of each
                        // possible reaction.
  int current_state = 0;

  std::vector<std::string> stop_codons = {"UAG", "UAA", "UGA"};
  // constants for WCcognate interaction in 1/sec
  std::map<std::string, float> WC1f;
  float WC1r = 85.f;
  float WC2f = 180.f;
  float WC2r = 0.22f;
  float WC3f = 260.f;
  float WC4f = 1000.f;
  float WC5f = 1000.f;
  float WCdiss = 60.f;
  float WC6f = 1000.f;
  float dec7f = 200.f;

  // constants for wobblecognate interaction in 1/sec
  std::map<std::string, float> wobble1f;
  float wobble1r = 85.f;
  float wobble2f = 190.f;
  float wobble2r = 1.f;
  float wobble3f = 25.f;
  float wobble4f = 1000.f;
  float wobble5f = 1000.f;
  float wobblediss = 1.1f;
  float wobble6f = 6.4f;

  // constants for nearcognate interaction in 1/sec
  std::map<std::string, float> near1f;
  float near1r = 85.f;
  float near2f = 190.f;
  float near2r = 80.f;
  float near3f = 0.4f;
  float near4f = 1000.f;
  float near5f = 1000.f;
  float neardiss = 1000.f;
  float near6f = 60.f;

  float totalconc = 1.9e-4f;

  // constants for noncognate interaction in 1/sec.
  // Non-cognates are assumed to not undergo any significant
  // interaction but to simply dissociate quickly.
  std::map<std::string, float> non1f;
  float non1r = 2e3;

  // based on yeast value of 226000 molecules per cell as determined
  // in von der Haar 2008 (PMID 18925958)
  float eEF2conc = 1.36e-5;
  // constants for translocation in 1/sec
  // 150 uM-1 s-1 = is from Fluitt et al 2007 (PMID 17897886)
  float trans1f = eEF2conc * 1.5e8f;
  float trans1r = 140;
  float trans2 = 250;
  float trans3 = 350;
  float trans4 = 1000;
  float trans5 = 1000;
  float trans6 = 1000;
  float trans7 = 1000;
  float trans8 = 1000;
  float trans9 = 1000;

  // TEMPORARY SOLUTION: in order to provide a mechanism easily change only one
  // propensity I've created a map of references to the propensities variables,
  // but this must be updated manually if in the future we add or remove
  // reactions.
  std::map<std::string, float*> propensities_map;
};

}  // namespace Simulations

#endif  // SIMULATIONS_CODONSIMULATOR_H
