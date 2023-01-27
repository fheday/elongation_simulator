/*
 * @file  ribosomesimulator.cpp
 * 
 * @brief class where a codon is represented and could be individually simulated
 *
 * @author Fabio Hedayioglu
 * Contact: fheday@gmail.com
 *
 */

#include "ribosomesimulator.h"
#include <float.h>
#include <algorithm>
#ifndef _MSC_VER
#include <eigen3/Eigen/Dense>
#else
#include <Eigen/Dense>
#endif
#include <numeric>
#include "concentrationsreader.h"

#if defined(COMIPLE_PYTHON_MODULE) || defined(RIBOSOMESIMULATOR)

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

PYBIND11_MODULE(ribosomesimulator, mod)
{
  mod.doc() = R"pbdoc(
    Module for simulating ribosomes outside a mRNA contruct. 
    This module simulates ribosomes decoding a given codon using a given concentration.
    )pbdoc";
  py::class_<Simulations::RibosomeSimulator>(mod, "ribosomesimulator")
      .def(py::init<>(),"Creates an empty simulator") // constructor
      .def("loadConcentrations",
           &Simulations::RibosomeSimulator::loadConcentrations, py::arg("file_name"),R"docstr(
             Loads a csv file containing the concentrations to be used in this simulation.
             
             file_name: string with the path to the file containing the concentrations.
             )docstr")
      .def("loadConcentrationsFromString",
           &Simulations::RibosomeSimulator::loadConcentrationsFromString, py::arg("string"),R"docstr(
             Loads a string similar to the csv file containing the concentrations to be used in this simulation.
             string: string with the csv values of the concentrations.
             )docstr")
      .def("setCodonForSimulation",
           &Simulations::RibosomeSimulator::setCodonForSimulation, R"docstr(
             Select the codon to be simulated. A simulator can simulate the decoding of only one codon.
             codon: 3-letter string with codon to be simulated. 
           )docstr")
      .def("setState", &Simulations::RibosomeSimulator::setState,py::arg("target_state"), R"docstr(
        Optional method: Set the ribosome's state accordingly to the reactions map.
        When creating a simaultion, the state is zero.
        target_state: State to set the ribosome.
      )docstr")
      .def("run_and_get_times",
           [](Simulations::RibosomeSimulator &rs) {
             double d = 0.0;
             double t = 0.0;
             rs.run_and_get_times(d, t);
             return std::make_tuple(d, t);
           }, R"docstr(
             Run the simulation and returns a tuple where:
             First term: total decoding time in seconds.
             Second term: total translocation time in seconds.
           )docstr")
      .def("run_repeatedly_get_average_time",&Simulations::RibosomeSimulator::run_repeatedly_get_average_time, R"docstr(
          Runs a simulation a given number of times and return the average translation time.
          repetitions: the number of times the simulation is being run.
          return: the average translation time of the simulations.
      )docstr")
      .def("setPropensities", &Simulations::RibosomeSimulator::setPropensities, py::arg("prop"), R"docstr(
        This method changes the reactions propensities of the codon selected for simulation.
        prop: dictionary with new propensities. An initial dictionary can be acquired by calling getPropensities().
      )docstr")
      .def("setNonCognate", &Simulations::RibosomeSimulator::setNonCognate, py::arg("nonCognatePropensity"), R"docstr(
        Set the propensity of non-cognates for the selected codon.
        To use this function correctly, we must have set the codon for simulation.
        nonCognatePropensity: the propensity of non-coganates in reactions/sec 
      )docstr")
      .def("getPropensities", &Simulations::RibosomeSimulator::getPropensities, R"docstr(
        This method returns a dictionary with the reactions labels and their propensities.
        This method should be used after the setCodonForSimulation.
        The dictionary returned by this method can be changed and used as an input parameter for 
        setPropensities, in order to change a specific reaction's propensity. 
      )docstr")
      .def("getPropensity", &Simulations::RibosomeSimulator::getPropensity, py::arg("reaction"), R"docstr(
        This method returns the propensity of the given reaction label.
        reaction: string with the propensity label.
        return: reaction's propensity in reactions/sec.
      )docstr")
      .def("setPropensity", &Simulations::RibosomeSimulator::setPropensity, R"docstr(
        given a reaction, sets its propensity.
        reaction: string with the reaction label
        propensity: float with new propensity value.
      )docstr")
      .def_readonly("dt_history", &Simulations::RibosomeSimulator::dt_history, R"docstr(
        Attribute with the time taken by each reaction. This numpy array is filled after a simulation has been run.
      )docstr")
      .def_readonly("ribosome_state_history",
                    &Simulations::RibosomeSimulator::ribosome_state_history, R"docstr(
                      Attribute with the current state of the ribosome at each point in the simulation.
                      Each entry in this numpy array corresponds to an entry on dt_history at the same line.
                    )docstr")
      .def_property_readonly("saccharomyces_cerevisiae_concentrations",
                             [](py::object) {
                               py::object conc_path = py::module::import("concentrations"); // load module
                               std::string file_name = "/Saccharomyces_cerevisiae.csv";     // file name
                               std::string conc_path_string;
                               for (auto item : conc_path.attr("__path__"))
                               { // iterate the path list
                                 //cast to string and concatenate with file to form proper path.
                                 conc_path_string = std::string(item.cast<py::str>());
                                 break;
                               }
                               return conc_path_string + file_name;
                             }, R"docstr(
                               This attribute can be use as a parameter when setting the concentrations file to the saccharomyces cerevisiae.
                               E.g: sim.loadConcentrations(sim.saccharomyces_cerevisiae_concentrations)
                             )docstr");
}
#endif

Simulations::RibosomeSimulator::RibosomeSimulator() : gen(rd()), dis(0, 1)
{
  // set initial state to 0
  current_state = 0;
  // initialize the random generator
  // std::random_device
  //     rd;                 // Will be used to obtain a seed for the random number engine
  // std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
  // std::uniform_real_distribution<> dis(0, 1);

}

/**
 * Loads the concentrations file.
 *
 * This method loads the csv file containing the tRNA concentrations
 * to be simulated.
 *
 * @param file_name string with the path to the file containing the concentrations.
 * @return void
 */
void Simulations::RibosomeSimulator::loadConcentrations(
    const std::string &file_name)
{
  concentrations_reader.loadConcentrations(file_name);
  buildReactionsMap();
}

void Simulations::RibosomeSimulator::loadConcentrationsFromString(const std::string& data) {
  concentrations_reader.loadConcentrationsFromString(data);
  buildReactionsMap();
}

void Simulations::RibosomeSimulator::buildReactionsMap()
{
  std::vector<csv_utils::concentration_entry> codons_concentrations;
  concentrations_reader.getContents(codons_concentrations);
  reactions_map.clear(); // make sure the map is clear.
  for (csv_utils::concentration_entry entry : codons_concentrations)
  {
    auto result =
        std::find(stop_codons.begin(), stop_codons.end(), entry.codon);
    if (result == end(stop_codons))
    {
      // Not a stop codon. Proceed.
      double nonconc = totalconc - entry.wc_cognate_conc -
                       entry.wobblecognate_conc - entry.nearcognate_conc;
      // constants for WCcognate interaction in 1/sec
      WC1f[entry.codon] = 1.4e8 * entry.wc_cognate_conc;

      // constants for wobblecognate interaction in 1/sec
      wobble1f[entry.codon] = 1.4e8 * entry.wobblecognate_conc;

      // constants for nearcognate interaction in 1/sec
      near1f[entry.codon] = 1.4e8 * entry.nearcognate_conc;

      // constants for noncognate interaction in 1/sec.
      // Non-cognates are assumed to not undergo any significant
      // interaction but to simply dissociate quickly.
      non1f[entry.codon] = 1.4e8 * nonconc;

      reactions_map[entry.codon] = createReactionsGraph(entry);
    }
  }
  if (!simulation_codon_3_letters.empty())
  {
    reactions_graph = reactions_map.at(simulation_codon_3_letters);
  }
}

void Simulations::RibosomeSimulator::setPropensities(
    std::map<std::string, double> prop)
{
  for (auto it : prop)
  {
    if (std::find(reactions_identifiers.begin(), reactions_identifiers.end(), it.first) != reactions_identifiers.end()) {
      // key exist.
      if (it.first == "non1f") non1f[simulation_codon_3_letters] = it.second;
      if (it.first == "near1f") near1f[simulation_codon_3_letters] = it.second;
      if (it.first == "wobble1f") wobble1f[simulation_codon_3_letters] = it.second;
      if (it.first == "WC1f") WC1f[simulation_codon_3_letters] = it.second;
      if (it.first == "non1r") non1r = it.second;
      if (it.first == "near1r") near1r = it.second;
      if (it.first == "near2f") near2f = it.second;
      if (it.first == "near2r") near2r = it.second;
      if (it.first == "near3f") near3f = it.second;
      if (it.first == "near4f") near4f = it.second;
      if (it.first == "near5f") near5f = it.second;
      if (it.first == "neardiss") neardiss = it.second;
      if (it.first == "near6f") near6f = it.second;
      if (it.first == "wobble1r") wobble1r = it.second;
      if (it.first == "wobble2f") wobble2f = it.second;
      if (it.first == "wobble2r") wobble2r = it.second;
      if (it.first == "wobble3f") wobble3f = it.second;
      if (it.first == "wobble4f") wobble4f = it.second;
      if (it.first == "wobble5f") wobble5f = it.second;
      if (it.first == "wobblediss") wobblediss = it.second;
      if (it.first == "wobble6f") wobble6f = it.second;
      if (it.first == "WC1r") WC1r = it.second;
      if (it.first == "WC2f") WC2f = it.second;
      if (it.first == "WC2r") WC2r = it.second;
      if (it.first == "WC3f") WC3f = it.second;
      if (it.first == "WC4f") WC4f = it.second;
      if (it.first == "WC5f") WC5f = it.second;
      if (it.first == "WCdiss") WCdiss = it.second;
      if (it.first == "WC6f") WC6f = it.second;
      if (it.first == "dec7f") dec7f = it.second;
      if (it.first == "trans1f") trans1f = it.second;
      if (it.first == "trans1r") trans1r = it.second;
      if (it.first == "trans2") trans2 = it.second;
      if (it.first == "trans3") trans3 = it.second;
      if (it.first == "trans4") trans4 = it.second;
      if (it.first == "trans5") trans5 = it.second;
      if (it.first == "trans6") trans6 = it.second;
      if (it.first == "trans7") trans7 = it.second;
      if (it.first == "trans8") trans8 = it.second;
      if (it.first == "trans9") trans9 = it.second;
    }
  }
}

void Simulations::RibosomeSimulator::setPropensity(std::string &reaction,
                                                   const double propensity)
{
  *propensities_map.at(reaction) = propensity;
}

void Simulations::RibosomeSimulator::setNonCognate(double noNonCog)
{
  non1f[simulation_codon_3_letters] = noNonCog;
}

double Simulations::RibosomeSimulator::getPropensity(std::string reaction)
{
  return *propensities_map.at(reaction);
}

std::map<std::string, double>
Simulations::RibosomeSimulator::getPropensities()
{
  std::map<std::string, double> result;
  std::vector<double> ks = {non1f[simulation_codon_3_letters],
                            near1f[simulation_codon_3_letters],
                            wobble1f[simulation_codon_3_letters],
                            WC1f[simulation_codon_3_letters],
                            non1r,
                            near1r,
                            near2f,
                            near2r,
                            near3f,
                            near4f,
                            near5f,
                            neardiss,
                            near6f,
                            wobble1r,
                            wobble2f,
                            wobble2r,
                            wobble3f,
                            wobble4f,
                            wobble5f,
                            wobblediss,
                            wobble6f,
                            WC1r,
                            WC2f,
                            WC2r,
                            WC3f,
                            WC4f,
                            WC5f,
                            WCdiss,
                            WC6f,
                            dec7f,
                            trans1f,
                            trans1r,
                            trans2,
                            trans3,
                            trans4,
                            trans5,
                            trans6,
                            trans7,
                            trans8,
                            trans9};

  for (std::size_t i = 0; i < ks.size(); i++)
  {
    result[reactions_identifiers[i]] = ks[i];
  }
  return result;
}

void Simulations::RibosomeSimulator::setCodonForSimulation(
    const std::string &codon)
{
  simulation_codon_3_letters = codon;
  reactions_graph = reactions_map.at(codon);
  // populate propensities map so we can change propensities later.
  propensities_map.clear();
  propensities_map.emplace("non1f", &non1f[codon]);
  propensities_map.emplace("near1f", &near1f[codon]);
  propensities_map.emplace("wobble1f", &wobble1f[codon]);
  propensities_map.emplace("WC1f", &WC1f[codon]);
  propensities_map.emplace("non1r", &non1r);
  propensities_map.emplace("near1r", &near1r);
  propensities_map.emplace("near2f", &near2f);
  propensities_map.emplace("near2r", &near2r);
  propensities_map.emplace("near3f", &near3f);
  propensities_map.emplace("near4f", &near4f);
  propensities_map.emplace("near5f", &near5f);
  propensities_map.emplace("neardiss", &neardiss);
  propensities_map.emplace("near6f", &near6f);
  propensities_map.emplace("wobble1r", &wobble1r);
  propensities_map.emplace("wobble2f", &wobble2f);
  propensities_map.emplace("wobble2r", &wobble2r);
  propensities_map.emplace("wobble3f", &wobble3f);
  propensities_map.emplace("wobble4f", &wobble4f);
  propensities_map.emplace("wobble5f", &wobble5f);
  propensities_map.emplace("wobblediss", &wobblediss);
  propensities_map.emplace("wobble6f", &wobble6f);
  propensities_map.emplace("WC1r", &WC1r);
  propensities_map.emplace("WC2f", &WC2f);
  propensities_map.emplace("WC2r", &WC2r);
  propensities_map.emplace("WC3f", &WC3f);
  propensities_map.emplace("WC4f", &WC4f);
  propensities_map.emplace("WC5f", &WC5f);
  propensities_map.emplace("WCdiss", &WCdiss);
  propensities_map.emplace("WC6f", &WC6f);
  propensities_map.emplace("dec7f", &dec7f);
  propensities_map.emplace("trans1f", &trans1f);
  propensities_map.emplace("trans1r", &trans1r);
  propensities_map.emplace("trans2", &trans2);
  propensities_map.emplace("trans3", &trans3);
  propensities_map.emplace("trans4", &trans4);
  propensities_map.emplace("trans5", &trans5);
  propensities_map.emplace("trans6", &trans6);
  propensities_map.emplace("trans7", &trans7);
  propensities_map.emplace("trans8", &trans8);
  propensities_map.emplace("trans9", &trans9);
}

void Simulations::RibosomeSimulator::run_and_get_times(
    double &decoding_time, double &translocation_time)
{
  dt_history.clear();
  ribosome_state_history.clear();
  current_state = 0;

  // initialize the random generator
  std::random_device
      rd;                 // Will be used to obtain a seed for the random number engine
  std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
  std::uniform_real_distribution<> dis(0, 1);

  double r1 = 0, r2 = 0;
  double tau = 0;
  std::vector<double> alphas;
  std::vector<int> next_state;
  while (true)
  {
    // update history
    dt_history.push_back(tau);
    ribosome_state_history.push_back(getState());
    // randomly generate parameter for calculating dt
    r1 = dis(gen) + DBL_MIN; // adding minumum double value in order to avoid
                             // division by zero and infinities.
    // randomly generate parameter for selecting reaction
    r2 = dis(gen) + DBL_MIN; // adding minumum double value in order to avoid
                             // division by zero and infinities.
    // calculate an
    getAlphas(alphas, next_state);
    if (alphas.empty())
    {
      translocation_time = 0;
      decoding_time = 0;
      // no available reactions, get times and quit.
      bool is_translocating = true;
      for (int i = static_cast<int>(ribosome_state_history.size() - 1); i >= 0;
           i--)
      {
        if (is_translocating)
        {
          translocation_time += dt_history[static_cast<std::size_t>(i)];
          if (ribosome_state_history[static_cast<std::size_t>(i)] < 23)
          {
            is_translocating = false;
          }
        }
        else
        {
          decoding_time += dt_history[static_cast<std::size_t>(i)];
        }
      }
      // however, if quitting because there is no more reactions, time is infinity
      if (ribosome_state_history.size() == 0)
      {
        translocation_time = std::numeric_limits<double>::infinity();
        decoding_time = std::numeric_limits<double>::infinity();
      }
      return;
    }
    double a0 = std::accumulate(alphas.begin(), alphas.end(), 0.0);
    // select next reaction to execute
    double cumsum = 0;
    int selected_alpha_vector_index = -1;
    // TODO(Heday): vectorization of this loop would increase performance
    do
    {
      selected_alpha_vector_index++;
      cumsum += alphas[static_cast<std::size_t>(selected_alpha_vector_index)];
    } while (cumsum < a0 * r2);
    // Apply reaction
    setState(next_state[static_cast<std::size_t>(selected_alpha_vector_index)]);
    // Update time
    tau = (1.0 / a0) * log(1.0 / r1);
  }
}

double Simulations::RibosomeSimulator::run_repeatedly_get_average_time(const int repetitions)
{

  double r1 = 0.0, r2 = 0.0, a0 = 0.0;
  double cumsum = 0.0;
  double tau = 0.0, clock = 0.0;
  std::size_t selected_alpha_vector_index = 0;
  std::array<double, 4> alphas;
  std::array<int, 4> next_state;
  double k;      // get alphas
  int index, ii; //get alphas
  for (std::size_t i = 0; i < static_cast<std::size_t>(repetitions); i++)
  {
    current_state = 0;
    while (current_state < 32)
    {
      // randomly generate parameter for calculating dt
      r1 = dis(gen) + DBL_MIN; // adding minumum double value in order to avoid
                               // division by zero and infinities.
      // randomly generate parameter for selecting reaction
      r2 = dis(gen) + DBL_MIN; // adding minumum double value in order to avoid
                               // division by zero and infinities.
      // calculate an
      auto &alphas_and_indexes = reactions_graph[static_cast<std::size_t>(
          current_state)]; // go the possible
                           // reactions of that
                           // state.
      a0 = 0;
      ii = 0;
      for (auto element : alphas_and_indexes)
      {
        std::tie(k, index) = element;
        a0 += k;
        alphas[ii] = k;
        next_state[ii++] = index;
      }

      if (ii == 0)
        break;
      // a0 = std::accumulate(alphas.begin(), alphas.end(), 0.0);
      // select next reaction to execute
      cumsum = 0;
      selected_alpha_vector_index = -1;
      // TODO(Heday): vectorization of this loop would increase performance
      do
      {
        // selected_alpha_vector_index++;
        cumsum += alphas[selected_alpha_vector_index++];
      } while (cumsum < a0 * r2);
      // Apply reaction
      setState(next_state[selected_alpha_vector_index - 1]);
      // Update time
      tau = (1.0 / a0) * log(1.0 / r1);
      clock += tau;
    }
  }
  return clock / repetitions;
}

std::vector<std::vector<std::tuple<std::reference_wrapper<double>, int>>>
Simulations::RibosomeSimulator::createReactionsGraph(
    const csv_utils::concentration_entry &codon)
{
  std::array<std::reference_wrapper<double>, 40> ks = {{non1f[codon.codon],
                                                        near1f[codon.codon],
                                                        wobble1f[codon.codon],
                                                        WC1f[codon.codon],
                                                        non1r,
                                                        near1r,
                                                        near2f,
                                                        near2r,
                                                        near3f,
                                                        near4f,
                                                        near5f,
                                                        neardiss,
                                                        near6f,
                                                        wobble1r,
                                                        wobble2f,
                                                        wobble2r,
                                                        wobble3f,
                                                        wobble4f,
                                                        wobble5f,
                                                        wobblediss,
                                                        wobble6f,
                                                        WC1r,
                                                        WC2f,
                                                        WC2r,
                                                        WC3f,
                                                        WC4f,
                                                        WC5f,
                                                        WCdiss,
                                                        WC6f,
                                                        dec7f,
                                                        trans1f,
                                                        trans1r,
                                                        trans2,
                                                        trans3,
                                                        trans4,
                                                        trans5,
                                                        trans6,
                                                        trans7,
                                                        trans8,
                                                        trans9}};

  Eigen::MatrixXi reactionMatrix[40];
  // build the vector of reactions.
  // [] x=0 -> non1f:(x'=1);
  reactionMatrix[0].resize(32, 1);
  reactionMatrix[0].fill(0);
  reactionMatrix[0](0, 0) = -1;
  reactionMatrix[0](1, 0) = 1;

  // [] x=0 -> near1f:(x'=2);
  reactionMatrix[1].resize(32, 1);
  reactionMatrix[1].fill(0);
  reactionMatrix[1](0, 0) = -1;
  reactionMatrix[1](2, 0) = 1;

  // [] x=0 -> wobble1f:(x'=9);
  reactionMatrix[2].resize(32, 1);
  reactionMatrix[2].fill(0);
  reactionMatrix[2](0, 0) = -1;
  reactionMatrix[2](9, 0) = 1;

  // [] x=0 -> WC1f:(x'=16);
  reactionMatrix[3].resize(32, 1);
  reactionMatrix[3].fill(0);
  reactionMatrix[3](0, 0) = -1;
  reactionMatrix[3](16, 0) = 1;

  // [] x=1 -> non1r:(x'=0);
  reactionMatrix[4].resize(32, 1);
  reactionMatrix[4].fill(0);
  reactionMatrix[4](1, 0) = -1;
  reactionMatrix[4](0, 0) = 1;

  // [] x=2 -> near1r:(x'=0);
  reactionMatrix[5].resize(32, 1);
  reactionMatrix[5].fill(0);
  reactionMatrix[5](2, 0) = -1;
  reactionMatrix[5](0, 0) = 1;

  // [] x=2 -> near2f:(x'=3);
  reactionMatrix[6].resize(32, 1);
  reactionMatrix[6].fill(0);
  reactionMatrix[6](2, 0) = -1;
  reactionMatrix[6](3, 0) = 1;

  // [] x=3 -> near2r:(x'=2);
  reactionMatrix[7].resize(32, 1);
  reactionMatrix[7].fill(0);
  reactionMatrix[7](3, 0) = -1;
  reactionMatrix[7](2, 0) = 1;

  // [] x=3 -> near3f:(x'=4);
  reactionMatrix[8].resize(32, 1);
  reactionMatrix[8].fill(0);
  reactionMatrix[8](3, 0) = -1;
  reactionMatrix[8](4, 0) = 1;

  // [] x=4 -> near4f:(x'=5);
  reactionMatrix[9].resize(32, 1);
  reactionMatrix[9].fill(0);
  reactionMatrix[9](4, 0) = -1;
  reactionMatrix[9](5, 0) = 1;

  // [] x=5 -> near5f:(x'=6);
  reactionMatrix[10].resize(32, 1);
  reactionMatrix[10].fill(0);
  reactionMatrix[10](5, 0) = -1;
  reactionMatrix[10](6, 0) = 1;

  // [] x=6 -> neardiss:(x'=0);
  reactionMatrix[11].resize(32, 1);
  reactionMatrix[11].fill(0);
  reactionMatrix[11](6, 0) = -1;
  reactionMatrix[11](0, 0) = 1;

  //   // [] x=6 -> near6f:(x'=7);
  //   reactionMatrix[12].resize(32, 1);
  //   reactionMatrix[12].fill(0);
  //   reactionMatrix[12](6, 0) = -1;
  //   reactionMatrix[12](7, 0) = 1;

  //   // [] x=6 -> dec7f:(x'=21);
  reactionMatrix[12].resize(32, 1);
  reactionMatrix[12].fill(0);
  reactionMatrix[12](6, 0) = -1;
  reactionMatrix[12](21, 0) = 1;

  // [] x=7 -> near7f:(x'=8);
  //   reactionMatrix[13].resize(32, 1);
  //   reactionMatrix[13].fill(0);
  //   reactionMatrix[13](7, 0) = -1;
  //   reactionMatrix[13](8, 0) = 1;

  // [] x=8 -> trans1f:(x'=23);
  //   reactionMatrix[14].resize(32, 1);
  //   reactionMatrix[14].fill(0);
  //   reactionMatrix[14](8, 0) = -1;
  //   reactionMatrix[14](23, 0) = 1;

  // [] x=9 -> wobble1r:(x'=0);
  reactionMatrix[13].resize(32, 1);
  reactionMatrix[13].fill(0);
  reactionMatrix[13](9, 0) = -1;
  reactionMatrix[13](0, 0) = 1;

  // [] x=9 -> wobble2f:(x'=10);
  reactionMatrix[14].resize(32, 1);
  reactionMatrix[14].fill(0);
  reactionMatrix[14](9, 0) = -1;
  reactionMatrix[14](10, 0) = 1;

  // [] x=10 -> wobble2r:(x'=9);
  reactionMatrix[15].resize(32, 1);
  reactionMatrix[15].fill(0);
  reactionMatrix[15](10, 0) = -1;
  reactionMatrix[15](9, 0) = 1;

  // [] x=10 -> wobble3f:(x'=11);
  reactionMatrix[16].resize(32, 1);
  reactionMatrix[16].fill(0);
  reactionMatrix[16](10, 0) = -1;
  reactionMatrix[16](11, 0) = 1;

  // [] x=11 -> wobble4f:(x'=12);
  reactionMatrix[17].resize(32, 1);
  reactionMatrix[17].fill(0);
  reactionMatrix[17](11, 0) = -1;
  reactionMatrix[17](12, 0) = 1;

  // [] x=12 -> wobble5f:(x'=13);
  reactionMatrix[18].resize(32, 1);
  reactionMatrix[18].fill(0);
  reactionMatrix[18](12, 0) = -1;
  reactionMatrix[18](13, 0) = 1;

  // [] x=13 -> wobblediss:(x'=0);
  reactionMatrix[19].resize(32, 1);
  reactionMatrix[19].fill(0);
  reactionMatrix[19](13, 0) = -1;
  reactionMatrix[19](0, 0) = 1;

  //   // [] x=13 -> wobble6f:(x'=14);
  //   reactionMatrix[22].resize(32, 1);
  //   reactionMatrix[22].fill(0);
  //   reactionMatrix[22](13, 0) = -1;
  //   reactionMatrix[22](14, 0) = 1;

  // [] x=13 -> dec7f:(x'=21);
  reactionMatrix[20].resize(32, 1);
  reactionMatrix[20].fill(0);
  reactionMatrix[20](13, 0) = -1;
  reactionMatrix[20](21, 0) = 1;

  //   // [] x=14 -> wobble7f:(x'=15);
  //   reactionMatrix[23].resize(32, 1);
  //   reactionMatrix[23].fill(0);
  //   reactionMatrix[23](14, 0) = -1;
  //   reactionMatrix[23](15, 0) = 1;

  //   // [] x=15 -> trans1f:(x'=23);
  //   reactionMatrix[24].resize(32, 1);
  //   reactionMatrix[24].fill(0);
  //   reactionMatrix[24](15, 0) = -1;
  //   reactionMatrix[24](23, 0) = 1;

  // [] x=16 -> WC1r:(x'=0);
  reactionMatrix[21].resize(32, 1);
  reactionMatrix[21].fill(0);
  reactionMatrix[21](16, 0) = -1;
  reactionMatrix[21](0, 0) = 1;

  // [] x=16 -> WC2f:(x'=17);
  reactionMatrix[22].resize(32, 1);
  reactionMatrix[22].fill(0);
  reactionMatrix[22](16, 0) = -1;
  reactionMatrix[22](17, 0) = 1;

  // [] x=17 -> WC2r:(x'=16);
  reactionMatrix[23].resize(32, 1);
  reactionMatrix[23].fill(0);
  reactionMatrix[23](17, 0) = -1;
  reactionMatrix[23](16, 0) = 1;

  // [] x=17 -> WC3f:(x'=18);
  reactionMatrix[24].resize(32, 1);
  reactionMatrix[24].fill(0);
  reactionMatrix[24](17, 0) = -1;
  reactionMatrix[24](18, 0) = 1;

  // [] x=18 -> WC4f:(x'=19);
  reactionMatrix[25].resize(32, 1);
  reactionMatrix[25].fill(0);
  reactionMatrix[25](18, 0) = -1;
  reactionMatrix[25](19, 0) = 1;

  // [] x=19 -> WC5f:(x'=20);
  reactionMatrix[26].resize(32, 1);
  reactionMatrix[26].fill(0);
  reactionMatrix[26](19, 0) = -1;
  reactionMatrix[26](20, 0) = 1;

  // [] x=20 -> WCdiss:(x'=0);
  reactionMatrix[27].resize(32, 1);
  reactionMatrix[27].fill(0);
  reactionMatrix[27](20, 0) = -1;
  reactionMatrix[27](0, 0) = 1;

  // [] x=20 -> WC6f:(x'=21);
  reactionMatrix[28].resize(32, 1);
  reactionMatrix[28].fill(0);
  reactionMatrix[28](20, 0) = -1;
  reactionMatrix[28](21, 0) = 1;

  // [] x=21 -> WC7f:(x'=22);
  reactionMatrix[29].resize(32, 1);
  reactionMatrix[29].fill(0);
  reactionMatrix[29](21, 0) = -1;
  reactionMatrix[29](22, 0) = 1;

  // [] x=22 -> trans1f:(x'=23);
  reactionMatrix[30].resize(32, 1);
  reactionMatrix[30].fill(0);
  reactionMatrix[30](22, 0) = -1;
  reactionMatrix[30](23, 0) = 1;

  // [] x=23 -> trans1r:(x'=22);
  reactionMatrix[31].resize(32, 1);
  reactionMatrix[31].fill(0);
  reactionMatrix[31](23, 0) = -1;
  reactionMatrix[31](22, 0) = 1;

  // [] x=23 -> trans2:(x'=24);
  reactionMatrix[32].resize(32, 1);
  reactionMatrix[32].fill(0);
  reactionMatrix[32](23, 0) = -1;
  reactionMatrix[32](24, 0) = 1;

  // [] x=24 -> trans3:(x'=25);
  reactionMatrix[33].resize(32, 1);
  reactionMatrix[33].fill(0);
  reactionMatrix[33](24, 0) = -1;
  reactionMatrix[33](25, 0) = 1;

  // [] x=25 -> trans4:(x'=26);
  reactionMatrix[34].resize(32, 1);
  reactionMatrix[34].fill(0);
  reactionMatrix[34](25, 0) = -1;
  reactionMatrix[34](26, 0) = 1;

  // [] x=26 -> trans5:(x'=27);
  reactionMatrix[35].resize(32, 1);
  reactionMatrix[35].fill(0);
  reactionMatrix[35](26, 0) = -1;
  reactionMatrix[35](27, 0) = 1;

  // [] x=27 -> trans6:(x'=28);
  reactionMatrix[36].resize(32, 1);
  reactionMatrix[36].fill(0);
  reactionMatrix[36](27, 0) = -1;
  reactionMatrix[36](28, 0) = 1;

  // [] x=28 -> trans7:(x'=29);
  reactionMatrix[37].resize(32, 1);
  reactionMatrix[37].fill(0);
  reactionMatrix[37](28, 0) = -1;
  reactionMatrix[37](29, 0) = 1;

  // [] x=29 -> trans8:(x'=30);
  reactionMatrix[38].resize(32, 1);
  reactionMatrix[38].fill(0);
  reactionMatrix[38](29, 0) = -1;
  reactionMatrix[38](30, 0) = 1;

  // [] x=30 -> trans9:(x'=31);
  reactionMatrix[39].resize(32, 1);
  reactionMatrix[39].fill(0);
  reactionMatrix[39](30, 0) = -1;
  reactionMatrix[39](31, 0) = 1;

  int ii = 0;
  std::vector<std::vector<std::tuple<std::reference_wrapper<double>, int>>> r_g;
  r_g.resize(32);
  std::fill(r_g.begin(), r_g.end(),
            std::vector<std::tuple<std::reference_wrapper<double>, int>>());
  // the vector reactions_graph (I know, not a good name. needs to be changed
  // at some point.), have the following format: reactions_graph[current
  // ribisome state] = [vector of tuples(reaction propensity, ribosome state)]
  // this way, if the ribosome state is, say, 0, we check the possible
  // reactions at reactions_graph[0]. if, say we select the reaction with the
  // tuple (0.3, 16), it means that the reaction propensity is 0.3 and it will
  // make the ribosome state go to 16. This is purely for the sake of
  // optimization. the loop below populates reactions_graph automatically. It
  // assumes that each reaction is first-degree.
  for (const Eigen::MatrixXi &m : reactionMatrix)
  {
    if (ks.at(static_cast<std::size_t>(ii)) > 0)
    {
      // populate the local index.
      Eigen::Index maxRow, maxCol, minRow, minCol;
      m.maxCoeff(&maxRow, &maxCol); // 1
      m.minCoeff(&minRow, &minCol); // -1
      if (ks.at(static_cast<std::size_t>(ii)) > 0)
      {
        r_g.at(static_cast<std::size_t>(minRow))
            .push_back({ks.at(static_cast<std::size_t>(ii)), maxRow});
      }
    }
    ii++;
  }
  return r_g;
}

int Simulations::RibosomeSimulator::getState() { return current_state; }
void Simulations::RibosomeSimulator::setState(int s) { current_state = s; }

void Simulations::RibosomeSimulator::getAlphas(
    std::vector<double> &as, std::vector<int> &reactions_index)
{
  as.clear();
  reactions_index.clear();
  auto alphas_and_indexes = reactions_graph[static_cast<std::size_t>(
      current_state)]; // go the possible
                       // reactions of that
                       // state.
  double k;
  int index;
  for (auto element : alphas_and_indexes)
  {
    std::tie(k, index) = element;
    as.push_back(k);
    reactions_index.push_back(index);
  }
}

void Simulations::RibosomeSimulator::getDecodingAlphas(
    std::vector<double> &as, std::vector<int> &reactions_index)
{
  as.clear();
  reactions_index.clear();
  auto alphas_and_indexes = reactions_graph[static_cast<std::size_t>(
      current_state)]; // go the possible
                       // reactions of that
                       // state.
  double k;
  int index;
  for (auto element : alphas_and_indexes)
  {
    std::tie(k, index) = element;
    if (index < 23)
    {
      as.push_back(k);
      reactions_index.push_back(index);
    }
  }
}
