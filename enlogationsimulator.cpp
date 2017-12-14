#ifdef COMIPLE_PYTHON_MODULE
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

#endif

#include <algorithm>
#include <fstream>
#include <iostream>
#include <list>
#include <numeric>
#include <set>
#include <tuple>
#include "enlogationsimulator.h"

using namespace Simulations;

#ifdef COMIPLE_PYTHON_MODULE
PYBIND11_PLUGIN(enlogationsimulator) {
  pybind11::module mod("enlogationsimulator", "auto-compiled c++ extension");

  py::class_<Gillespie>(mod, "gillespie")
      .def(py::init<>())  // constructor
      .def("setIterationLimit", &Gillespie::setIterationLimit)
      .def("run", &Gillespie::run);

  py::class_<EnlogationSimulator, Gillespie>(mod, "enlogationsimulator")
      .def(py::init<>())  // constructor
      .def("setAverageTimesFileName",
           &EnlogationSimulator::setAverageTimesFileName)
      .def("setConcentrationsFileName",
           &EnlogationSimulator::setConcentrationsFileName)
      .def("setMRnaFileName", &EnlogationSimulator::setMRnaFileName)
      .def("setInitiationRate", &EnlogationSimulator::setInitiationRate)
      .def("setTimeLimit", &Gillespie::setTimeLimit)
      .def("setTerminationRate", &EnlogationSimulator::setTerminationRate)
      .def("setIterationLimit", &EnlogationSimulator::setIterationLimit)
      .def("updateRibosomeHistory", &EnlogationSimulator::updateRibosomeHistory,
           py::arg("clear_population_history") = false)
      .def("getEnlogationDuration", &EnlogationSimulator::getEnlogationDuration)
      .def("calculateAverageTimes", &EnlogationSimulator::calculateAverageTimes)
      .def_readonly("ribosome_positions_history",
                    &EnlogationSimulator::ribosome_positions_history)
      .def_readonly("dt_history", &EnlogationSimulator::dt_history)
      .def_readonly("total_time", &EnlogationSimulator::total_time)
      .def_readonly("n_times_occupied", &EnlogationSimulator::n_times_occupied)
      .def_readonly("average_times",
                    &EnlogationSimulator::codons_average_occupation_time);

  return mod.ptr();
}
#endif

EnlogationSimulator::EnlogationSimulator() {}

void EnlogationSimulator::setInitiationRate(double ir) {
  if (ir > 0) {
    initiation_rate = ir;
  } else {
    throw std::runtime_error("invalid initiation rate: " + std::to_string(ir));
  }
  intializeMRNAReader();
}

void EnlogationSimulator::setTerminationRate(double tr) {
  if (tr > 0) {
    termination_rate = tr;
  } else {
    throw std::runtime_error("invalid termination rate: " + std::to_string(tr));
  }
  intializeMRNAReader();
}

void EnlogationSimulator::setMRnaFileName(std::string file_name) {
  std::ifstream ist{file_name};

  if (!ist) {
    throw std::runtime_error("can't open input file: " + file_name);
  } else {
    mRNA_file_name = file_name;
  }
  intializeMRNAReader();
}

void EnlogationSimulator::intializeMRNAReader() {
  if (!mRNA_file_name.empty() and !average_times_file_name.empty() &&
      initiation_rate > 0 && termination_rate > 0) {
    // we can proceed with the mRNAReader object.
    //    mrna_reader.loadRateCalculatorFile(average_times_file_name);
    mrna_reader.loadmRNAFile(mRNA_file_name);
    mrna_reader.setInitiationRate(initiation_rate);
    mrna_reader.setTerminationRate(termination_rate);
    mrna_reader.generateInitialPopulation();
    mrna_reader.generateReactions();
    Gillespie::setInitialPopulation(mrna_reader.initial_population);
    Gillespie::setReactionsSet(mrna_reader.reactions_set);
  }
}

void EnlogationSimulator::setAverageTimesFileName(std::string file_name) {
  std::ifstream ist{file_name};

  if (!ist) {
    throw std::runtime_error("can't open input file: " + file_name);
  } else {
    average_times_file_name = file_name;
  }
  intializeMRNAReader();
}

void EnlogationSimulator::setConcentrationsFileName(std::string file_name) {
  std::ifstream ist{file_name};

  if (!ist) {
    throw std::runtime_error("can't open input file: " + file_name);
  } else {
    concentrations_file_name = file_name;
    // when setting the concentrations file name, we can also
    // initialize the RibosomeSimulator object.
    ribosome_simulator.loadConcentrations(file_name);
  }
}

double EnlogationSimulator::getReactionTime(double& a0, double& r1,
                                            std::string& codon) {
  double result = 0;
  if (codon == "tra") {
    if (translocation_times.empty()) {
      // calculate the time.
      double translocating;
      ribosome_simulator.setCodonForSimulation("AUG");
      ribosome_simulator.run_and_get_times(result, translocating);
      // add the translocating time to our pool.
      translocation_times.push_back(translocating);
    }
    result = translocation_times.back();
    translocation_times.pop_back();
  } else if (codon == "ter" || codon == "ini" ||
             std::find(stop_codons.begin(), stop_codons.end(), codon) !=
                 end(stop_codons)) {
    result = Gillespie::getReactionTime(a0, r1, codon);
  } else {
    // calculate the time.
    double translocating;
    ribosome_simulator.setCodonForSimulation(codon);
    ribosome_simulator.run_and_get_times(result, translocating);
    // add the translocating time to our pool.
    translocation_times.push_back(translocating);
  }
  return result;
}

/**
 * @brief reads the population_history to calculate the ribosome positions, then
 * populate the positions_set vector with the ribosome positions (in codon
 * number, starting from 0) in each iteration.
 *
 * @param clear_population_history p_clear_population_history:...if true, erase
 * the population_history vector. This might be useful in simulations with large
 * number of iterations and/or large mRNA sequences, since these factors will
 * considerably increase the memory usage of this object. This options is for
 * memory optimization only.
 */
void EnlogationSimulator::updateRibosomeHistory(bool clear_population_history) {
  ribosome_positions_history.clear();
  for (Eigen::MatrixXi population : population_history) {
    std::set<int> positions_set;
    std::vector<int> positions_vector;
    for (int i = 0; i < population.cols(); i++) {
      if (population(2, i) > 0 || population(3, i) > 0) {
        positions_set.insert(i);
      }
    }
    positions_vector.clear();
    positions_vector.insert(positions_vector.end(), positions_set.begin(),
                            positions_set.end());
    ribosome_positions_history.push_back(positions_vector);
  }
  if (clear_population_history) population_history.clear();
}

/**
 * @brief Returns a tuple where the first element is a vector with the
 * enlogation duration of the ribosomes that terminated in the simulation, and
 * the second element is a vector with the iteration where such ribosomes
 * started enlogating. This method should be called after updateRibosomeHistory,
 * since it uses the positions_vector to do its job.
 *
 */

std::tuple<std::vector<double>, std::vector<int>, std::vector<int>>
EnlogationSimulator::getEnlogationDuration() {
  std::list<int> rib_initiation_iteration;
  std::vector<double> result;
  std::vector<int> initiation_iteration, termination_iteration;
  std::vector<double> clock;
  clock.clear();
  result.clear();
  termination_iteration.clear();
  initiation_iteration.clear();
  double cumsum = 0;
  for (unsigned int i = 0; i < dt_history.size(); i++) {
    cumsum += dt_history[i];
    clock.push_back(cumsum);
  }
  bool previous_zero_occupied = false, previous_last_occupied = false;
  bool current_zero_occupied = false, current_last_occupied = false;
  int last_codon_position =
      static_cast<int>((mrna_reader.mRNA_sequence.size() / 3) - 1);
  for (int i = 1;
       static_cast<std::size_t>(i) < ribosome_positions_history.size(); i++) {
    std::vector<int>& ribosomes_positions =
        ribosome_positions_history[static_cast<std::size_t>(i)];
    current_zero_occupied =
        std::find(ribosomes_positions.begin(), ribosomes_positions.end(), 0) !=
        ribosomes_positions.end();
    current_last_occupied =
        std::find(ribosomes_positions.begin(), ribosomes_positions.end(),
                  last_codon_position) != ribosomes_positions.end();
    if (!previous_zero_occupied && current_zero_occupied) {
      // new ribosome.
      rib_initiation_iteration.push_back(i);
    } else if (previous_last_occupied && !current_last_occupied) {
      // ribosome left.
      result.push_back(
          clock[static_cast<std::size_t>(i - 1)] -
          clock[static_cast<std::size_t>(rib_initiation_iteration.front())]);
      initiation_iteration.push_back(rib_initiation_iteration.front());
      termination_iteration.push_back(i - 1);
      rib_initiation_iteration.pop_front();
    }
    // update variables for next iteration.
    previous_last_occupied = current_last_occupied;
    previous_zero_occupied = current_zero_occupied;
  }
  return std::make_tuple(result, initiation_iteration, termination_iteration);
}

void EnlogationSimulator::calculateAverageTimes() {
  std::size_t number_codons = mrna_reader.mRNA_sequence.size() / 3;
  // initialize the total_time vector.
  total_time = std::vector<double>(number_codons);
  std::fill(total_time.begin(), total_time.end(), 0);
  // initialize the n_times_occupied vector.
  n_times_occupied = std::vector<int>(number_codons);
  std::fill(n_times_occupied.begin(), n_times_occupied.end(), 0);
  // iteration where we last seen the codon being occupied.
  std::vector<int> last_index_occupied(number_codons);
  std::fill(last_index_occupied.begin(), last_index_occupied.end(), -1);
  int iteration_number = 0;
  for (std::vector<int> ribosome_vector : ribosome_positions_history) {
    for (int position : ribosome_vector) {
      total_time[static_cast<std::size_t>(position)] +=
          dt_history[static_cast<std::size_t>(iteration_number)];
      if (last_index_occupied[static_cast<std::size_t>(position)] == -1 ||
          last_index_occupied[static_cast<std::size_t>(position)] !=
              iteration_number - 1) {
        // we are facing a re-entering ribosome. We need to add the previous
        // occupation.
        n_times_occupied[static_cast<std::size_t>(position)]++;
      }
      // update the last time this position was occupied.
      last_index_occupied[static_cast<std::size_t>(position)] =
          iteration_number;
    }
    iteration_number++;
  }
  // the above procedure does not count for the last time a position has been
  // occupied: it ignores it.  we could try to fix this in a number of ways, but
  // I guess it wouldn't matter much in the big picture.

  // now we calculate the averages.
  codons_average_occupation_time.clear();
  codons_average_occupation_time = std::vector<double>(number_codons);
  for (std::size_t codon_position = 0; codon_position < number_codons;
       codon_position++) {
    codons_average_occupation_time[codon_position] =
        n_times_occupied[codon_position] > 0
            ? total_time[codon_position] / n_times_occupied[codon_position]
            : 0;
  }
}
