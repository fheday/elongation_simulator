#include "translation.h"

#include <float.h>
#include <algorithm>
#include <deque>
#include <fstream>
#include <iostream>
#include <list>
#include <random>
#include "initiationterminationcodon.h"
#include "mrna_reader.h"

#define RIBOSOME_SIZE 10

#ifdef COMIPLE_PYTHON_MODULE

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

PYBIND11_MODULE(translation, mod) {
  py::class_<Simulations::Translation>(mod, "translation")
      .def(py::init<>())  // constructor
      .def("loadMRNA", &Simulations::Translation::loadMRNA)
      .def("loadConcentrations", &Simulations::Translation::loadConcentrations)
      .def("setInitiationRate", &Simulations::Translation::setInitiationRate)
      .def("setTerminationRate", &Simulations::Translation::setTerminationRate)
      .def("setIterationLimit", &Simulations::Translation::setIterationLimit)
      .def("setTimeLimit", &Simulations::Translation::setTimeLimit)
      .def("setIterationLimit", &Simulations::Translation::setIterationLimit)
      .def("setFinishedRibosomes",
           &Simulations::Translation::setFinishedRibosomes)
      .def("run", &Simulations::Translation::run,
           py::call_guard<py::gil_scoped_release>())
      .def("setInitiationRate", &Simulations::Translation::setInitiationRate)
      .def("setInitiationRate", &Simulations::Translation::setInitiationRate)
      .def("getEnlogationDuration",
           &Simulations::Translation::getEnlogationDuration)
      .def("calculateAverageTimes",
           &Simulations::Translation::calculateAverageTimes)
      .def("setPrepopulate", &Simulations::Translation::setPrepopulate)
      .def("getInitiationEnlongationTermination",
           &Simulations::Translation::getInitiationEnlongationTermination)
      .def("setLogCodonStates", &Simulations::Translation::setLogCodonStates)
      .def("getLogCodonStates", &Simulations::Translation::getLogCodonStates)
      .def("setWCPropensities", &Simulations::Translation::setWCPropensities)
      .def("setWooblePropensities",
           &Simulations::Translation::setWooblePropensities)
      .def("setNearCognatePropensities",
           &Simulations::Translation::setNearCognatePropensities)
      .def("setNonCogPropensities",
           &Simulations::Translation::setNonCogPropensities)
      .def("setTranslocationPropensities",
           &Simulations::Translation::setTranslocationPropensities)
      .def("getPropensities", &Simulations::Translation::getPropensities)

      .def_readonly("mrna_file_name", &Simulations::Translation::mrna_file_name)
      .def_readonly("concentrations_file_name",
                    &Simulations::Translation::concentrations_file_name)
      .def_readonly("dt_history", &Simulations::Translation::dt_history)
      .def_readonly("ribosome_positions_history",
                    &Simulations::Translation::ribosome_positions_history)
      .def_readonly("initiationRate",
                    &Simulations::Translation::initiation_rate)
      .def_readonly("terminationRate",
                    &Simulations::Translation::termination_rate)

      .def_readonly("initiations_durations",
                    &Simulations::Translation::initiations_durations)
      .def_readonly("enlongations_durations",
                    &Simulations::Translation::enlongations_durations)
      .def_readonly("terminations_durations",
                    &Simulations::Translation::terminations_durations)

      .def_readonly("total_time", &Simulations::Translation::total_time)
      .def_readonly("n_times_occupied",
                    &Simulations::Translation::n_times_occupied)
      .def_readonly("average_times",
                    &Simulations::Translation::codons_average_occupation_time);
}

#endif

void Simulations::Translation::loadConcentrations(
    const std::string& file_name) {
  std::ifstream ist{file_name};

  if (!ist) {
    throw std::runtime_error("can't open input file: " + file_name);
  } else {
    concentrations_file_name = file_name;
    initializeMRNAReader();
  }
}

void Simulations::Translation::loadMRNA(const std::string& file_name) {
  std::ifstream ist{file_name};

  if (!ist) {
    throw std::runtime_error("can't open input file: " + file_name);
  } else {
    mrna_file_name = file_name;
    initializeMRNAReader();
  }
}

void Simulations::Translation::initializeMRNAReader() {
  if (!concentrations_file_name.empty() && !mrna_file_name.empty() &&
      initiation_rate > 0 && termination_rate > 0) {
    // we have the concentrations and mrna file names. we can proceed.
    mRNA_utils::mRNAReader mrr;
    mrr.loadmRNAFile(mrna_file_name);
    // fill codon vector.
    int n_codons = mrr.sizeInCodons();
    std::unique_ptr<Simulations::InitiationTerminationCodon> initiation_codon(
        new Simulations::InitiationTerminationCodon(initiation_rate, true));
    initiation_codon->codon = mrr.getCodon(0);
    initiation_codon->index = 0;
    initiation_codon->setState(0);
    codons_vector.push_back(std::move(initiation_codon));
    for (int i = 1; i < n_codons - 1; i++) {
      std::unique_ptr<Simulations::EnlongationCodon> c(
          new Simulations::EnlongationCodon());
      c->loadConcentrations(concentrations_file_name);
      c->codon = mrr.getCodon(i);
      c->setCodon(c->codon);
      c->index = i;
      codons_vector.push_back(std::move(c));
    }

    // termination codon
    std::unique_ptr<Simulations::InitiationTerminationCodon> termination_codon(
        new Simulations::InitiationTerminationCodon(termination_rate, false));
    termination_codon->codon = mrr.getCodon(n_codons - 1);

    termination_codon->index = n_codons - 1;
    codons_vector.push_back(std::move(termination_codon));

    // link codons.
    for (unsigned int i = 1; i < codons_vector.size() - 1; i++) {
      codons_vector[i]->setNextCodon(codons_vector[i + 1].get());
      codons_vector[i]->setPreviousCodon(codons_vector[i - 1].get());
    }
    codons_vector[codons_vector.size() - 1]->setPreviousCodon(
        codons_vector[codons_vector.size() - 2].get());
    codons_vector[0]->setNextCodon(codons_vector[1].get());
  }
}

void Simulations::Translation::setWCPropensities(std::array<double, 10> prop) {
  for (std::size_t i = 1; i < codons_vector.size() - 1; i++) {
    codons_vector[i]->setWCPropensities(prop);
  }
}

void Simulations::Translation::setWooblePropensities(
    std::array<double, 10> prop) {
  for (std::size_t i = 1; i < codons_vector.size() - 1; i++) {
    codons_vector[i]->setWooblePropensities(prop);
  }
}

void Simulations::Translation::setNearCognatePropensities(
    std::array<double, 10> prop) {
  for (std::size_t i = 1; i < codons_vector.size() - 1; i++) {
    codons_vector[i]->setNearCognatePropensities(prop);
  }
}

void Simulations::Translation::setNonCogPropensities(
    std::array<double, 2> prop) {
  for (std::size_t i = 1; i < codons_vector.size() - 1; i++) {
    codons_vector[i]->setNonCogPropensities(prop);
  }
}

void Simulations::Translation::setTranslocationPropensities(
    std::array<double, 10> prop) {
  for (std::size_t i = 1; i < codons_vector.size() - 1; i++) {
    codons_vector[i]->setTranslocationPropensities(prop);
  }
}

std::vector<std::map<std::string, double>>
Simulations::Translation::getPropensities() {
  auto result = std::vector<std::map<std::string, double>>();
  result.push_back(std::map<std::string, double>());  // codon 0 will be empty.
  for (std::size_t i = 1; i < codons_vector.size() - 1; i++) {
    result.push_back(codons_vector[i]->getPropensities());
  }
  return result;
}

void Simulations::Translation::setInitiationRate(double ir) {
  if (ir > 0) {
    initiation_rate = ir;
  }
  initializeMRNAReader();
}

void Simulations::Translation::setTerminationRate(double tr) {
  if (tr > 0) {
    termination_rate = tr;
  }
  initializeMRNAReader();
}

void Simulations::Translation::setPrepopulate(bool prep) {
  pre_populate = prep;
}

/**
 * @brief Set a iteration limit for the Gillespie simulation.
 *
 * @param i integer with the maximum number of iterations. The algorithm halts
 * before this condition is met if there are no possible reations left to be
 * performed.
 */
void Simulations::Translation::setIterationLimit(int i) {
  if (i > 0) {
    iteration_limit = i;
  }
}

/**
 * @brief Set a time limit for the Gillespie simulation. This time is in
 * seconds, and it is compared against the simulation's clock.
 *
 * @param t time limit in seconds.
 */

void Simulations::Translation::setTimeLimit(double t) {
  if (t > 0) {
    time_limit = t;
  }
}

/**
 * @brief Set a limit of the number of ribosomes that sucessfully initiate and
 * terminates the mRNA.
 *
 * @param n_ribosomes p_n_ribosomes:The simulation will end after this number of
 * ribosomes terminates the mRNA.
 */
void Simulations::Translation::setFinishedRibosomes(int n_ribosomes) {
  if (n_ribosomes > 0) {
    finished_ribosomes_limit = n_ribosomes;
  }
}

void Simulations::Translation::getAlphas() {
  alphas.clear();
  codon_index.clear();
  reaction_index.clear();
  std::vector<double> a;
  std::vector<int> r_i;

  // populate the vectors.
  std::vector<int> ribosome_positions = ribosome_positions_history.back();
  std::size_t ribosome_index;
  // add initiation if needed.
  if (ribosome_positions.empty() ||
      (ribosome_positions[0] > (RIBOSOME_SIZE - 1) &&
       codons_vector[0]->isAvailable())) {
    // need to add initalization.
    codons_vector[0]->getAlphas(a, r_i);
    alphas.insert(alphas.end(), a.begin(), a.end());
    std::fill_n(std::back_inserter(codon_index), a.size(), 0);
    reaction_index.insert(reaction_index.end(), r_i.begin(), r_i.end());
  }

  for (unsigned i = 0; i < ribosome_positions.size(); i++) {
    ribosome_index = static_cast<std::size_t>(ribosome_positions[i]);
    codons_vector[ribosome_index]->getAlphas(a, r_i);
    alphas.insert(alphas.end(), a.begin(), a.end());
    std::fill_n(std::back_inserter(codon_index), a.size(), ribosome_index);
    reaction_index.insert(reaction_index.end(), r_i.begin(), r_i.end());
  }
}

void Simulations::Translation::run() {
  dt_history = std::vector<double>(
      iteration_limit > 0 ? static_cast<std::size_t>(iteration_limit) : 100000);
  dt_history.clear();

  ribosome_positions_history = std::vector<std::vector<int>>(
      iteration_limit > 0 ? static_cast<std::size_t>(iteration_limit) : 100000);
  ribosome_positions_history.clear();

  // initialize the random generator
  std::random_device
      rd;  // Will be used to obtain a seed for the random number engine
  std::mt19937 gen(rd());  // Standard mersenne_twister_engine seeded with rd()
  std::uniform_real_distribution<> dis(DBL_MIN, 1);

  double r1 = 0, r2 = 0;
  double tau = 0, clock = 0.0;
  int i = 0;

  int finished_ribosomes = 0, pre_filled_ribosomes = 0;
  std::vector<int> rib_positions((codons_vector.size() / RIBOSOME_SIZE) + 1);
  rib_positions.clear();
  // pre-fill codons based on the rates.
  if (pre_populate) {
    std::vector<double> a;
    std::vector<int> r_i;
    codons_vector[0]->getAlphas(a, r_i);
    double initiation_time = 1 / a[0];  // propensity
    std::size_t last_index = codons_vector.size() - 1;
    double time_sum = 0;
    codons_vector[last_index]->setOccupied(true);
    codons_vector[last_index]->setAvailable(false);
    codons_vector[last_index]->setState(0);
    pre_filled_ribosomes++;
    for (int i = static_cast<int>(codons_vector.size() - 2); i >= 0; i--) {
      if (last_index - static_cast<std::size_t>(i) <= RIBOSOME_SIZE - 1) {
        codons_vector[static_cast<std::size_t>(i)]->setAvailable(false);
      } else {
        // update ribosome.
        codons_vector[static_cast<std::size_t>(i)]->getAlphas(a, r_i);
        time_sum += 1 / a[0];
      }
      if (time_sum >= initiation_time) {
        // put a ribosome here.
        codons_vector[static_cast<std::size_t>(i)]->setOccupied(true);
        codons_vector[static_cast<std::size_t>(i)]->setAvailable(false);
        codons_vector[static_cast<std::size_t>(i)]->setState(0);
        if (i == 0) {
          codons_vector[static_cast<std::size_t>(i)]->setState(23);
        }
        time_sum = 0;  // reset timer.
        last_index = static_cast<std::size_t>(
            i);  // mark this as last inserted ribosome.
        pre_filled_ribosomes++;
      }
    }
  }
  for (unsigned int i = 0; i < codons_vector.size(); i++) {
    if (codons_vector[i]->isOccupied()) {
      rib_positions.push_back(static_cast<int>(i));
    }
  }
  finished_ribosomes -=
      pre_filled_ribosomes;  // we should ignore these ribosomes.
  int moved_codon = -1, current_codon = -1;
  bool initiation = false, termination = false, moved = true;

  while ((iteration_limit > 0 && i < iteration_limit) ||
         (time_limit > 0 && clock < time_limit) ||
         (finished_ribosomes_limit > 0 &&
          finished_ribosomes_limit > finished_ribosomes)) {
    // get the vector with the positions of all ribosomes
    if (i > 0 && moved) {
      if (termination) {
        // terminated. remove last position.
        rib_positions.pop_back();
      } else if (initiation) {
        // initiated.
        rib_positions.insert(rib_positions.begin(), 0);
      } else if (moved) {
        auto i = std::find(rib_positions.begin(), rib_positions.end(),
                           moved_codon - 1);
        auto pos = std::distance(rib_positions.begin(), i);
        rib_positions[static_cast<std::size_t>(pos)] = moved_codon;
      }
    }
    if (!moved) {
      // no ribosome movement. just update dt_history.
      dt_history.back() += tau;
    } else {
      // ribosome movement detected. create new entry in the history.
      dt_history.push_back(tau);
      ribosome_positions_history.push_back(rib_positions);
    }

    moved = false;
    initiation = false;
    termination = false;
    // randomly generate parameter for calculating dt
    r1 = dis(gen);
    // randomly generate parameter for selecting reaction
    r2 = dis(gen);
    // calculate an
    getAlphas();
    if (alphas.empty()) {
      // no available reactions, quit loop prematurely.
      std::cout << "no available reactions. quitting.\n";
      break;
    }
    double a0 = std::accumulate(alphas.begin(), alphas.end(), 0.0);
    int selected_alpha_vector_index = -1;
    // The commented code below is the vectorized version of the reaction
    // selection. upper_bound stops when it finds the first position that is
    // greater than the last parameter. we do an additional operation to find
    // the index of that position. as it seems so far, this is almost equivalent
    // in speed to the non-vectorized version.

    //         std::vector<double> cumsum(alphas.size());
    //         std::partial_sum(alphas.begin(), alphas.end(), cumsum.begin());
    //         selected_alpha_vector_index = std::distance(cumsum.begin(),
    //         std::upper_bound(cumsum.begin(), cumsum.end(), a0 * r2));

    // select next reaction to execute
    double cumsum = 0;
    do {
      selected_alpha_vector_index++;
      cumsum += alphas[static_cast<std::size_t>(selected_alpha_vector_index)];
    } while (cumsum < a0 * r2);
    current_codon =
        codon_index[static_cast<std::size_t>(selected_alpha_vector_index)];
    // Apply reaction
    codons_vector[static_cast<std::size_t>(current_codon)]->executeReaction(
        reaction_index[static_cast<std::size_t>(selected_alpha_vector_index)]);

    if (codon_index[static_cast<std::size_t>(selected_alpha_vector_index)] ==
            0 &&
        codons_vector[0]->getState() == 23) {
      // initiated.
      codons_vector[0]->setAvailable(false);
      codons_vector[0]->setOccupied(true);
      initiation = true;
      moved = true;
    }
    moved_codon = current_codon + 1;
    // 2- Any codon with state == 31 means the ribosome already moved to the
    // next codon (or left the mRNA). update states.
    if (codons_vector[static_cast<std::size_t>(current_codon)]->getState() ==
        31) {
      codons_vector[static_cast<std::size_t>(current_codon)]->setState(0);
      codons_vector[static_cast<std::size_t>(current_codon)]->setAvailable(
          false);
      codons_vector[static_cast<std::size_t>(current_codon)]->setOccupied(
          false);
      moved = true;
      moved_codon = current_codon + 1;
      if (static_cast<std::size_t>(moved_codon + 1) < codons_vector.size()) {
        codons_vector[static_cast<std::size_t>(moved_codon)]->setOccupied(true);
        codons_vector[static_cast<std::size_t>(moved_codon)]->setAvailable(
            false);
      }
      // update free codons due to the size of the ribosome.
      if (static_cast<std::size_t>(moved_codon) > RIBOSOME_SIZE - 1) {
        // we need to do some tidying up after the ribosome.
        if (static_cast<std::size_t>(moved_codon) < codons_vector.size()) {
          // update freed space left by the ribosome's movement.
          codons_vector[static_cast<std::size_t>(moved_codon - RIBOSOME_SIZE)]
              ->setAvailable(true);
        } else if (static_cast<std::size_t>(moved_codon) ==
                   codons_vector.size()) {
          // ribosome terminated. free codons positions occupied by it.
          termination = true;
          for (std::size_t i = static_cast<std::size_t>(codons_vector.size() -
                                                        RIBOSOME_SIZE),
                           total = codons_vector.size();
               i < total; ++i) {
            codons_vector[i]->setAvailable(true);
          }
          finished_ribosomes++;
        }
      }
    }
    // Update time
    tau = (1.0 / a0) * log(1.0 / r1);
    if (is_logging_codon_state) {
      // add state reaction to the codon's history
      codons_vector[static_cast<std::size_t>(current_codon)]
          ->addReactionToHistory(reaction_index[static_cast<std::size_t>(
                                     selected_alpha_vector_index)],
                                 tau);
    }
    clock += tau;
    i++;  // update iteration number.
  }
}

/**
 * @brief Returns a tuple where the first element is a vector with the
 * enlogation duration of the ribosomes that terminated in the simulation, and
 * the second element is a vector with the iteration where such ribosomes
 * started enlogating. This method should be called after updateRibosomeHistory,
 * since it uses the positions_vector to do its job.
 *
 */
std::tuple<std::vector<double>, std::vector<int>>
Simulations::Translation::getEnlogationDuration() {
  if (enlongations_durations.empty()) {
    getInitiationEnlongationTermination();
  }
  return std::make_tuple(enlongations_durations, initiation_iteration);
}

void Simulations::Translation::getInitiationEnlongationTermination() {
  initiations_durations.clear();
  enlongations_durations.clear();
  terminations_durations.clear();
  initiation_iteration.clear();

  std::deque<int> indexes;  // array with the index number of the ribosomes
  indexes.clear();
  std::list<int> initiations, enlongations, terminations;
  std::size_t ribosomes_to_ignore = ribosome_positions_history[0].size();
  std::size_t last_position = codons_vector.size() - 1,
              previous_size = ribosomes_to_ignore;
  for (std::size_t i = 0; i < ribosomes_to_ignore; i++) {
    indexes.push_back(static_cast<int>(indexes.size()));
    initiations_durations.push_back(0);
    enlongations_durations.push_back(0);
    terminations_durations.push_back(0);
    initiation_iteration.push_back(0);
  }
  for (std::size_t i = 1; i < ribosome_positions_history.size(); i++) {
    std::vector<int>& rib_positions = ribosome_positions_history[i];
    for (std::size_t j = 0; j < rib_positions.size(); j++) {
      std::size_t pos = static_cast<std::size_t>(rib_positions[j]);
      if (pos == 0) {
        // initiating.
        if (rib_positions.size() > previous_size) {
          // adjust offset for the iteration;
          indexes.push_front(static_cast<int>(indexes.size()));
          // new ribosome initiating.
          initiations_durations.push_back(0);
          enlongations_durations.push_back(0);
          terminations_durations.push_back(0);
          initiation_iteration.push_back(static_cast<int>(i));
        }
        initiations_durations[static_cast<std::size_t>(indexes[j])] +=
            dt_history[i];
      } else if (pos == last_position) {
        // started terminating
        terminations_durations[static_cast<std::size_t>(indexes[j])] +=
            dt_history[i];
      } else {
        // enlongating codon.
        enlongations_durations[static_cast<std::size_t>(indexes[j])] +=
            dt_history[i];
      }
    }
    previous_size = rib_positions.size();
  }
  if (ribosomes_to_ignore > 0) {
    // remove pre-filled ribosomes.
    initiations_durations.erase(
        initiations_durations.begin(),
        initiations_durations.begin() + static_cast<int>(ribosomes_to_ignore));
    enlongations_durations.erase(
        enlongations_durations.begin(),
        enlongations_durations.begin() + static_cast<int>(ribosomes_to_ignore));
    terminations_durations.erase(
        terminations_durations.begin(),
        terminations_durations.begin() + static_cast<int>(ribosomes_to_ignore));
    initiation_iteration.erase(
        initiation_iteration.begin(),
        initiation_iteration.begin() + static_cast<int>(ribosomes_to_ignore));
  }
  // maybe some of these ribosomes did not terminated. remove them from the
  // list.
  unsigned int ribosomes_to_remove = 0;
  for (std::size_t i = terminations_durations.size() - 1; i > 0; i--) {
    if (terminations_durations[i] == 0.0) {
      ribosomes_to_remove++;
    } else {
      break;
    }
  }

  for (unsigned int i = 0; i < ribosomes_to_remove; i++) {
    initiations_durations.pop_back();
    enlongations_durations.pop_back();
    terminations_durations.pop_back();
    initiation_iteration.pop_back();
  }
}

void Simulations::Translation::calculateAverageTimes() {
  std::size_t number_codons = codons_vector.size();
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
  for (auto ribosome_vector : ribosome_positions_history) {
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

void Simulations::Translation::setLogCodonStates(bool log) {
  is_logging_codon_state = log;
}

std::vector<std::tuple<std::vector<int>, std::vector<double>>>
Simulations::Translation::getLogCodonStates() {
  std::vector<std::tuple<std::vector<int>, std::vector<double>>> result(
      codons_vector.size());
  std::vector<int> state;
  std::vector<double> dt;
  for (std::size_t i = 0; i < codons_vector.size(); i++) {
    std::tie(state, dt) = codons_vector[i]->getHistory();
    result[i] = std::tie(state, dt);
  }
  return result;
}
