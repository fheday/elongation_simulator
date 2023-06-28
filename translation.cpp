/*
 * @file  translation.cpp
 * 
 * @brief class implementing the translation process
 *
 * @author Fabio Hedayioglu
 * Contact: fheday@gmail.com
 *
 */

#include "translation.h"

#include <algorithm>
#include <deque>
#include <float.h>
#include <numeric>

#include "initiationterminationcodon.h"
#include "mrna_reader.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <list>
#include <random>


#define RIBOSOME_SIZE 10

#if defined(COMIPLE_PYTHON_MODULE) || defined(TRANSLATIONSIMULATOR)

#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
namespace py = pybind11;

void init_simulation_manager(py::module &); // declare simulation manager
void init_simulation_processor(py::module &); //declare simulation processor

PYBIND11_MODULE(translation, mod) {
    mod.doc() = R"pbdoc(
    Module for simulating mRNA translation. 
    This module simulates elongation, given a tRNA concentration, initation rate, termination rate and a stop condition.
    After a simulation has been run, logs and statistics are also available.
    )pbdoc";
  py::class_<Simulations::Translation>(mod, "translation")
      .def(py::init<>(), R"docstr(
        Creates an empty simulator.
      )docstr") // constructor
      .def("loadMRNA", (void (Simulations::Translation::*) (const std::string &)) &Simulations::Translation::loadMRNA, R"docstr(
        Reads a Fasta file. Assumes there is only one gene in the file.
        file_name: string with the FASTA file to be read. All occurrences of 'T' will be replaced by 'U'.
      )docstr")
      .def("loadMRNA", (void (Simulations::Translation::*) (const std::string &, const std::string &)) &Simulations::Translation::loadMRNA, R"docstr(
        Reads a specific gene in a Fasta file.
        file_name: string with the FASTA file to be read. All occurrences of 'T' will be replaced by 'U'.
        gene_name: name of the gene to be read.
      )docstr")
      .def("inputMRNA", &Simulations::Translation::inputMRNA, R"docstr(
        Allows to pass the mRNA to be simulated as a string.
        user_mrna: String with the gene sequence. All occurrences of 'T' will be replaced by 'U'.
      )docstr")
      .def("loadConcentrations", &Simulations::Translation::loadConcentrations, R"docstr(
        Loads a csv file containing the concentrations to be used in this simulation.
        file_name: string with the path to the file containing the concentrations.
      )docstr")
      .def("loadConcentrationsFromString", &Simulations::Translation::loadConcentrationsFromString, R"docstr(
        Loads a csv string containing the concentrations to be used in this simulation.
        data: string containing the concentrations. this could be the content of the csv concentrations file.
      )docstr")
      .def("setInitiationRate", &Simulations::Translation::setInitiationRate, R"docstr(
        Sets the initiation rate in initations/sec.
        init_rate: float with the initiation rate.
      )docstr")
      .def("setTerminationRate", &Simulations::Translation::setTerminationRate, R"docstr(
        Sets the termination rate in initations/sec.
        term_rate: float with the termination rate.
      )docstr")
      .def("setIterationLimit", &Simulations::Translation::setIterationLimit, R"docstr(
        Stop condition. Only one can be set.
        Sets number of iterations (ribosome movements) wich will halt the simulation.
        i: integer with the maximum number of iterations to simulate.
      )docstr")
      .def("setTimeLimit", &Simulations::Translation::setTimeLimit, R"docstr(
        Stop condition. Only one can be set.
        Set the time (in cell time, seconds) by wich the simulation will halt.
        t: float maximum time where the simulation will halt.
      )docstr")
      .def("setFinishedRibosomes",
           &Simulations::Translation::setFinishedRibosomes, R"docstr(
             Stop condition. Only one can be set.
             Set the maximum number of terminating ribosomes. The simulation will halt when this number of ribosomes finishes.
             n_ribosomes: integer with the maximum number of terminating ribosomes to simualte.
           )docstr")
      .def("setSimulateToSteadyState", &Simulations::Translation::setSimulateToSteadyState, R"docstr(
        Stop condition. This one has to be set with one of these: setSteadyStateTime or setSteadyStateTerminations.
        When setSimulateToSteadyState is set to True, the simulator will first have to reach a steady state situation.
        A steady state situation is when (rate of observed initiations/rate of observed terminations) is between 0.9 and 1.1
        It is worth noticing that depending on the parameters set for the simulation, this could be unreacheable.
      )docstr")
      .def("setSteadyStateTime", &Simulations::Translation::setSteadyStateTime, R"docstr(
        Once the simulation enters the steady state, it will run for at most the time set in this method.
        time: float with the time (in seconds) the simulation will run after reaching steady state.
      )docstr")
      .def("setSteadyStateTerminations", &Simulations::Translation::setSteadyStateTerminations, R"docstr(
        Once the simulation enters the steady state, it will run for at most the number of ribosomes informed in this method terminates.
        terminations: integer with the max number of ribosomes to terminate after reaching steady state.
      )docstr")
      .def("setHistorySize", &Simulations::Translation::setHistorySize, R"docstr(
        Sets the maximum size of the history log (ribosome positions and time) in entries (lines). If the simulation generates more entries, the old ones
        will be removed. Default size = 100000 entries.
        size: integer with the maximum size of the log in number of entries.
      )docstr")
      .def("run", &Simulations::Translation::run,
           py::call_guard<py::gil_scoped_release>(), R"docstr(
             Runs the simulation.
             For this to work, it needs to have set up:
             * loadConcentrations
             * loadMRNA or inputMRNA
             * setInitiationRate
             * setTerminationRate
             * setFinishedRibosomes or setIterationLimit or setTimeLimit 
               or setSimulateToSteadyState and either setSteadyStateTerminations or setSteadystateTime
           )docstr")
      .def("getElongationDuration",
           &Simulations::Translation::getElongationDuration, R"docstr(
             Analysis the simulation log (dt_history and ribosome_positions_history) and calculates the times (in seconds) each ribosome takes from initiation to termination.
             returns two lists: 
             the first list contains the durations
             the second contains the indexes in dt_history and ribosome_positions_history when the ribosome initiated
           )docstr")
      .def("setPrepopulate", &Simulations::Translation::setPrepopulate, R"docstr(
        This is an optional method to be called BEFORE the simulation.
        It tries to create an approximated configuration where ribosomes would be located in the mRNA
        in a steady state, and use this configutation as the starting state for the simulation. 
      )docstr")
      .def("getRibosomesPositions",
           &Simulations::Translation::getRibosomesPositions, R"docstr(
             Return the codon number of all ribosomes in the current simulation state.
           )docstr")
      .def("setRibosomePositions",
           &Simulations::Translation::setRibosomePositions, R"docstr(
             Set ribosome positions in the mRNA strip. Used before starting the
             simulation.
             positions: list of integers with the positions of the ribosomes in the mRNA.
           )docstr")
      .def("setLogCodonStates", &Simulations::Translation::setLogCodonStates, R"docstr(
        If set to true, in addition to logging the ribosome positions in the mRNA, it also log the ribosome's internal states.
      )docstr")
      .def("getLogCodonStates", &Simulations::Translation::getLogCodonStates, R"docstr(
        This method returns a list of lists with the states of the ribosomes in each codon as follows:
        Each element of the list represents a codon in the mRNA.
        For each codon there will be two lists: 
        - the first one with the state, 
        - the second one, with the total time spent in that state.
      )docstr")
      .def("setPropensities", &Simulations::Translation::setPropensities, R"docstr(
        This method changes the reactions propensities of all codons.
        prop: a vector with the same size as the mRNA where each entry consists of a dictionary 
        with new propensities. The original propensities vector can obtained by calling getPropensities().
      )docstr")
      .def("setNonCognate", &Simulations::Translation::setNoNonCognate, R"docstr(
        Optimization: this option will disable the non-cognate pathway in the ribosomes.
      )docstr")
      .def("getPropensities", &Simulations::Translation::getPropensities, R"docstr(
        This method returns a vector with the same size as the mRNA where each entry consists of 
        a dictionary with the reactions labels and their propensities.
        The vector of dictionaries returned by this method can be changed and used as an input parameter for 
        setPropensities, in order to change a specific reaction's propensity.
      )docstr")

      .def_readonly("mrna_file_name", &Simulations::Translation::mrna_file_name, R"docstr(
        Atribute: string with the path to mRNA file, set with loadMRNA.
      )docstr")
      .def_readonly("concentrations_file_name",
                    &Simulations::Translation::concentrations_file_name, R"docstr(
                      Atribute: string with the path to the concentrations file, set with loadConcentrations.
                    )docstr")
      .def_readonly("dt_history", &Simulations::Translation::dt_history, R"docstr(
        Attribute: The time taken by each reaction. This numpy array is filled after a simulation has been run.
      )docstr")
      .def_readonly("ribosome_positions_history",
                    &Simulations::Translation::ribosome_positions_history, R"docstr(
                      Attribute: List of lists where each entry is created when a ribosome moves, inititates or terminates elongation.
                      each entry constitutes of the positions of all ribosomes in the mRNA. 
                      All codons are zero-based. Each entry in this list has a corresponding dt_history with the time taken between 
                      two consecutive entries of this list.
                    )docstr")
      .def_readonly("initiationRate",
                    &Simulations::Translation::initiation_rate, R"docstr(
                      Attribute: the initiation rate set by setInitiationRate.
                    )docstr")
      .def_readonly("terminationRate",
                    &Simulations::Translation::termination_rate, R"docstr(
                      Attribute: the termination rate set by setInitiationRate.
                    )docstr")
      .def_property_readonly("elongations_durations",[](Simulations::Translation &sim){
                      return sim.getElongationDuration();
                    }, R"docstr(
                      Attribute: Durations of each completed elongation. The first list contains the duration in seconds of the time
                      taken by a ribosome from initiation to termination, the second list contains the entry number where the ribosome
                      initiated. This attribute is only populated AFTER getElongationDuration method is run.
                    )docstr")
      .def_readonly("total_time", &Simulations::Translation::total_time, R"docstr(
        Attribute: total time ribosomes spent in each codon. Populated after calling getAverageTimes method.
      )docstr")
      .def_readonly("n_times_occupied",
                    &Simulations::Translation::n_times_occupied, R"docstr(
                      Attribute: vector with number of times each codon is occupied. Populated after calling getAverageTimes method.
                    )docstr")
      .def_property_readonly("average_times", [](Simulations::Translation &sim) {
             if (sim.codons_average_occupation_time.empty() && !sim.ribosome_positions_history.empty()) sim.getAverageTimes();
             return sim.codons_average_occupation_time;
           }
          , R"docstr(
            Attribute: vector with average time each codon is occupied by a ribosome. Populated after calling getAverageTimes method.
          )docstr")
      .def_property_readonly("colliding_ribosomes", [](Simulations::Translation &sim) {
             sim.getRibosomeCollisions();
             return sim.colliding_ribosomes;
           }
          , R"docstr(
            Attribute: vector with the positions of the colliding ribosomes. By definition a ribosome is only colliding if another one is stopping it moving forward.
          )docstr")
      .def_property_readonly("stalled_ribosomes", [](Simulations::Translation &sim) {
             sim.getRibosomeCollisions();
             return sim.stalled_ribosomes;
           }
          , R"docstr(
            Attribute: vector with the positions of stalling ribosomes. By definition a stalled ribosome can move forward but is blocking another ribosome to move forward.
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
                             )docstr");;

    init_simulation_manager(mod); //include simulation manager to package
    init_simulation_processor(mod); //include simulation processor to package
}

#endif

void Simulations::Translation::loadConcentrations(
    const std::string &file_name) {
  std::ifstream ist{file_name};

  if (!ist) {
    throw std::runtime_error("can't open input file: " + file_name);
  } else {
    concentrations_file_name = file_name;
    concentrations_source = "file_name";
    initializeMRNAReader();
  }
}

void Simulations::Translation::loadConcentrationsFromString(const std::string & data) {
  concentrations_string = data;
  concentrations_source = "string";
	initializeMRNAReader();
}


void Simulations::Translation::loadMRNA(const std::string &file_name) {
  std::ifstream ist{file_name};

  if (!ist) {
    throw std::runtime_error("can't open input file: " + file_name);
  } else {
    mrna_file_name = file_name;
    initializeMRNAReader();
  }
}

void Simulations::Translation::loadMRNA(const std::string &file_name, const std::string &g_n) {
  std::ifstream ist{file_name};

  if (!ist) {
    throw std::runtime_error("can't open input file: " + file_name);
  } else {
    mrna_file_name = file_name;
    gene_name = g_n;
    initializeMRNAReader();
  }
}

void Simulations::Translation::inputMRNA(std::string user_mRNA) {
  mrna_input = user_mRNA;
  initializeMRNAReader();
}


void Simulations::Translation::initializeMRNAReader() {
  if (concentrations_source != "None" && (!mrna_file_name.empty() || !mrna_input.empty() ) &&
      is_initiation_set && is_termination_set) {
    // we have the concentrations and mrna file names. we can proceed.
    mRNA_utils::mRNAReader mrr;
    if (gene_name.empty() && !mrna_file_name.empty()){ 
      mrr.loadmRNAFile(mrna_file_name);
    } else if (!gene_name.empty() && !mrna_file_name.empty()){
      mrr.loadGene(mrna_file_name, gene_name);
    } else {
      mrr.readMRNA(mrna_input);
    }
    std::vector<std::string> stop_codons = {"UAG", "UAA",
                                          "UGA"};  // list of stop codons.
    // fill codon vector.
    int n_codons = mrr.sizeInCodons();
    std::unique_ptr<Simulations::InitiationTerminationCodon> initiation_codon(
        new Simulations::InitiationTerminationCodon(initiation_rate, true));
    initiation_codon->codon = mrr.getCodon(0);
    initiation_codon->index = 0;
    initiation_codon->setState(0);
    codons_vector.push_back(std::move(initiation_codon));
    for (int i = 1; i < n_codons - 1; i++) {
      //check if codon is a stop codon.
      if(std::find(stop_codons.begin(), stop_codons.end(), mrr.getCodon(i)) != stop_codons.end()) {
        //stop codon found. we shouldn't have it here.
        std::cout<<"Stop codon found before the end of gene. Not simulating.\n";
        is_mRNA_valid =false;
        break;
      }
      std::unique_ptr<Simulations::ElongationCodon> c(
          new Simulations::ElongationCodon());
      if (concentrations_source == "file_name"){
        c->loadConcentrations(concentrations_file_name);
      } else if (concentrations_source == "string") {
        c->loadConcentrationsFromString(concentrations_string);
      }
      
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

void Simulations::Translation::setPropensities(std::vector<std::map<std::string, double>> prop) {
  changed_propensities = true;
  for (std::size_t i = 1; i < codons_vector.size() - 1; i++) {
    codons_vector[i]->setPropensities(prop[i]);
  }
}

void Simulations::Translation::setNoNonCognate(bool noNonCog) {
  no_noCognate = noNonCog;
  for (std::size_t i = 1; i < codons_vector.size() - 1; i++) {
    codons_vector[i]->setNoNonCognate(noNonCog);
  }
}

std::vector<std::map<std::string, double>>
Simulations::Translation::getPropensities() {
  auto result = std::vector<std::map<std::string, double>>();
  result.push_back(std::map<std::string, double>()); // codon 0 will be empty.
  for (std::size_t i = 1; i < codons_vector.size() - 1; i++) {
    result.push_back(codons_vector[i]->getPropensities());
  }
  result.push_back(std::map<std::string, double>()); // last codon will also be empty.

  return result;
}

void Simulations::Translation::setInitiationRate(double ir) {
  if (ir >= 0) {
    initiation_rate = ir;
  }
  is_initiation_set = true;
  initializeMRNAReader();
}

void Simulations::Translation::setTerminationRate(double tr) {
  if (tr >= 0) {
    termination_rate = tr;
  }
  is_termination_set = true;
  initializeMRNAReader();
}

void Simulations::Translation::setPrepopulate(bool prep) {
  pre_populate = prep;
}

void Simulations::Translation::setHistorySize(std::size_t size) {
  history_size = size;
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

void Simulations::Translation::setSimulateToSteadyState(bool ss) {
  simulate_to_steady_state = ss;
};

void Simulations::Translation::setSteadyStateTime(float t) {
  // only add this condition if terminations is not set.
  if (t >= 0 && steady_state_terminations < 0) steady_state_time = t;
}

void Simulations::Translation::setSteadyStateTerminations(int t) {
  // only add this condition if time is not set.
  if (t >= 0 && steady_state_time < 0) steady_state_terminations = t;
}

void Simulations::Translation::getAlphas(utils::circular_buffer<std::vector<int>>& ribosome_positions_history_circ_buffer) {
  std::size_t global_index = 0;

  // add initiation if needed.
  if (initiation_rate && codons_vector[0]->isAvailable()) {
    // need to add initalization.
    for (global_index = 0; global_index < codons_vector[0]->alphas.size();
         global_index++) {
      alphas[global_index] = codons_vector[0]->alphas[global_index];
      codon_index[global_index] = 0;
      reaction_index[global_index] = codons_vector[0]->alphas[global_index];
    }
  }
  for (auto ribosome_index : ribosome_positions_history_circ_buffer.peek_back()) {
    for (std::size_t index = 0;
         index < codons_vector[ribosome_index]->alphas.size(); index++) {
      alphas[global_index] = codons_vector[ribosome_index]->alphas[index];
      codon_index[global_index] = ribosome_index;
      reaction_index[global_index] =
          codons_vector[ribosome_index]->reactions_index[index];
      global_index++;
    }
  }
  global_size = global_index; // update global size.
}

void Simulations::Translation::insertRibosome(std::size_t position,
                                              bool set_neighborhood = false) {
  codons_vector[position]->setOccupied(true);
  codons_vector[position]->setAvailable(false);
  codons_vector[position]->setState(0);
  if (position == 0) {
    codons_vector[position]->setState(23);
  }
  if (set_neighborhood) {
    for (int i = 0;
         i < RIBOSOME_SIZE && static_cast<std::size_t>(i) <= position; i++) {
      codons_vector[position - static_cast<std::size_t>(i)]->setAvailable(
          false);
    }
  }
}

void Simulations::Translation::run() {
  utils::circular_buffer<double> dt_history_circ_buffer(history_size);
  utils::circular_buffer<std::vector<int>> ribosome_positions_history_circ_buffer(history_size);
  // initialize the random generator
  std::random_device
      rd; // Will be used to obtain a seed for the random number engine
  std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
  std::uniform_real_distribution<> dis(DBL_MIN, 1);

  double r1 = 0, r2 = 0;
  double tau = 0, clock = 0.0;
  int i = 0;

  int finished_ribosomes = 0, pre_filled_ribosomes = 0;
  // pre-allocate space for some vectors.
  std::size_t max_ribosomes = (codons_vector.size() / RIBOSOME_SIZE) + 1;
  utils::circular_buffer<int> rib_positions(max_ribosomes);
  alphas.resize(4 * max_ribosomes);
  codon_index.resize(4 * max_ribosomes);
  reaction_index.resize(4 * max_ribosomes);

  //check if the mRNA is greater than 2 ribosomes. if not, print error. do not simulate.
  if (codons_vector.size() < 2 * RIBOSOME_SIZE) {
    std::cout<<"mRNA too short. not simulating.\n";
    return;
  }
  //check if mRNA is valid
  if (is_mRNA_valid == false) {
    std::cout<<"Can't elongate invalid mRNA. Quitting.\n";
    return;
  }

  // pre-fill codons based on the rates.
  // do not pre-fill if ribosome propensities have been changed.
  // do not pre-fill if initiation rate is zero.
  if (pre_populate && !changed_propensities && initiation_rate > 0) {
    std::size_t last_index = 0;
    double time_sum = 0;
    std::map<std::string, double> estimated_codon_time{
        {"AAA", 0.033433631062507}, {"AAC", 0.076520502567291}, 
        {"AAG", 0.213781446218490}, {"AAU", 0.795929551124572},
        {"ACA", 0.507978320121765}, {"ACC", 0.791741490364074},  
        {"ACG", 0.631365478038787}, {"ACU", 0.231777504086494},
        {"AGA", 0.032752707600593}, {"AGC", 0.159446164965629},
        {"AGG", 0.687197327613830}, {"AGU", 0.927559494972229},
        {"AUA", 0.629360139369964}, {"AUC", 0.721506834030151},
        {"AUG", 0.223424538969994}, {"AUU", 0.112727925181389},
        {"CAA", 0.086375549435616}, {"CAC", 0.176103606820106},
        {"CAG", 0.820655703544617}, {"CAU", 0.953253030776978},
        {"CCA", 0.274570554494858}, {"CCC", 2.00137615203857},
        {"CCG", 0.814045131206512}, {"CCU", 0.841190457344055},
        {"CGA", 0.892297625541687}, {"CGC", 0.809189140796661},
        {"CGG", 1.16730070114136},  {"CGU", 0.0852445140481},
        {"CUA", 0.473841965198517}, {"CUC", 1.00331914424896},
        {"CUG", 1.51964175701141},  {"CUU", 1.32332634925842},
        {"GAA", 0.064450852572918}, {"GAC", 0.103697247803211},
        {"GAG", 0.717609882354736}, {"GAU", 0.839430928230285},
        {"GCA", 0.544444143772125}, {"GCC", 0.876713514328003},
        {"GCG", 0.904575228691101}, {"GCU", 0.316948801279068},
        {"GGA", 0.276408851146698}, {"GGC", 0.054870706051588},
        {"GGG", 0.662749826908112}, {"GGU", 0.762714862823486},
        {"GUA", 0.659511029720306}, {"GUC", 0.897049844264984},
        {"GUG", 0.753913283348084}, {"GUU", 0.139686241745949},
        {"UAC", 0.175240993499756}, {"UAU", 0.953546106815338},
        {"UCA", 0.623058080673218}, {"UCC", 0.881772577762603},
        {"UCG", 0.836475849151611}, {"UCU", 0.131028860807419},
        {"UGC", 0.170239016413689}, {"UGG", 0.174226880073547},
        {"UGU", 0.944100022315979}, {"UUA", 0.187631607055664},
        {"UUC", 0.131067156791687}, {"UUG", 0.287777960300446},
        {"UUU", 0.882003247737884}
    };
    double initiation_time = 0;
    for (std::size_t i = 0; i < 10; ++i) initiation_time += estimated_codon_time[codons_vector[i]->codon];
    
    insertRibosome(last_index, true);
    for (std::size_t i = RIBOSOME_SIZE; i < codons_vector.size(); ++i) {
      if (i - last_index < RIBOSOME_SIZE)
        continue;

      time_sum += estimated_codon_time[codons_vector[i]->codon];
      if (time_sum >= initiation_time) {
        // put a ribosome here.
        insertRibosome(i, true);
        time_sum = 0;   // reset timer.
        last_index = i; // mark this as last inserted ribosome.
      }
    }
  }
  for (int i = static_cast<int>(codons_vector.size()) - 1; i >= 0; i--) {
    if (codons_vector[i]->isOccupied()) {
      rib_positions.put(static_cast<int>(i));
      pre_filled_ribosomes++;
    }
  }
  ribosome_positions_history_circ_buffer.put(rib_positions.get_vector(true));
  dt_history_circ_buffer.put(0.0);
  finished_ribosomes -=
      pre_filled_ribosomes; // we should ignore these ribosomes.
      
  std::size_t moved_codon = 0, current_codon = 0;
  bool initiation = false, termination = false, moved = true;

  double cumsum = 0, a0 = 0;
  std::size_t selected_alpha_vector_index = 0;
  bool steady_state = false;
  bool steady_state_stop_condition_met = false;
  float initiations = 0, terminations = 0;
  int n_total_initiations = 0, n_total_terminations = 0;
  int steady_state_performed_terminations = 0;
  float steady_state_elapsed_time = 0;
  float last_initiation = -1, last_termination = -1;
  while ((iteration_limit > 0 && i < iteration_limit) ||
         (time_limit > 0 && clock < time_limit) ||
         (finished_ribosomes_limit > 0 &&
          finished_ribosomes_limit > finished_ribosomes) ||
          (simulate_to_steady_state && !steady_state_stop_condition_met)) {
    moved = false;
    initiation = false;
    termination = false;
    // randomly generate parameter for calculating dt
    r1 = dis(gen);
    // randomly generate parameter for selecting reaction
    r2 = dis(gen);
    // calculate an
    getAlphas(ribosome_positions_history_circ_buffer);
    if (global_size == 0) {
      // no available reactions, quit loop prematurely.
      std::cout << "no available reactions. quitting.\n";
      break;
    }
    a0 = std::accumulate(alphas.begin(), alphas.begin() + global_size, 0.0);
    selected_alpha_vector_index = 0;
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
    cumsum = alphas[selected_alpha_vector_index];
    while (cumsum < a0 * r2) {
      selected_alpha_vector_index++;
      cumsum += alphas[selected_alpha_vector_index];
    };
    current_codon = codon_index[selected_alpha_vector_index];
    // Apply reaction
    codons_vector[current_codon]->executeReaction(
        static_cast<int>(reaction_index[selected_alpha_vector_index]));
    // Update time
    tau = (1.0 / a0) * log(1.0 / r1);

    if (current_codon == 0 && codons_vector[0]->getState() == 23) {
      // initiated.
      codons_vector[0]->setAvailable(false);
      codons_vector[0]->setOccupied(true);
      initiation = true;
      moved = true;
    }
    // 2- Any codon with state == 31 means the ribosome already moved to the
    // next codon (or left the mRNA). update states.
    if (codons_vector[current_codon]->getState() == 31) {
      codons_vector[current_codon]->setState(0);
      codons_vector[current_codon]->setAvailable(false);
      codons_vector[current_codon]->setOccupied(false);
      moved = true;
      moved_codon = current_codon + 1;
      if (moved_codon < codons_vector.size()) {
        codons_vector[moved_codon]->setOccupied(true);
        codons_vector[moved_codon]->setAvailable(false);
      }
      // update free codons due to the size of the ribosome.
      // we need to do some tidying up after the ribosome.
      if ((moved_codon > RIBOSOME_SIZE - 1) &&
          (moved_codon < codons_vector.size())) {
        // update freed space left by the ribosome's movement.
        codons_vector[moved_codon - RIBOSOME_SIZE]->setAvailable(true);
      } else if (moved_codon == codons_vector.size()) {
        // ribosome terminated. free codons positions occupied by it.
        termination = true;
        for (std::size_t i = codons_vector.size() - RIBOSOME_SIZE;
             i < codons_vector.size(); ++i) {
          codons_vector[i]->setAvailable(true);
        }
        finished_ribosomes++;
      }
    }

    // update ribosome position.
    // check if there was movement.
    if (moved) {
      if (termination) {
        // terminated. remove last position.
        rib_positions.get();
      } else if (initiation) {
        // initiated.
        rib_positions.put(0);
      } else {
        rib_positions.replace(static_cast<int>(moved_codon) - 1, static_cast<int>(moved_codon));
      }
      // ribosome movement detected. create new entry in the history.
      dt_history_circ_buffer.put(tau);
      ribosome_positions_history_circ_buffer.put(rib_positions.get_vector(true));
    } else {
      // no ribosome movement. just update dt_history.
      dt_history_circ_buffer.peek_back() += tau;
    }
    if (is_logging_codon_state) {
      // add state reaction to the codon's history
      codons_vector[current_codon]->addReactionToHistory(
          reaction_index[selected_alpha_vector_index], tau);
    }
    clock += tau;
    
    if ( simulate_to_steady_state ) {
      if (steady_state) steady_state_elapsed_time += tau; // if steady state enabled, add to stady state time.
      if ( initiation ) {
        if ( last_initiation > 0 ) {
          initiations +=1 / (clock - last_initiation);
          n_total_initiations++;
        }
        last_initiation = clock;
      } else if ( termination ) {
        if ( last_termination > 0 ) {
          if (!steady_state) {
            terminations += 1/(clock - last_termination); // update termination rate
          } else {
            //add terminations in steady state.
            steady_state_performed_terminations++;
          }
          n_total_terminations++;
        }
        last_termination = clock;
      }
      //check if termination condition is met.
      if (n_total_terminations > 2 && n_total_initiations >2) { // start checking from the 3rd initiation and termination onwards.
        if (!steady_state) {
          float rate = (initiations/n_total_initiations)/(terminations/n_total_terminations);
          if (rate >=0.9 && rate <=1.1) { //average initiation rate =+- 90% termination rate.
            steady_state = true; // steady state achieved.
            if (steady_state_terminations < 0 && steady_state_time < 0) {
              // there is no additional stop conditions. done.
              steady_state_stop_condition_met = true;
            }
          }
        } else {
          // steady state already achieved. check next condition
          if (steady_state_terminations == steady_state_performed_terminations || (steady_state_time > 0 && 
              (steady_state_elapsed_time > steady_state_time))) {
            steady_state_stop_condition_met = true;
          }
        }
      }
    }
    i++; // update iteration number.
  }
  // copy log data to the object-wide log system.
  dt_history = dt_history_circ_buffer.get_vector(false); 
  ribosome_positions_history = ribosome_positions_history_circ_buffer.get_vector(false);
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
Simulations::Translation::getElongationDuration() {
  if (elongations_durations.empty() && !ribosome_positions_history.empty()) {
    getInitiationElongationTermination();
  }
  return std::make_tuple(elongations_durations, initiation_iteration);
}

void Simulations::Translation::getInitiationElongationTermination() {
  initiations_durations.clear();
  elongations_durations.clear();
  terminations_durations.clear();
  initiation_iteration.clear();

  std::deque<int> indexes; // array with the index number of the ribosomes
  indexes.clear();
  std::list<int> initiations, elongations, terminations;
  std::size_t ribosomes_to_ignore = ribosome_positions_history[0].size();
  std::size_t last_position = codons_vector.size() - 1,
              previous_size = ribosomes_to_ignore;
  for (std::size_t i = 0; i < ribosomes_to_ignore; i++) {
    indexes.push_back(static_cast<int>(indexes.size()));
    initiations_durations.push_back(0);
    elongations_durations.push_back(0);
    terminations_durations.push_back(0);
    initiation_iteration.push_back(0);
  }
  for (std::size_t i = 1; i < ribosome_positions_history.size(); i++) {
    std::vector<int> &rib_positions = ribosome_positions_history[i];
    for (std::size_t j = 0; j < rib_positions.size(); j++) {
      std::size_t pos = static_cast<std::size_t>(rib_positions[j]);
      if (pos == 0) {
        // initiating.
        if (rib_positions.size() > previous_size) {
          // adjust offset for the iteration;
          indexes.push_front(static_cast<int>(indexes.size()));
          // new ribosome initiating.
          initiations_durations.push_back(0);
          elongations_durations.push_back(0);
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
        // elongating codon.
        elongations_durations[static_cast<std::size_t>(indexes[j])] +=
            dt_history[i];
      }
    }
    previous_size = rib_positions.size();
  }
  if (ribosomes_to_ignore > 0) {
    // remove pre-filled ribosomes.
    initiations_durations.erase(initiations_durations.begin(),
                                initiations_durations.begin() +
                                    static_cast<int>(ribosomes_to_ignore));
    elongations_durations.erase(elongations_durations.begin(),
                                elongations_durations.begin() +
                                    static_cast<int>(ribosomes_to_ignore));
    terminations_durations.erase(terminations_durations.begin(),
                                 terminations_durations.begin() +
                                     static_cast<int>(ribosomes_to_ignore));
    initiation_iteration.erase(initiation_iteration.begin(),
                               initiation_iteration.begin() +
                                   static_cast<int>(ribosomes_to_ignore));
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
    elongations_durations.pop_back();
    terminations_durations.pop_back();
    initiation_iteration.pop_back();
  }
}

/*
 * @brief Return the codon number of all ribosomes in the current simualtion
 * state.
 */
std::vector<int> Simulations::Translation::getRibosomesPositions() {
  std::vector<int> result;
  for (std::size_t i = 0; i < codons_vector.size(); i++) {
    if (codons_vector[i]->isOccupied()) {
      result.emplace_back(i);
    }
  }
  return result;
}

/*
 * @brief set ribosome positions in the mRNA strip. Used before starting the
 * simulation.
 */
void Simulations::Translation::setRibosomePositions(
    std::vector<int> positions) {
  // validate: did the user passed ribosomes?
  if (positions.size() == 0) {
    throw std::out_of_range("No ribosomes in the vector...");
  }
  // validate: check if all ribosomes are inside mRNA
  if (static_cast<std::size_t>(*std::max_element(
          positions.begin(), positions.end())) >= codons_vector.size()) {
    throw std::out_of_range("Ribosome positioned after the end of mRNA.");
  }
  if (*min_element(positions.begin(), positions.end()) < 0) {
    throw std::out_of_range("Invalid (negative) position informed.");
  }
  std::sort(positions.begin(), positions.end()); // sort positions.
  // validate: minimum distance between ribosomes = RIBOSOME_SIZE.
  insertRibosome(positions[0], true);
  for (std::size_t i = 1; i < positions.size(); i++) {
    if (positions[i] - positions[i - 1] < RIBOSOME_SIZE) {
      throw std::out_of_range("Ribosome " + std::to_string(positions[i]) +
                              " too close to ribosome " +
                              std::to_string(positions[i]));
    } else {
      insertRibosome(positions[i], true);
    }
  }
}

void Simulations::Translation::getAverageTimes() {
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
  for (std::size_t i = 0; i < codons_vector.size(); i++) result[i] = codons_vector[i]->getHistory();
  return result;
}

void Simulations::Translation::getRibosomeCollisions()
{
  if (!is_collisions_calculated && !ribosome_positions_history.empty())
  {
    // calculate collisions
    for (auto ribosomes_positions : ribosome_positions_history)
    {
      std::vector<int> collision_entry;
      std::vector<int> stall_entry;
      if (ribosomes_positions.size() == 0) continue;
      for (std::size_t i = 0; i < ribosomes_positions.size() - 1; i++)
      {
        if (ribosomes_positions[i + 1] - ribosomes_positions[i] == 10)
        {
          // colliding ribosome: collision with the next ribosome detected.
          collision_entry.emplace_back(ribosomes_positions[i]);
        }
        else if (!collision_entry.empty() && ribosomes_positions[i] - collision_entry.back() == 10)
        {
          // stalled ribosome: no collision with next ribosome,
          // but collision with previous ribosome detected.
          stall_entry.emplace_back(ribosomes_positions[i]);
        }
      }
      // check last entry. it can only stall.
      if (!collision_entry.empty() && ribosomes_positions.back() - collision_entry.back() == 10)
      {
        // stalled ribosome: no collision with next ribosome,
        // but collision with previous ribosome detected.
        stall_entry.emplace_back(ribosomes_positions.back());
      }
      colliding_ribosomes.emplace_back(collision_entry);
      stalled_ribosomes.emplace_back(stall_entry);
    }
    is_collisions_calculated = true;
  }
}
