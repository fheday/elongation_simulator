/*
 * @file  elongation_simulation_manager.cpp
 * 
 * @brief manager class to load and execute parallel simulations
 *
 * @author Fabio Hedayioglu
 * Contact: fheday@gmail.com
 *
 */

#include "elongation_simulation_manager.h"
#include <string>
#include <thread>
#include <numeric>
#include <iostream>
#include <fstream>
#include "jsoncpp/json/json.h"
#include "elongation_simulation_processor.h"
#include "threadPool.h"

#if defined(COMIPLE_PYTHON_MODULE) || defined(TRANSLATIONSIMULATOR)

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

void init_simulation_manager(py::module &mod) {
  mod.doc() = R"pbdoc(
    Module for simulating many mRNA translations in parallel. 
    This module simulates many elongations in either different mRNAs our the same mRNA. 
    It make use of parallel processing. In order to instantiate this class, you must
    have simulation configuration file, generated by elongation.simulationbuilder.Gui()
    It will have all the information needed to perform multiple simulations in parallel.

    If the execution is cancelled before simulating all mRNAs, wen the simulation is re-run,
    It will keep the already finalised simulations and run the ones not yet run. 
    If you want to re-run an already run simulation, just delete its result file.
    )pbdoc";
  py::class_<Elongation_manager::SimulationManager>(mod, "SimulationManager")
      .def(py::init<std::string>(), py::arg("configuration_file_path"),R"docstr(
        Creates a simulator manager for the simulation configuration passed as parameter.
        configuration_file_path: string with the path to the simulation configuration file.
      )docstr") // constructor
      .def("getConcentrationFilePath", &Elongation_manager::SimulationManager::get_concentration_file_path,R"docstr(
        Method to return the string with the simulation file path. This can be used to check the parameters of the simulation set.
      )docstr")
      .def("getPrePopulate", &Elongation_manager::SimulationManager::get_pre_populate,R"docstr(
        Method to query if the simulation will pre-populate the mRNAs before runing the simulation. 
        This can be used to check the parameters of the simulation set.
      )docstr")
      .def("getStopConditionType", [](Elongation_manager::SimulationManager& sm){
        std::string result;
        switch (sm.get_stop_condition_type())
        {
        case Elongation_manager::ITERATION:
          result = "ITERATION";
          break;
        case Elongation_manager::TIME:
          result = "TIME";
          break;
        case Elongation_manager::RIBOSOMES:
          result = "RIBOSOMES";
          break;
        case Elongation_manager::STEADY_STATE_TIME:
          result = "TIME AFTER STEADY STATE";
          break;
        case Elongation_manager::STEADY_STATE_RIBOSOMES:
          result = "TERMINATED RIBOSOMES AFTER STEADY STATE";
          break;
        }
        return result;
        },R"docstr(
          Method check the simulations stop condition. This can be used to check the parameters of the simulation set.
        )docstr")
      .def("getStopConditionValue", &Elongation_manager::SimulationManager::get_stop_condition_value,R"docstr(
        Method to return value of the stop condition, if applicable. 
        E.g.: if stop condition is, say 5 terminating ribosomes. This method will return 5.
        This can be used to check the parameters of the simulation set.
      )docstr")
      .def("getHistorySize", &Elongation_manager::SimulationManager::get_history_size,R"docstr(
        Method to return maximum number on entries in the history of each simulation. 
        This can be used to check the parameters of the simulation set.
      )docstr")
      .def("start", &Elongation_manager::SimulationManager::start, py::call_guard<py::gil_scoped_release>(),R"docstr(
        Run the simulations in parallel. The number of simulations run at the same time is the number of threads in the system.
        When a simulation is finished, a json file with its results is saved.
      )docstr")
      .def("set_save_collisions", &Elongation_manager::SimulationManager::set_save_collisions, py::arg("save_collision"), R"docstr(
        If set to True, will also save the ribosomes collisions. By definition, say there are ribosomes at codons 30, 20 and 10,
        the algorithm considers ribosomes at codons 20 and colliding (there are no ribosomes preventing the moviment of the one at codon 30,
        therefore it is not colliding).
        save_collision: Boolean to set if collisions should also be saved.
      )docstr")
      .def("set_remove_ribosome_positions", &Elongation_manager::SimulationManager::set_remove_ribosome_positions, py::arg("remove_ribosome_positions"), R"docstr(
        If set to True, do not save the history of ribosome positions.
        This could be usefull for simulations with large logs where we are only interessed in the collisions.
        remove_ribosome_positions: Boolean to set if ribosome positions should not be saved.
      )docstr");
}

#endif


Elongation_manager::SimulationManager::SimulationManager(std::string cfp) {
  configuration_file_path = cfp;
  // calculate directory of configuration file.
  const size_t last_slash_idx = cfp.rfind('/');
  if (std::string::npos != last_slash_idx)
  {
      directory = cfp.substr(0, last_slash_idx) + '/';
  }

  std::ifstream config_doc(cfp, std::ifstream::binary);
  Json::Value root; // the json document.
  config_doc >> root;
  concentration_file_path = directory +
      root.get("concentration_file", "missing").asString();
  pre_populate = root.get("pre_populate", false).asBool();

  std::string fasta_file_path;
  std::string gene_string;
  float init_rate = 0, term_rate = 0, gene_copy_number = 0;
  for (unsigned int i = 0; i < root["mRNA_entries"].size(); ++i) {
    fasta_file_path = directory + root["mRNA_entries"][i]["fasta_file"].asString();
    gene_string = root["mRNA_entries"][i]["gene"].asString();
    init_rate = root["mRNA_entries"][i]["initiation_rate"].asFloat();
    term_rate = root["mRNA_entries"][i]["termination_rate"].asFloat();
    gene_copy_number =
        root["mRNA_entries"][i]["transcript_copy_number"].asFloat();
    simulations_configurations.push_back({
      fasta_file_path, gene_string, init_rate, term_rate, gene_copy_number});
  }

  if (root.isMember("iteration_limit")) {
    stop_condition_type = ITERATION;
    stop_condition_value = root.get("iteration_limit", "-1").asFloat();
  } else if (root.isMember("time_limit")) {
    stop_condition_type = TIME;
    stop_condition_value = root.get("time_limit", "-1").asFloat();
  } else if (root.isMember("finished_ribosomes")) {
    stop_condition_type = RIBOSOMES;
    stop_condition_value = root.get("finished_ribosomes", "-1").asFloat();
  } else if (root.isMember("steady_state_time")) {
    stop_condition_type = STEADY_STATE_TIME;
    stop_condition_value = root.get("steady_state_time", "-1").asFloat();
  } else if (root.isMember("steady_state_ribosomes")) {
    stop_condition_type = STEADY_STATE_RIBOSOMES;
    stop_condition_value = root.get("steady_state_ribosomes", "-1").asFloat();
  }
  // check for propensities modifiers
  if (root.isMember("propensity_modifiers")){
    //there are modifiers.
    auto propensities_json_list = root["propensity_modifiers"].getMemberNames();

    for (auto reaction_name : propensities_json_list)
      reactions_modifiers[reaction_name] = root["propensity_modifiers"][reaction_name].asFloat();
  }

  history_size = root.get("history_size", 10000).asUInt();
  config_doc.close(); // close file
  if (!is_simulation_valid())
    std::cout << "Error in configuration\n";
};

bool file_exists(std::string file_path) {
    // check if files exist and if numerical values are valid.
    std::ifstream conf_file_path{file_path};

    if (!conf_file_path) {
      return false;
    }

    return true;
  };


bool Elongation_manager::SimulationManager::is_simulation_valid() {
  if (!file_exists(configuration_file_path)) {
    std::cout<<"Can't open configuration file.\n";
    return false;
  }
    
  if (!file_exists(concentration_file_path)) {
    std::cout<<"Can't open concentration file.\n";
    return false;
  }
  if (stop_condition_value <= 0 || history_size <= 0) {
    std::cout<< "Invalid stop condition or history size.\n";
    return false;
  }
  std::string fasta_file_path, gene_name;
  float initiation_rate, termination_rate, gene_copy_number;
  for (auto const& tup : simulations_configurations) {
    std::tie(fasta_file_path, gene_name, initiation_rate, termination_rate,
             gene_copy_number) = tup;
    if (!file_exists(fasta_file_path)) {
      std::cout<<"Can't open fasta file: "<< fasta_file_path<<"\n";
      return false;
    }
    if (gene_name.empty()) {
      std::cout<<"Empty gene name.\n";
      return false;
    }
    if (initiation_rate < 0 || termination_rate < 0 || gene_copy_number <= 0) {
      std::cout<<"Invalid initiation, termination rate or gene copy number.\n";
      return false;
    }
  }
  return true;
}

std::string
Elongation_manager::SimulationManager::get_concentration_file_path() {
  return concentration_file_path;
}

std::string
Elongation_manager::SimulationManager::get_configuration_file_path() {
  return configuration_file_path;
}

std::map<std::string, float> 
Elongation_manager::SimulationManager::get_reactions_modifiers() {
  return reactions_modifiers;
}

bool Elongation_manager::SimulationManager::get_pre_populate() {
  return pre_populate;
}

std::vector<std::tuple<std::string, std::string, float, float, float>> &
Elongation_manager::SimulationManager::get_simulations_configurations() {
  return simulations_configurations;
}

Elongation_manager::stop_condition_enum
Elongation_manager::SimulationManager::get_stop_condition_type() {
  return stop_condition_type;
}

float Elongation_manager::SimulationManager::get_stop_condition_value() {
  return stop_condition_value;
}

std::size_t Elongation_manager::SimulationManager::get_history_size() {
  return history_size;
}

void Elongation_manager::SimulationManager::set_save_collisions(bool col) {
  save_collisions = col;
}

void Elongation_manager::SimulationManager::set_remove_ribosome_positions(bool rib_pos) {
  remove_ribosome_positions = rib_pos;
}

bool Elongation_manager::SimulationManager::start(bool verbose, unsigned int n_threads) {
  auto number_of_simulations = simulations_configurations.size();
  auto do_simulation = [&](std::string concentration_file_path,
                          bool pre_populate, std::string fasta_file,
                          std::string gene, float init_rate, float term_rate,
                          stop_condition_enum stop_condition, float stop_value,
                          std::size_t log_size) {
    if (verbose) std::cout<<" Simulating file: "<<fasta_file<<", gene = "<<gene<<"\n";
    // prepare and run the simulation.
    Simulations::Translation ts;
    ts.loadConcentrations(concentration_file_path);
    ts.loadMRNA(fasta_file, gene);
    ts.setInitiationRate(init_rate);
    ts.setTerminationRate(term_rate);
    if (stop_condition == ITERATION) {
      ts.setIterationLimit(static_cast<int>(stop_value));
    } else if (stop_condition == TIME) {
      ts.setTimeLimit(stop_value);
    } else if (stop_condition == RIBOSOMES) {
      ts.setFinishedRibosomes(static_cast<int>(stop_value));
    } else if (stop_condition == STEADY_STATE_TIME) {
      ts.setSimulateToSteadyState(true);
      ts.setSteadyStateTime(stop_value);
    } else if (stop_condition == STEADY_STATE_RIBOSOMES) {
      ts.setSimulateToSteadyState(true);
      ts.setSteadyStateTerminations(static_cast<int>(stop_value));;
    }
    ts.setPrepopulate(pre_populate); // simulations pre-populate the mRNA
                                     // by default. This can be changed in
                                     // the future.
    ts.setHistorySize(log_size); // update history size
    
    if (!reactions_modifiers.empty()) {
      //modify reactions
      //get reactions

      // auto original_reactions = ts.getPropensities();
      std::vector<std::map<std::string, double>> changed_propensities_vector = ts.getPropensities();
      for (auto &codon : changed_propensities_vector){
        // std::map<std::string, double> new_propensities;
        
        for (auto item : reactions_modifiers){
          if ( codon.find(item.first) == codon.end() ) {
            continue; //not found. skip.
          } else {
            codon[item.first] *= item.second; //update reaction to desired value.
          }
        }
      }
      // int i = 0;
      // for (auto codon: changed_propensities_vector) {
      //   std::cout<<"CODON NUMBER: "<<i++<<"-----------------\n";
      //   for (auto item: codon){
      //     std::cout<<"codon[" << item.first<<"] = "<< item.second<<"\n";
      //   }
      // }

      //update reactions.
      ts.setPropensities(changed_propensities_vector);

      // auto new_propensities_vector = ts.getPropensities();
      // i = 0;
      // for (auto codon:new_propensities_vector) {
      //   for (auto item: codon) {
      //     std::cout<<"old value:"<<original_reactions[i][item.first]<<"...     ";
      //     std::cout<<"new value: "<<item.first << " = "<<item.second<< "ratio : "<<item.second/ original_reactions[i][item.first]<<"\n";
      //   }
      //   i++;
      // }

    }
    //execute simulation.
    ts.run();
    return ts;
  };
  // create thread pool of threads
  ThreadPool thread_pool(n_threads);
  std::vector<std::future<bool>> sims;
  for (std::size_t i = 0; i < number_of_simulations; i++) {
    std::string fasta_path, gene_name;
    float init_rate, term_rate, copy_number;
    std::tie(fasta_path, gene_name, init_rate, term_rate, copy_number) =
        simulations_configurations[i];
    if (file_exists(directory + gene_name + ".json")) continue; // don't simulate again.
    sims.emplace_back(thread_pool.enqueue([=]{auto result = do_simulation(concentration_file_path,
                                     pre_populate, fasta_path, gene_name,
                                     init_rate, term_rate, stop_condition_type,
                                     stop_condition_value, history_size);
                                     save_sim(result);
                                     return true;}));

  }
  for(auto && result: sims){
    auto sim = result.get();
    if (sim) {
      std::cout<<"done.\n";
    }
  }
  // int j = 1;
  // for (auto sim_item  = sims.begin(); sim_item != sims.end();) {
  //   auto sim = sim_item->get();
  //   save_sim(sim);
  //   sims.erase(sim_item);
  //   j++;
  // }

  return true;
}

bool Elongation_manager::SimulationManager::save_sim(Simulations::Translation& sim) {
  Json::Value newjson;
  newjson["fasta_file"] = sim.mrna_file_name;
  newjson["initiation_rate"] = sim.initiation_rate;
  newjson["termination_rate"] = sim.termination_rate;
  std::vector<double> clock(sim.dt_history.size());
  std::partial_sum(sim.dt_history.begin(), sim.dt_history.end(), clock.begin(), std::plus<double>());
  for (auto time:clock) newjson["clock"].append(time);
  Json::Value ribosomes_history;
  auto generate_json_vector_of_vector = [&](auto &data_vector) {
        Json::Value json_value;
        for (auto entry:data_vector){
            Json::Value entry_vector;
            for (auto element:entry) entry_vector.append(element);
            json_value.append(entry_vector);
        }
        return json_value;
    };

  newjson["elongating_ribosomes"] = generate_json_vector_of_vector(sim.ribosome_positions_history);
  
  // now post-processing.
  std::string json_file_name = directory + sim.gene_name + ".json";
  if (save_collisions || remove_ribosome_positions) {
    Simulations::SimulationProcessor processor = Simulations::SimulationProcessor(newjson, json_file_name);
    if (save_collisions) {
      processor.calculateRibosomeCollisions();
    }
    if (remove_ribosome_positions) {
      processor.removeRibosomePositions();
    }
    if (save_collisions && remove_ribosome_positions) {
      // we can save space.
      processor.packData();
    }
    processor.save(); // save simualtion file with changes.
  } else {
    // no need for post processing. save file here.
    // write in a nice readible way
    Json::StreamWriterBuilder builder;
    builder["indentation"] = "   ";   // four spaces identation
    std::unique_ptr<Json::StreamWriter> writer(builder.newStreamWriter());
    std::ofstream config_doc_writer(json_file_name,
                                    std::ifstream::binary);
    writer->write(newjson, &config_doc_writer);
    
  }

  return true;
}