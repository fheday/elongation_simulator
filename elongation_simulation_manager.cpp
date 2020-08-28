#include "elongation_simulation_manager.h"
#include <thread>
#include <iostream>
#include <fstream>
#include "json/json.h"

#if defined(COMIPLE_PYTHON_MODULE) || defined(TRANSLATIONSIMULATOR)

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

void init_simulation_manager(py::module &mod) {
  py::class_<Elongation_manager::SimulationManager>(mod, "SimulationManager")
      .def(py::init<std::string>()) // constructor
      .def("getConcentrationFilePath", &Elongation_manager::SimulationManager::get_concentration_file_path)
      .def("getConfigurationFilePath", &Elongation_manager::SimulationManager::get_configuration_file_path)
      .def("getPrePopulate", &Elongation_manager::SimulationManager::get_pre_populate)
      // .def("getSimulationsConfigurations", &Elongation_manager::SimulationManager::get_simulations_configurations)
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
        }
        return result;
        })
      .def("getStopConditionValue", &Elongation_manager::SimulationManager::get_stop_condition_value)
      .def("getHistorySize", &Elongation_manager::SimulationManager::get_history_size)
      .def("start", &Elongation_manager::SimulationManager::start, py::call_guard<py::gil_scoped_release>());

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
  if (!file_exists(configuration_file_path))
    return false;
  if (!file_exists(concentration_file_path))
    return false;
  if (stop_condition_value <= 0 || history_size <= 0)
    return false;
  std::string fasta_file_path, gene_name;
  float initiation_rate, termination_rate, gene_copy_number;
  for (auto const& tup : simulations_configurations) {
    std::tie(fasta_file_path, gene_name, initiation_rate, termination_rate,
             gene_copy_number) = tup;
    if (!file_exists(fasta_file_path))
      return false;
    if (gene_name.empty())
      return false;
    if (initiation_rate < 0 || termination_rate < 0 || gene_copy_number <= 0)
      return false;
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

bool Elongation_manager::SimulationManager::start() {
  auto number_of_simulations = simulations_configurations.size();
  auto do_simulation = [&](std::string concentration_file_path,
                          bool pre_populate, std::string fasta_file,
                          std::string gene, float init_rate, float term_rate,
                          stop_condition_enum stop_condition, float stop_value,
                          std::size_t log_size) {
    // prepare and run the simulation.
    Simulations::Translation ts;
    ts.loadConcentrations(concentration_file_path);
    ts.loadMRNA(fasta_file, gene);
    ts.setInitiationRate(init_rate);
    ts.setTerminationRate(term_rate);
    if (stop_condition == ITERATION) {
      ts.setIterationLimit(stop_value);
    } else if (stop_condition == TIME) {
      ts.setTimeLimit(stop_value);
    } else if (stop_condition == RIBOSOMES) {
      ts.setFinishedRibosomes(stop_value);
    }
    ts.setPrepopulate(pre_populate); // simulations pre-populate the mRNA
                                     // by default. This can be changed in
                                     // the future.
    ts.setHistorySize(log_size);
    ts.run();
    return ts;
  };
  // create thread pool of threads
  auto n_threads = std::thread::hardware_concurrency();
  std::vector<std::future<Simulations::Translation>> sims;
  bool has_finished_tasks = false;
  std::size_t finished_index = 0;
  for (std::size_t i = 0; i < number_of_simulations; i++) {
    std::string fasta_path, gene_name;
    float init_rate, term_rate, copy_number;
    std::tie(fasta_path, gene_name, init_rate, term_rate, copy_number) =
        simulations_configurations[i];
    if (file_exists(directory + gene_name + ".json")) continue; // don't simulate again.
    sims.push_back(std::async(do_simulation, concentration_file_path,
                                     pre_populate, fasta_path, gene_name,
                                     init_rate, term_rate, stop_condition_type,
                                     stop_condition_value, history_size));

    if (sims.size() % n_threads == 0){
      do {
        has_finished_tasks = false;
        for (std::size_t sim_index = 0; sim_index < sims.size(); sim_index++) {
          if (sims[sim_index].wait_for(std::chrono::seconds(0)) == std::future_status::ready) {
            has_finished_tasks = true;
            finished_index = sim_index;
            break;
          }
        }
      } while(has_finished_tasks == false);
      if (has_finished_tasks) {
          sims[finished_index].wait();
          auto sim = sims[finished_index].get();
          save_sim(sim);
          sims.erase(sims.begin() + finished_index);
      }
    }
  }
  int j = 1;
  for (auto sim_item  = sims.begin(); sim_item != sims.end();) {
    auto sim = sim_item->get();
    save_sim(sim);
    sims.erase(sim_item);
    j++;
  }

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
  for (auto entry:sim.ribosome_positions_history){
    Json::Value entry_vector;
    for (auto ribosome:entry) entry_vector.append(ribosome);
    ribosomes_history.append(entry_vector);
  }
  newjson["elongating_ribosomes"] = ribosomes_history;
  
  std::ifstream file_stream(configuration_file_path, std::ifstream::binary);
  // write in a nice readible way
  Json::StreamWriterBuilder builder;
  // builder["commentStyle"] = "None"; // no comments.
  builder["indentation"] = "   ";   // four spaces identation
  std::unique_ptr<Json::StreamWriter> writer(builder.newStreamWriter());
  std::ofstream config_doc_writer(directory + sim.gene_name+".json",
                                  std::ifstream::binary);
  writer->write(newjson, &config_doc_writer);
  return true;
}