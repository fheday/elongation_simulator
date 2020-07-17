#include "elongation_simulation_manager.h"
#include "jsoncpp.cpp"
#include <fstream>
#include <thread>

Elongation_manager::SimulationManager::SimulationManager(std::string cfp) {
  configuration_file_path = cfp;
  std::ifstream config_doc(cfp, std::ifstream::binary);
  Json::Value root; // the json document.
  config_doc >> root;
  concentration_file_path =
      root.get("concentration_file", "missing").asString();
  pre_populate = root.get("pre_populate", false).asBool();

  std::string fasta_file_path;
  std::string gene_string;
  float init_rate = 0, term_rate = 0, gene_copy_number = 0;
  for (unsigned int i = 0; i < root["mRNA_entries"].size(); ++i) {
    fasta_file_path = root["mRNA_entries"][i]["fasta_file"].asString();
    gene_string = root["mRNA_entries"][i]["gene"].asString();
    init_rate = root["mRNA_entries"][i]["initiation_rate"].asFloat();
    term_rate = root["mRNA_entries"][i]["termination_rate"].asFloat();
    gene_copy_number =
        root["mRNA_entries"][i]["transcript_copy_number"].asFloat();
    simulations_configurations.push_back(std::tuple(
        fasta_file_path, gene_string, init_rate, term_rate, gene_copy_number));
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

bool Elongation_manager::SimulationManager::is_simulation_valid() {
  auto file_exists = [](std::string file_path) {
    // check if files exist and if numerical values are valid.
    std::ifstream conf_file_path{file_path};

    if (!conf_file_path) {
      throw std::runtime_error("can't open input file: " + file_path);
      return false;
    }

    return true;
  };
  if (!file_exists(configuration_file_path))
    return false;
  if (!file_exists(concentration_file_path))
    return false;
  if (stop_condition_value <= 0 || history_size <= 0)
    return false;
  // if (simulations_configurations.empty()) return false;
  for (auto [fasta_file_path, gene_name, initiation_rate, termination_rate,
             gene_copy_number] : simulations_configurations) {
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
  auto do_simulation = [](std::string concentration_file_path,
                          bool pre_populate, std::string fasta_file,
                          std::string gene, float init_rate, float term_rate,
                          stop_condition_enum stop_condition, float stop_value,
                          std::size_t log_size) {
    // prepare and run the simulation.
    Simulations::Translation ts;
    ts.loadConcentrations(concentration_file_path);
    ts.loadMRNA(fasta_file, gene);
    ts.setInitiationRate(init_rate);
    std::cout << "init rate = " << init_rate << "\n";
    ts.setTerminationRate(term_rate);
    std::cout << "term rate = " << term_rate << "\n";
    if (stop_condition == ITERATION) {
      ts.setIterationLimit(stop_value);
      std::cout << "ITERATION\n";
    } else if (stop_condition == TIME) {
      ts.setTimeLimit(stop_value);
      std::cout << "TIME\n";
    } else if (stop_condition == RIBOSOMES) {
      ts.setFinishedRibosomes(stop_value);
      std::cout << "RIBOSOMES\n";
    }
    ts.setPrepopulate(pre_populate); // simulations pre-populate the mRNA
                                     // by default. This can be changed in
                                     // the future.
    ts.setHistorySize(log_size);
    ts.run();
    return ts;
  };

  for (std::size_t i = 0; i < number_of_simulations; i++) {
    std::string fasta_path, gene_name;
    float init_rate, term_rate, copy_number;
    std::tie(fasta_path, gene_name, init_rate, term_rate, copy_number) =
        simulations_configurations[i];
    simulations.push_back(std::async(do_simulation, concentration_file_path,
                                     pre_populate, fasta_path, gene_name,
                                     init_rate, term_rate, stop_condition_type,
                                     600, history_size));
  }

  for (auto &item : simulations) {
    results.push_back(item.get());
  }
  return true;
}

void Elongation_manager::SimulationManager::save_results() {
  std::ifstream config_doc_read(configuration_file_path, std::ifstream::binary);
  Json::Value root; // the json document.
  config_doc_read >> root;
  config_doc_read.close();

  for (auto &item : results) {
    auto ribosome_positions = item.ribosome_positions_history;
    // std::vector<float> clock(item.dt_history.size());
    // std::partial_sum(item.dt_history.begin(), item.dt_history.end(),
    //                  clock.begin(), std::plus<float>());
    Json::Value ribosome_positions_array;
    Json::Value dt_array;
    for (std::size_t i = 0; i < ribosome_positions.size(); i++) {
      Json::Value ribosome_entry;
      for (auto ribosome : ribosome_positions[i])
        ribosome_entry.append(ribosome);
      ribosome_positions_array.append(ribosome_entry);
      dt_array.append(item.dt_history[i]);
    }

    root["results"][item.gene_name]["ribosome_positions"] =
        ribosome_positions_array;
    root["results"][item.gene_name]["dt"] = dt_array;
  }
  // write in a nice readible way
  Json::StreamWriterBuilder builder;
  // builder["commentStyle"] = "None"; // no comments.
  // builder["indentation"] = "   ";   // four spaces identation
  std::unique_ptr<Json::StreamWriter> writer(builder.newStreamWriter());
  std::ofstream config_doc_writer(configuration_file_path,
                                  std::ifstream::binary);
  writer->write(root, &config_doc_writer);
}