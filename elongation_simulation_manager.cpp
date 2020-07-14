#include "elongation_simulation_manager.h"
#include <fstream>
#include "jsoncpp.cpp"

Elongation_manager::SimulationManager::SimulationManager(std::string cfp)
{
    configuration_file_path = cfp;
    std::ifstream config_doc(cfp, std::ifstream::binary);
    Json::Value root;
    config_doc >> root;
    concentration_file_path = root.get("concentration_file", "missing").asString();
    pre_populate = root.get("pre_populate", false).asBool();

    std::string fasta_file_path;
    std::string gene_string;
    float init_rate = 0, term_rate = 0, gene_copy_number = 0;
    for (unsigned int i = 0; i < root["mRNA_entries"].size(); ++i)
    {
        fasta_file_path = root["mRNA_entries"][i]["fasta_file"].asString();
        gene_string = root["mRNA_entries"][i]["gene"].asString();
        init_rate = root["mRNA_entries"][i]["initiation_rate"].asFloat();
        term_rate = root["mRNA_entries"][i]["termination_rate"].asFloat();
        gene_copy_number = root["mRNA_entries"][i]["transcript_copy_number"].asFloat();
        simulations_configurations.push_back(std::tuple(fasta_file_path, gene_string, init_rate, term_rate, gene_copy_number));
    }

    if (root.isMember("iteration_limit"))
    {
        stop_condition_type = ITERATION;
        stop_condition_value = root.get("iteration_limit", "-1").asFloat();
    }
    else if (root.isMember("time_limit"))
    {
        stop_condition_type = TIME;
        stop_condition_value = root.get("time_limit", "-1").asFloat();
    }
    else if (root.isMember("finished_ribosomes"))
    {
        stop_condition_type = RIBOSOMES;
        stop_condition_value = root.get("finished_ribosomes", "-1").asFloat();
    }
    history_size = root.get("history_size", 10000).asUInt();
    config_doc.close(); // close file
    if (!is_simulation_valid()) std::cout<<"Error in configuration\n";
};

bool Elongation_manager::SimulationManager::is_simulation_valid()
{
    auto file_exists = [](std::string file_path)
    {
        // check if files exist and if numerical values are valid.
        std::ifstream conf_file_path{file_path};

        if (!conf_file_path)
        {
            throw std::runtime_error("can't open input file: " + file_path);
            return false;
        }

        return true;
    };
    if (!file_exists(configuration_file_path)) return false;
    if (!file_exists(concentration_file_path)) return false;
    if (stop_condition_value <= 0 || history_size <= 0) return false;
    // if (simulations_configurations.empty()) return false;
    for (auto [fasta_file_path, gene_name, initiation_rate, termination_rate, gene_copy_number] : simulations_configurations) {
        if (!file_exists(fasta_file_path)) return false;
        if (gene_name.empty()) return false;
        if (initiation_rate < 0 || termination_rate < 0 || gene_copy_number <= 0) return false;
    }
    return true;
}

std::string Elongation_manager::SimulationManager::get_concentration_file_path()
{
    return concentration_file_path;
}

std::string Elongation_manager::SimulationManager::get_configuration_file_path()
{
    return configuration_file_path;
}

bool Elongation_manager::SimulationManager::get_pre_populate()
{
    return pre_populate;
}

std::vector<std::tuple<std::string, std::string, float, float, float>> &Elongation_manager::SimulationManager::get_simulations_configurations()
{
    return simulations_configurations;
}

Elongation_manager::stop_condition_enum Elongation_manager::SimulationManager::get_stop_condition_type()
{
    return stop_condition_type;
}

float Elongation_manager::SimulationManager::get_stop_condition_value()
{
    return stop_condition_value;
}

std::size_t Elongation_manager::SimulationManager::get_history_size()
{
    return history_size;
}