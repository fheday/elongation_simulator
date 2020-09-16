#include "elongation_simulation_processor.h"
#include <fstream>
#include <iostream>
#include "json/json.h"

#if defined(COMIPLE_PYTHON_MODULE) || defined(TRANSLATIONSIMULATOR)

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

void init_simulation_processor(py::module &mod) {
  py::class_<Simulations::SimulationProcessor>(mod, "SimulationProcessor")
      .def(py::init<std::string>()) // constructor
      .def("getClock", &Simulations::SimulationProcessor::getClock)
      .def("getRibosomes", &Simulations::SimulationProcessor::getRibosomes)
      .def("removeRibosomePositions", &Simulations::SimulationProcessor::removeRibosomePositions)
      .def("calculateRibosomeCollisions", &Simulations::SimulationProcessor::calculateRibosomeCollisions)
      .def("getCollidingRibosomes", &Simulations::SimulationProcessor::getCollidingRibosomes)
      .def("getStalledRibosomes", &Simulations::SimulationProcessor::getStalledRibosomes)
      .def("save", &Simulations::SimulationProcessor::save);
}
#endif


Simulations::SimulationProcessor::SimulationProcessor(std::string file_name) {
    configuration_file_name = file_name;
    std::ifstream config_doc(file_name, std::ifstream::binary);
    Json::Value root; // the json document.
    config_doc >> root;
    if (root.isMember("fasta_file")) fasta_file = root.get("fasta_file", "").asString();
    if (root.isMember("initiation_rate")) initiation_rate = root.get("initiation_rate", -1).asFloat();
    if (root.isMember("termination_rate")) termination_rate = root.get("termination_rate", -1).asFloat();
    if (!root.isMember("clock")) {
        std::cout<<"Error. simulation file has no clock information.\n";
        return;
    }
    if (!root.isMember("elongating_ribosomes")) {
        std::cout<<"Error. simulation file has no ribosomes position information.\n";
        return;
    }
    if (root.isMember("clock")){
        for (unsigned int i = 0; i < root["clock"].size(); ++i) {
            clock.emplace_back(root["clock"][i].asFloat());
        }
    }
    if (root.isMember("elongating_ribosomes")){
        for (unsigned int i = 0; i < root["elongating_ribosomes"].size(); ++i) {
            std::vector<int> entry(root["elongating_ribosomes"][i].size());
            for (unsigned int j = 0; j < root["elongating_ribosomes"][i].size(); ++j) {
                entry[j] = root["elongating_ribosomes"][i][j].asInt();
            }
            elongating_ribosomes.emplace_back(entry);
        }
    }

}

std::vector<float>& Simulations::SimulationProcessor::getClock() {
    return clock;
}

std::vector<std::vector<int>>& Simulations::SimulationProcessor::getRibosomes() {
    return elongating_ribosomes;
}

void Simulations::SimulationProcessor::removeRibosomePositions() {
    elongating_ribosomes.clear();
}

void Simulations::SimulationProcessor::calculateRibosomeCollisions() {
    for (auto ribosomes_positions:elongating_ribosomes) {
        std::vector<int> collision_entry;
        std::vector<int> stall_entry;
        for (std::size_t i = 0; i < ribosomes_positions.size() - 1; i++) {
            if (ribosomes_positions[i + 1] - ribosomes_positions[i] == 10) {
                // colliding ribosome: collision with the next ribosome detected.
                collision_entry.emplace_back(ribosomes_positions[i]);
            } else if (!collision_entry.empty() && ribosomes_positions[i] - collision_entry.back() == 10) {
                //stalled ribosome: no collision with next ribosome,
                //but collision with previous ribosome detected.
                stall_entry.emplace_back(ribosomes_positions[i]);
            }
        }
        // check last entry. it can only stall.
        if (!collision_entry.empty() && ribosomes_positions.back() - collision_entry.back() == 10) {
            //stalled ribosome: no collision with next ribosome,
            //but collision with previous ribosome detected.
            stall_entry.emplace_back(ribosomes_positions.back());
        }
        colliding_ribosomes.emplace_back(collision_entry);
        stalled_ribosomes.emplace_back(stall_entry);
    }
}

std::vector<std::vector<int>>& Simulations::SimulationProcessor::getCollidingRibosomes() {
    return colliding_ribosomes;
}

std::vector<std::vector<int>>& Simulations::SimulationProcessor::getStalledRibosomes() {
    return stalled_ribosomes;
}

void Simulations::SimulationProcessor::save() {
    Json::Value newjson;
    newjson["fasta_file"] = fasta_file;
    newjson["initiation_rate"] = initiation_rate;
    newjson["termination_rate"] = termination_rate;
    
    auto generate_json_vector_of_vector = [&](auto data_vector) {
        Json::Value json_value;
        for (auto entry:data_vector){
            Json::Value entry_vector;
            for (auto element:entry) entry_vector.append(element);
            json_value.append(entry_vector);
        }
        return json_value;
    };


    for (auto time:clock) newjson["clock"].append(time);

    newjson["elongating_ribosomes"] = generate_json_vector_of_vector(elongating_ribosomes);
    newjson["colliding_ribosomes"] = generate_json_vector_of_vector(colliding_ribosomes);
    newjson["stalling_ribosomes"] = generate_json_vector_of_vector(stalled_ribosomes);

    // std::ifstream file_stream(configuration_file_name, std::ifstream::binary);
    // write in a nice readible way
    Json::StreamWriterBuilder builder;
    builder["commentStyle"] = "None"; // no comments.
    builder["indentation"] = "   ";   // four spaces identation
    std::unique_ptr<Json::StreamWriter> writer(builder.newStreamWriter());
    std::ofstream config_doc_writer(configuration_file_name,
                                        std::ifstream::binary);
    writer->write(newjson, &config_doc_writer);
    return;
}
