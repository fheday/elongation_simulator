#include "elongation_simulation_processor.h"
#include <fstream>
#include <iostream>
#include "json/json.h"

#if defined(COMIPLE_PYTHON_MODULE) || defined(TRANSLATIONSIMULATOR)

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

void init_simulation_processor(py::module &mod) 
 {
  py::class_<Simulations::SimulationProcessor>(mod, "SimulationProcessor")
      .def(py::init<std::string>()) // constructor
      .def("getClock", &Simulations::SimulationProcessor::getClock)
      .def("getRibosomes", &Simulations::SimulationProcessor::getRibosomes);
}
#endif


Simulations::SimulationProcessor::SimulationProcessor(std::string file_name) {
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

