#include <fstream>
#include "enlogationsimulator.h"

using namespace Simulations;

EnlogationSimulator::EnlogationSimulator()
{

}

void EnlogationSimulator::set_initiation_rate(double ir)
{
    if (ir > 0) {
        initiation_rate = ir;
    } else {
        throw std::runtime_error("invalid initiation rate: " + std::to_string(ir));
    }
    intializeMRNAReader();
}

void EnlogationSimulator::set_termination_rate(double tr)
{
    if (tr > 0) {
        termination_rate = tr;
    } else {
        throw std::runtime_error("invalid termination rate: " + std::to_string(tr));
    }
    intializeMRNAReader();
}

void EnlogationSimulator::set_mRna_file_name(std::string file_name)
{
    std::ifstream ist{file_name};

    if (!ist) {
        throw std::runtime_error("can't open input file: "+ file_name);
    } else {
        mRNAFileName = file_name;
    }
    intializeMRNAReader();
}

void EnlogationSimulator::intializeMRNAReader()
{
    if (!mRNAFileName.empty() and !average_times_file_name.empty() && initiation_rate > 0 && termination_rate > 0) {
        // we can proceed with the mRNAReader object.
        mrna_reader.loadRateCalculatorFile(average_times_file_name);
        mrna_reader.loadmRNAFile(mRNAFileName);
        mrna_reader.setInitiationRate(initiation_rate);
        mrna_reader.setTerminationRate(termination_rate);
        mrna_reader.generateInitialPopulation();
        mrna_reader.generateReactions();
        Gillespie::setInitialPopulation(mrna_reader.initial_population);
        Gillespie::setReactionsSet(mrna_reader.reactions_set);
    }
}

void EnlogationSimulator::setAverageTimesFileName(std::string file_name)
{
    std::ifstream ist{file_name};

    if (!ist) {
        throw std::runtime_error("can't open input file: "+ file_name);
    } else {
        average_times_file_name = file_name;
    }
    intializeMRNAReader();

}


void EnlogationSimulator::set_concentrations_file_name(std::string file_name)
{
    std::ifstream ist{file_name};

    if (!ist) {
        throw std::runtime_error("can't open input file: "+ file_name);
    } else {
        concentrationsFileName = file_name;
        // when setting the concentrations file name, we can also
        // initialize the RibosomeSimulator object.
        ribosome_simulator.loadConcentrations(file_name);
        ribosome_simulator.setIterationLimit(2000);
        ribosome_simulator.setNumberOfRibosomes(1);
    }
}


double EnlogationSimulator::getReactionTime(double a0, double r1, std::string codon)
{
    double result = 0;
    if (codon == "tra") {
        result = translocation_times.back();
        translocation_times.pop_back();
    } else {
        //calculate the time.
        double translocating;
        ribosome_simulator.setCodonForSimulation(codon);
        ribosome_simulator.run_and_get_times(result, translocating);
        //add the translocating time to our pool.
        translocation_times.push_back(translocating);
    }
    return result;
}

