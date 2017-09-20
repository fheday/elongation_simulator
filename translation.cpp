#include "translation.h"

#include "mrna_reader.h"
#include "initiationterminationcodon.h"
#include <iostream>
#include <fstream>
#include <random>
#include <float.h>
/*
<%

cfg['compiler_args'] = ['-O3', '-ffast-math']
cfg['sources'] = ['mrna_reader.cpp', 'initiationterminationcodon.cpp', 'concentrationsreader.cpp', 'ribosomesimulator.cpp', 'mrnaelement.cpp', 'gillespie.cpp', 'enlongation_codon.cpp', 'reactionsset.cpp', 'ratecalculator.cpp']
cfg['parallel'] = True

setup_pybind11(cfg)
%>
*/
#ifndef CMAKE_BUILD
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

#endif


#ifndef CMAKE_BUILD
PYBIND11_MODULE(translation, mod){

    py::class_<Simulations::Translation> (mod, "translation")
    .def(py::init<>()) //constructor
    .def("loadMRNA", &Simulations::Translation::loadMRNA)
    .def("loadConcentrations", &Simulations::Translation::loadConcentrations)
    .def("setInitiationRate", &Simulations::Translation::setInitiationRate)
    .def("setTerminationRate", &Simulations::Translation::setTerminationRate)
    .def("setIterationLimit", &Simulations::Translation::setIterationLimit)
    .def("setTimeLimit", &Simulations::Translation::setTimeLimit)
    .def("run", &Simulations::Translation::run)
    .def("setInitiationRate", &Simulations::Translation::setInitiationRate)
    .def("setInitiationRate", &Simulations::Translation::setInitiationRate)
    .def("getEnlogationDuration", &Simulations::Translation::getEnlogationDuration)
    .def("calculateAverageTimes", &Simulations::Translation::calculateAverageTimes)


    .def_readonly("mrna_file_name", &Simulations::Translation::mrna_file_name)
    .def_readonly("concentrations_file_name", &Simulations::Translation::concentrations_file_name)
    .def_readonly("dt_history", &Simulations::Translation::dt_history)
    .def_readonly("ribosome_positions_history", &Simulations::Translation::ribosome_positions_history)
    
    .def_readonly("total_time", &Simulations::Translation::total_time)
    .def_readonly("n_times_occupied", &Simulations::Translation::n_times_occupied)
    .def_readonly("average_times", &Simulations::Translation::codons_average_occupation_time);

}
#endif

void Simulations::Translation::loadConcentrations(std::string file_name)
{
    std::ifstream ist{file_name};

    if (!ist) {
        throw std::runtime_error("can't open input file: "+ file_name);
    } else {
        concentrations_file_name = file_name;
        initializeMRNAReader();
    }
}

void Simulations::Translation::loadMRNA(std::string file_name)
{
    std::ifstream ist{file_name};

    if (!ist) {
        throw std::runtime_error("can't open input file: "+ file_name);
    } else {
        mrna_file_name = file_name;
        initializeMRNAReader();
    }
}

void Simulations::Translation::initializeMRNAReader()
{
    if (!concentrations_file_name.empty() && !mrna_file_name.empty() && initiation_rate >0 && termination_rate > 0) {
        //we have the concentrations and mrna file names. we can proceed.
        mRNA_utils::mRNAReader mrr;
        mrr.loadmRNAFile(mrna_file_name);
        //fill codon vector.
        int n_codons = mrr.mRNA_sequence.size()/3;
        std::unique_ptr<Simulations::InitiationTerminationCodon> initiation_codon(new Simulations::InitiationTerminationCodon(initiation_rate, true));
        initiation_codon->codon = mrr.getCodon(0);
        initiation_codon->index = 0;
        codons_vector.push_back(std::move(initiation_codon));
        for (int i = 1; i < n_codons - 1; i++) {
            std::unique_ptr<Simulations::EnlongationCodon> c(new Simulations::EnlongationCodon());
            c->loadConcentrations(concentrations_file_name);
            c->codon = mrr.getCodon(i);
            c->setCodon(c->codon);
            c->index = i;
            codons_vector.push_back(std::move(c));
        }
        
        //termination codon
        std::unique_ptr<Simulations::InitiationTerminationCodon> termination_codon(new Simulations::InitiationTerminationCodon(termination_rate, false));
        termination_codon->codon = mrr.getCodon(n_codons - 1);
        
        termination_codon->index = n_codons - 1;
        codons_vector.push_back(std::move(termination_codon));
    }
}


void Simulations::Translation::setInitiationRate(double ir)
{
    if (ir > 0){
        initiation_rate = ir;
    }
    initializeMRNAReader();
}

void Simulations::Translation::setTerminationRate(double tr)
{
    if (tr > 0) {
        termination_rate = tr;
    }
    initializeMRNAReader();
}

/**
 * @brief Set a iteration limit for the Gillespie simulation.
 * 
 * @param i integer with the maximum number of iterations. The algorithm halts before this condition is met if there are no possible reations left to be performed.
 */
void Simulations::Translation::setIterationLimit(int i)
{
    if (i > 0) iteration_limit = i;
}

/**
 * @brief Set a time limit for the Gillespie simulation. This time is in seconds, and it is compared against the simulation's clock.
 * 
 * @param t time limit in seconds.
 */

void Simulations::Translation::setTimeLimit(double t)
{
    if (t > 0) time_limit = t;
}


void Simulations::Translation::getAlphas()
{
    alphas = std::vector<double>();
    alphas.clear();
    codon_index = std::vector<int>();
    codon_index.clear();
    codon_index = std::vector<int>();
    reaction_index.clear();
    
    //populate the vectors.
    for (unsigned int i = 0; i < codons_vector.size(); i++){
        if ((i==0 && codons_vector[0]->isAvailable == true) || codons_vector[i]->isOccupied) {
            std::vector<double> as;
            std::vector<int> reactions_indexes;
            codons_vector[i]->getAlphas(as, reactions_indexes);
            //check: in case of translocation, the next ribosome must be AVAILABLE.
            for (int j = 0; j < as.size(); j++) {
                if (as[j] <=22) {
                    //still decoding. add.
                    alphas.push_back(as[j]);
                    codon_index.push_back(i);
                    reaction_index.push_back(reactions_indexes[j]);
                } else {
                    if (as[j] == 24 && i <= codons_vector.size() - 1 && !codons_vector[i + 1]->isAvailable) {
                        continue;
                    } 
                    if (i == codons_vector.size() - 1 || codons_vector[i + 1]->isAvailable){
                        //translocating. only add if next codon is available Or if termination codon.
                        alphas.push_back(as[j]);
                        codon_index.push_back(i);
                        reaction_index.push_back(reactions_indexes[j]);
                    }
                }
            }
        }
    }
}


void Simulations::Translation::run()
{
    dt_history = std::vector<double>(iteration_limit > 0 ? iteration_limit:100000);
    dt_history.clear();
    
    ribosome_positions_history = std::vector<std::vector<int>>(iteration_limit > 0 ? iteration_limit: 100000);
    ribosome_positions_history.clear();
    
    //Eigen::MatrixXi updated_populations;
    
    // initialize the random generator
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(0, 1);

    double r1 = 0, r2 = 0;
    double tau = 0, clock = 0.0;
    int i = 0;

    while ((iteration_limit > 0  && i < iteration_limit) || (time_limit > 0 && clock < time_limit))
    {
        //get the vector with the positions of all ribosomes
        std::vector<int> rib_positions;
        for (unsigned int i = 0; i < codons_vector.size(); i++){
            if (codons_vector[i]->isOccupied) rib_positions.push_back(i);
        }
        if (rib_positions.size() > 0 && rib_positions.size() == ribosome_positions_history.back().size() && std::equal(rib_positions.begin(), rib_positions.end(), ribosome_positions_history.back().begin()) ){
            //no ribosome movement. just update dt_history.
            dt_history.back() +=tau;
        } else {
            // ribosome movement detected. create new entry in the history.
            dt_history.push_back(tau);
            ribosome_positions_history.push_back(rib_positions);
        }
        // randomly generate parameter for calculating dt
        r1 = dis(gen) + DBL_MIN; // adding minumum double value in order to avoid division by zero and infinities.
        // randomly generate parameter for selecting reaction
        r2 = dis(gen)+ DBL_MIN; // adding minumum double value in order to avoid division by zero and infinities.
        // calculate an
        getAlphas();
        if (alphas.empty())
        {
            // no available reactions, quit loop prematurely.
            break;
        }
        double a0 = std::accumulate(alphas.begin(), alphas.end(), 0);
        // select next reaction to execute
        double cumsum = 0;
        int selected_alpha_vector_index = -1;
        // TODO: vectorization of this loop would increase performance
        do {
            selected_alpha_vector_index++;
            cumsum += alphas[selected_alpha_vector_index]; 
        } while (cumsum < a0 * r2);
        //Apply reaction
        codons_vector[codon_index[selected_alpha_vector_index]]->executeReaction(reaction_index[selected_alpha_vector_index]);
        // 2- Any codon with state == 31 means the ribosome already moved to the next codon (or left the mRNA). update states.
        if (codons_vector[codon_index[selected_alpha_vector_index]]->getState() == 31) {
            codons_vector[codon_index[selected_alpha_vector_index]]->setState(0);
            codons_vector[codon_index[selected_alpha_vector_index]]->isOccupied = false;
            if ((unsigned) (codon_index[selected_alpha_vector_index] + 1) <  codons_vector.size()){
                //there is a next codon. move on to it.
                codons_vector[codon_index[selected_alpha_vector_index]+ 1]->isOccupied = true;
                codons_vector[codon_index[selected_alpha_vector_index]+ 1]->isAvailable = false;
            }
            //update free codons due to the size of the ribosome.
            if ((unsigned) codon_index[selected_alpha_vector_index] > 8 && (unsigned) codon_index[selected_alpha_vector_index] < codons_vector.size() - 1){
                codons_vector[codon_index[selected_alpha_vector_index] - 9]->isAvailable = true;
            } else if ((unsigned) codon_index[selected_alpha_vector_index] == codons_vector.size() - 1) {
                //ribosome terminated. free codons positions occupied by it.
                for (unsigned int i = codons_vector.size() - 10; i < codons_vector.size(); i++){
                    codons_vector[i]->isAvailable = true;
                }
            }
        }

        //Update time
        tau = (1.0/a0) * log(1.0/r1);
        clock += tau;
        i++; // update iteration number.
    }

}

/**
 * @brief Returns a tuple where the first element is a vector with the enlogation duration of the ribosomes that terminated in the simulation, and the second element is a vector with the iteration where such ribosomes started enlogating. This method should be called after updateRibosomeHistory, since it uses the positions_vector to do its job.
 * 
 */

std::tuple<std::vector<double>, std::vector<int>, std::vector<int>> Simulations::Translation::getEnlogationDuration()
{
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
    int last_codon_position = codons_vector.size() - 1;
    for (unsigned int i = 1; i < ribosome_positions_history.size(); i++)
    {
        std::vector<int>& ribosomes_positions = ribosome_positions_history[i];
        current_zero_occupied = std::find(ribosomes_positions.begin(), ribosomes_positions.end(), 0) !=ribosomes_positions.end();
        current_last_occupied = std::find(ribosomes_positions.begin(), ribosomes_positions.end(), last_codon_position) !=ribosomes_positions.end();
        if (!previous_zero_occupied && current_zero_occupied){
            //new ribosome.
            rib_initiation_iteration.push_back(i);
        } else if (previous_last_occupied && !current_last_occupied){
            //ribosome left.
            result.push_back(clock[i - 1] - clock[rib_initiation_iteration.front()]);
            initiation_iteration.push_back(rib_initiation_iteration.front());
            termination_iteration.push_back(i - 1);
            rib_initiation_iteration.pop_front();
        }
        //update variables for next iteration.
        previous_last_occupied = current_last_occupied;
        previous_zero_occupied = current_zero_occupied;
    }
    return std::make_tuple(result, initiation_iteration, termination_iteration);
}


void Simulations::Translation::calculateAverageTimes()
{
    int number_codons = codons_vector.size();
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
    for (std::vector<int> ribosome_vector:ribosome_positions_history) {
        for (int position:ribosome_vector) {
            total_time[position] += dt_history[iteration_number];
            if (last_index_occupied[position] == -1 || last_index_occupied[position] != iteration_number - 1){
                //we are facing a re-entering ribosome. We need to add the previous occupation.
                n_times_occupied[position]++;
            }
            //update the last time this position was occupied.
            last_index_occupied[position] = iteration_number;
        }
        iteration_number++;
    }
    //the above procedure does not count for the last time a position has been occupied: it ignores it.
    //we could try to fix this in a number of ways, but I guess it wouldn't matter much in the big picture.

    //now we calculate the averages.
    codons_average_occupation_time.clear();
    codons_average_occupation_time = std::vector<double>(number_codons);
    for (int codon_position = 0; codon_position < number_codons; codon_position++) {
        codons_average_occupation_time[codon_position] = n_times_occupied[codon_position] > 0 ? total_time[codon_position] / n_times_occupied[codon_position]: 0;
    }
}
