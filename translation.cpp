#include "translation.h"

#include "mrna_reader.h"
#include "initiationterminationcodon.h"
#include <iostream>
#include <fstream>
#include <random>
#include <float.h>
#include <algorithm>
#include <deque>
/*
<%

cfg['compiler_args'] = ['-Ofast', '-ffast-math', '-flto', '-march=native']
cfg['sources'] = ['mrna_reader.cpp', 'initiationterminationcodon.cpp', 'concentrationsreader.cpp', 'ribosomesimulator.cpp', 'mrnaelement.cpp', 'gillespie.cpp', 'enlongation_codon.cpp', 'reactionsset.cpp', 'ratecalculator.cpp']
cfg['parallel'] = True

setup_pybind11(cfg)
%>
*/
#ifdef COMIPLE_PYTHON_MODULE
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;
#endif


#ifdef COMIPLE_PYTHON_MODULE
PYBIND11_MODULE(translation, mod){

    py::class_<Simulations::Translation> (mod, "translation")
    .def(py::init<>()) //constructor
    .def("loadMRNA", &Simulations::Translation::loadMRNA)
    .def("loadConcentrations", &Simulations::Translation::loadConcentrations)
    .def("setInitiationRate", &Simulations::Translation::setInitiationRate)
    .def("setTerminationRate", &Simulations::Translation::setTerminationRate)
    .def("setIterationLimit", &Simulations::Translation::setIterationLimit)
    .def("setTimeLimit", &Simulations::Translation::setTimeLimit)
    .def("setIterationLimit", &Simulations::Translation::setIterationLimit)
    .def("setFinishedRibosomes", &Simulations::Translation::setFinishedRibosomes)
    .def("run", &Simulations::Translation::run, py::call_guard<py::gil_scoped_release>())
    .def("setInitiationRate", &Simulations::Translation::setInitiationRate)
    .def("setInitiationRate", &Simulations::Translation::setInitiationRate)
    .def("getEnlogationDuration", &Simulations::Translation::getEnlogationDuration)
    .def("calculateAverageTimes", &Simulations::Translation::calculateAverageTimes)
    .def("setPrepopulate",&Simulations::Translation::setPrepopulate)
    .def("getInitiationEnlongationTermination", &Simulations::Translation::getInitiationEnlongationTermination)


    .def_readonly("mrna_file_name", &Simulations::Translation::mrna_file_name)
    .def_readonly("concentrations_file_name", &Simulations::Translation::concentrations_file_name)
    .def_readonly("dt_history", &Simulations::Translation::dt_history)
    .def_readonly("ribosome_positions_history", &Simulations::Translation::ribosome_positions_history)
    .def_readonly("initiationRate", &Simulations::Translation::initiation_rate)
    .def_readonly("terminationRate", &Simulations::Translation::termination_rate)
    
    .def_readonly("initiations_durations", &Simulations::Translation::initiations_durations)
    .def_readonly("enlongations_durations", &Simulations::Translation::enlongations_durations)
    .def_readonly("terminations_durations", &Simulations::Translation::terminations_durations)
    
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
        initiation_codon->setState(0);
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

void Simulations::Translation::setPrepopulate(bool prep)
{
    pre_populate = prep;
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


/**
* @brief Set a limit of the number of ribosomes that sucessfully initiate and terminates the mRNA.
* 
* @param n_ribosomes p_n_ribosomes:The simulation will end after this number of ribosomes terminates the mRNA.
*/
void Simulations::Translation::setFinishedRibosomes(int n_ribosomes)
{
    if (n_ribosomes > 0) finished_ribosomes_limit = n_ribosomes;
}



void Simulations::Translation::getAlphas()
{
    alphas.clear();
    codon_index.clear();
    reaction_index.clear();
    std::vector<double> a;
    std::vector<int>r_i;

    //populate the vectors.
    std::vector<int> ribosome_positions = ribosome_positions_history.back();
    std::size_t ribosome_index;
    //add initiation if needed.
    if (ribosome_positions.empty() || (ribosome_positions[0] > 9 && codons_vector[0]->isAvailable)) {
        //need to add initalization.
        codons_vector[0]->getAlphas(a, r_i);
        alphas.insert(alphas.end(), a.begin(), a.end());
        std::fill_n(std::back_inserter(codon_index), a.size(), 0);
        reaction_index.insert(reaction_index.end(), r_i.begin(), r_i.end());
    }

    for (unsigned i = 0; i < ribosome_positions.size(); i++){
        ribosome_index = ribosome_positions[i];
        codons_vector[ribosome_index]->getAlphas(a, r_i);
        alphas.insert(alphas.end(), a.begin(), a.end());
        std::fill_n(std::back_inserter(codon_index), a.size(), ribosome_index);
        reaction_index.insert(reaction_index.end(), r_i.begin(), r_i.end());
    }
}


void Simulations::Translation::run()
{
    dt_history = std::vector<double>(iteration_limit > 0 ? iteration_limit:100000);
    dt_history.clear();
    
    ribosome_positions_history = std::vector<std::vector<int>>(iteration_limit > 0 ? iteration_limit: 100000);
    ribosome_positions_history.clear();
    
    // initialize the random generator
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(DBL_MIN, 1);

    double r1 = 0, r2 = 0;
    double tau = 0, clock = 0.0;
    int i = 0;
    
    int finished_ribosomes = 0, pre_filled_ribosomes = 0;
    std::vector<int> rib_positions((codons_vector.size()/10) + 1);
    rib_positions.clear();
    //pre-fill codons based on the rates.
    if (pre_populate){
        std::vector<double> a;
        std::vector<int>r_i;
        codons_vector[0]->getAlphas(a, r_i);
        double initiation_time = 1/a[0]; //propensity
        int last_index = codons_vector.size() - 1;
        double time_sum = 0;
        codons_vector[last_index]->isOccupied = true;
        codons_vector[last_index]->isAvailable = false;
        codons_vector[last_index]->setState(0);
        pre_filled_ribosomes++;
        for (int i = codons_vector.size() - 2; i >= 0; i--) {
            if (last_index - i <= 9){
                codons_vector[i]->isAvailable = false;
            } else {
                //update ribosome.
                codons_vector[i]->getAlphas(a, r_i);
                time_sum += 1/a[0];
            }
            if (time_sum >= initiation_time) {
                //put a ribosome here.
                codons_vector[i]->isOccupied = true;
                codons_vector[i]->isAvailable = false;
                codons_vector[i]->setState(0);
                if (i == 0) codons_vector[i]->setState(23);
                time_sum = 0; //reset timer.
                last_index = i; //mark this as last inserted ribosome.
                pre_filled_ribosomes++;
            }
        }
    }
    for (unsigned int i = 0; i < codons_vector.size(); i++) {
        if (codons_vector[i]->isOccupied)  rib_positions.push_back(i);
    }
    
    for (unsigned int i = 0; i < rib_positions.size(); i++) {
        if ((unsigned) rib_positions[i] < codons_vector.size() - 1) {
            codons_vector[rib_positions[i]]->updateAlphas(codons_vector[rib_positions[i] + 1]->isAvailable); //update reactions.
        }
        //also need to update reactions of each previous codon! Otherwise, a ribosome will NOT collide.
        if (rib_positions[i]  > 8) codons_vector[rib_positions[i] - 9]->updateAlphas(false); //update reactions.
    }
    
    finished_ribosomes -= pre_filled_ribosomes; // we should ignore these ribosomes.
    int moved_codon = -1, current_codon=-1;
    bool initiation = false, termination = false, moved = true;
    
    while ((iteration_limit > 0  && i < iteration_limit) || (time_limit > 0 && clock < time_limit) || (finished_ribosomes_limit > 0 && finished_ribosomes_limit >= finished_ribosomes))
    {
        //get the vector with the positions of all ribosomes
        if (i > 0 && moved) {
            if (termination){
                //terminated. remove last position.
                rib_positions.pop_back();
            } else if (initiation) {
                //initiated.
                rib_positions.insert(rib_positions.begin(), 0);
            } else if (moved) {
                auto i = std::find(rib_positions.begin(), rib_positions.end(), moved_codon - 1);
                auto pos = std::distance(rib_positions.begin(), i);
                rib_positions[pos] = moved_codon;
            }
        }
        if (!moved){
            //no ribosome movement. just update dt_history.
            dt_history.back() +=tau;
        } else {
            // ribosome movement detected. create new entry in the history.
            dt_history.push_back(tau);
            ribosome_positions_history.push_back(rib_positions);
        }

        moved = false;
        initiation = false;
        termination = false;
        moved_codon = -1;
        // randomly generate parameter for calculating dt
        r1 = dis(gen);
        // randomly generate parameter for selecting reaction
        r2 = dis(gen);
        // calculate an
        getAlphas();
        if (alphas.empty())
        {
            // no available reactions, quit loop prematurely.
            std::cout<<"no available reactions. quitting.\n";
            break;
        }
        double a0 = std::accumulate(alphas.begin(), alphas.end(), 0.0);
        int selected_alpha_vector_index = -1;
        // The commented code below is the vectorized version of the reaction selection.
        // Here is the catch: the vectorized version always calculates all reactions, and from them, select one.
        // the non-vectorized one (do... while loop) calculates only the ones up to the selected one.
//         std::vector<double> cumsum(alphas.size());
//         std::partial_sum(alphas.begin(), alphas.end(), cumsum.begin());
//         selected_alpha_vector_index = std::distance(cumsum.begin(), std::upper_bound(cumsum.begin(), cumsum.end(), a0 * r2));
        
        // select next reaction to execute
        double cumsum = 0;
        do {
            selected_alpha_vector_index++;
            cumsum += alphas[selected_alpha_vector_index]; 
        } while (cumsum < a0 * r2);
        current_codon = codon_index[selected_alpha_vector_index];
        //Apply reaction
        codons_vector[current_codon]->executeReaction(reaction_index[selected_alpha_vector_index]);
        if ((unsigned) (current_codon + 1) <  codons_vector.size()) codons_vector[current_codon]->updateAlphas(codons_vector[current_codon + 1]->isAvailable); // update current codon.

        if (codon_index[selected_alpha_vector_index] == 0 && codons_vector[0]->getState() == 23){
            //initiated.
            codons_vector[0]->isAvailable=false;
            codons_vector[0]->isOccupied=true;
            initiation = true;
            moved = true;
        }
        moved_codon = current_codon + 1;
        // 2- Any codon with state == 31 means the ribosome already moved to the next codon (or left the mRNA). update states.
        if (codons_vector[current_codon]->getState() == 31) {
            codons_vector[current_codon]->setState(0);
            codons_vector[current_codon]->isAvailable = false;
            codons_vector[current_codon]->isOccupied = false;
            codons_vector[current_codon]->updateAlphas(false); //cannot translocate.
            moved = true;
            moved_codon = current_codon + 1;
            if ((unsigned) (moved_codon + 1) <  codons_vector.size()) {
                codons_vector[moved_codon]->isOccupied = true;
                codons_vector[moved_codon]->isAvailable = false;
                codons_vector[moved_codon]->updateAlphas(codons_vector[moved_codon + 1]->isAvailable); // update current codon.
            }
            //update free codons due to the size of the ribosome.
            if ((unsigned) moved_codon > 9){
                //we need to do some tidying up after the ribosome.
                if ((unsigned) moved_codon < codons_vector.size()){
                    //update freed space left by the ribosome's movement.
                    codons_vector[moved_codon - 10]->isAvailable = true;
                    codons_vector[moved_codon - 10]->updateAlphas(codons_vector[moved_codon - 9]->isAvailable);
                    if ((unsigned) moved_codon > 10) {
                        codons_vector[moved_codon - 11]->updateAlphas(true); // necessarily available.
                    }
                } else if ((unsigned) moved_codon == codons_vector.size()) {
                    //ribosome terminated. free codons positions occupied by it.
                    termination = true;
                    codons_vector[current_codon]->updateAlphas(true); // update codon: it can always translocate.
                    for (unsigned int i = codons_vector.size() - 10; i < codons_vector.size(); i++) codons_vector[i]->isAvailable = true;
                    if (codons_vector.size() > 10) codons_vector[codons_vector.size() - 11]->updateAlphas(codons_vector[codons_vector.size() - 10]->isAvailable);
                    //if the simulation has a terminated ribosome limit, update the variable.
                    if (finished_ribosomes_limit > 0){
                        finished_ribosomes++;
                    }
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

std::tuple<std::vector<double>, std::vector<int>> Simulations::Translation::getEnlogationDuration()
{
    if (enlongations_durations.empty()) getInitiationEnlongationTermination();
    return std::make_tuple(enlongations_durations, initiation_iteration);
}

void Simulations::Translation::getInitiationEnlongationTermination()
{
    initiations_durations.clear();
    enlongations_durations.clear();
    terminations_durations.clear();
    initiation_iteration.clear();

    std::deque<int> indexes;
    indexes.clear();
    std::list<int> initiations, enlongations, terminations;
    int ribosomes_to_ignore = ribosome_positions_history[0].size();
    unsigned int last_position = codons_vector.size() - 1, previous_size = ribosomes_to_ignore;
    for (int i = 0; i < ribosomes_to_ignore; i++) {
        indexes.push_back(indexes.size());
        initiations_durations.push_back(0);
        enlongations_durations.push_back(0);
        terminations_durations.push_back(0);
        initiation_iteration.push_back(0);
    }
    for (unsigned int i = 1; i < ribosome_positions_history.size(); i++) {
        std::vector<int> &rib_positions = ribosome_positions_history[i];
        for (unsigned int j = 0; j <rib_positions.size(); j++){
            unsigned int pos = rib_positions[j];
            if (pos == 0){
                //initiating.
                if (rib_positions.size() > previous_size){
                    //adjust offset for the iteration;
                    indexes.push_front(indexes.size());
                    //new ribosome initiating.
                    initiations_durations.push_back(0);
                    enlongations_durations.push_back(0);
                    terminations_durations.push_back(0);
                    initiation_iteration.push_back(i);
                }
                initiations_durations[indexes[j]] += dt_history[i];
            } else if (pos == last_position){
                //started terminating
                terminations_durations[indexes[j]] += dt_history[i];
            } else {
                //enlongating codon.
                enlongations_durations[indexes[j]] += dt_history[i];
            }
        }
        previous_size = rib_positions.size();
    }
    if (ribosomes_to_ignore > 0){
        //remove pre-filled ribosomes.
        initiations_durations.erase(initiations_durations.begin(), initiations_durations.begin() + ribosomes_to_ignore);
        enlongations_durations.erase(enlongations_durations.begin(), enlongations_durations.begin() + ribosomes_to_ignore);
        terminations_durations.erase(terminations_durations.begin(), terminations_durations.begin() + ribosomes_to_ignore);
        initiation_iteration.erase(initiation_iteration.begin(), initiation_iteration.begin() + ribosomes_to_ignore);
    }
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
    for (auto ribosome_vector:ribosome_positions_history) {
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
