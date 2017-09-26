/*
<%
cfg['dependencies'] = ['reactionsset.h', 'gillespie.h']
cfg['compiler_args'] = ['-O3']
cfg['sources'] = ['reactionsset.cpp', 'gillespie.cpp', 'concentrationsreader.cpp']

setup_pybind11(cfg)
%>
*/

#ifndef CMAKE_BUILD
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

#endif

#include "ribosomesimulator.h"
#include <eigen3/Eigen/Dense>
#include "concentrationsreader.h"
#include <numeric>
#include <float.h>
#include <float.h>

using namespace Simulations;

#ifndef CMAKE_BUILD
PYBIND11_MODULE(ribosomesimulator, mod){

    py::class_<Gillespie> (mod, "gillespie")
    .def(py::init<>()) //constructor
    .def("setIterationLimit", &Gillespie::setIterationLimit)
    .def("run", &Gillespie::run);

    py::class_<RibosomeSimulator, Gillespie> (mod, "ribosomesimulator")
    .def(py::init<>()) //constructor
    .def("setNumberOfRibosomes", &RibosomeSimulator::setNumberOfRibosomes)
    .def("setCodonForSimulation", &RibosomeSimulator::setCodonForSimulation)
    .def("run_and_get_times", [](RibosomeSimulator &rs) {double d=0.0; double t=0.0; rs.run_and_get_times(d, t); return std::make_tuple(d, t); });
    
}
#endif


RibosomeSimulator::RibosomeSimulator()
{
    // initialize the random generator
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    gen = std::mt19937(rd()); //Standard mersenne_twister_engine seeded with rd()
    dis = std::uniform_real_distribution<>(0, 1);
    //create initial population.
    Eigen::MatrixXi population(32, 1);
    population.fill(0);
    population(0,0) = 1;
    current_population = population;

}

void Simulations::RibosomeSimulator::loadConcentrations(std::string file_name)
{
    csv_utils::ConcentrationsReader cr;
    cr.loadConcentrations(file_name);
    std::vector<csv_utils::concentration_entry> concentrations_vector;
    concentrations_reader = cr;
    // copied code from RibosomeSimulator(csv_utils::concentrations_reader& cr)
    //TODO: NEEDS TO IMPROVE software engineering here.
    std::vector<std::string> stop_codons = {"UAG", "UAA", "UGA"};
    std::vector<csv_utils::concentration_entry> codons_concentrations;
    cr.getContents(codons_concentrations);
    for (csv_utils::concentration_entry entry:codons_concentrations) {
        auto result = std::find(stop_codons.begin(), stop_codons.end(), entry.codon);
        if (result == end(stop_codons)) {
            //Not a stop codon. Proceed.
            reactions_map[entry.codon] = createReactionsGraph(entry);
        }
    }
}

void RibosomeSimulator::setCodonForSimulation(const std::string& codon)
{
    reactions_graph = reactions_map.at(codon);
}

void RibosomeSimulator::run_and_get_times(double& decoding_time, double& translocation_time)
{
    std::vector<double> dt_history;
    dt_history.clear();
    std::vector<int> ribosome_state_history;
    ribosome_state_history.clear();

    // initialize the random generator
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(0, 1);

    double r1 = 0, r2 = 0;
    double tau = 0, clock = 0.0;
    std::vector<double> alphas;
    std::vector<int> next_state;
    while (true)
    {
        // update history
        dt_history.push_back(tau);
        ribosome_state_history.push_back(getState());
        // randomly generate parameter for calculating dt
        r1 = dis(gen) + DBL_MIN; // adding minumum double value in order to avoid division by zero and infinities.
        // randomly generate parameter for selecting reaction
        r2 = dis(gen)+ DBL_MIN; // adding minumum double value in order to avoid division by zero and infinities.
        // calculate an
        getAlphas(alphas, next_state);
        if (alphas.empty())
        {
            translocation_time = 0;
            decoding_time = 0;
            // no available reactions, get times and quit.
            bool is_translocating = true;
            for (int i = ribosome_state_history.size() -1; i >=0; i--){
                if (is_translocating) {
                    translocation_time+= dt_history[i];
                    if (ribosome_state_history[i]<23) is_translocating = false;
                } else {
                    decoding_time += dt_history[i];
                }
            }
            return;
        }
        double a0 = std::accumulate(alphas.begin(), alphas.end(), 0.0);
        // select next reaction to execute
        double cumsum = 0;
        int selected_alpha_vector_index = -1;
        // TODO: vectorization of this loop would increase performance
        do {
            selected_alpha_vector_index++;
            cumsum += alphas[selected_alpha_vector_index]; 
        } while (cumsum < a0 * r2);
        //Apply reaction
        setState(next_state[selected_alpha_vector_index]);
        //Update time
        tau = (1.0/a0) * log(1.0/r1);
        clock += tau;
    }
    
}

void RibosomeSimulator::getDecodingAndTranslocationTimes(double& decoding_time, double& translocation_time)
{
    int translocation_index = 0, i = 0;
    decoding_time = 0; // avoid propagating errors.
    translocation_time = 0; // avoid propagating errors.
    for (Eigen::MatrixXi population:population_history) {
        if (population(24,0) == 1) translocation_index = i;
        i++;
    }
    for (int ii = 0; ii < translocation_index; ii++) {
        decoding_time += dt_history[ii];
    }
    for (unsigned int ii = translocation_index; ii < dt_history.size(); ii++) translocation_time += dt_history[ii];
}


std::vector<std::vector<std::tuple<double, int>>> RibosomeSimulator::createReactionsGraph(const csv_utils::concentration_entry& codon)
{
    double totalconc = 1.9e-4;
    double nonconc = totalconc - codon.wc_cognate_conc - codon.wobblecognate_conc - codon.nearcognate_conc;
    // based on yeast value of 226000 molecules per cell as determined
    // in von der Haar 2008 (PMID 18925958)
    double eEF2conc = 1.36e-5;
    // constants for WCcognate interaction in 1/sec
    double WC1f = 1.4e8*codon.wc_cognate_conc;
    double WC1r = 85;
    double WC2f = 190;
    double WC2r = 0.23;
    double WC3f = 260;
    double WC4f = 1000;
    double WC5f = 1000;
    double WCdiss = 60;
    double WC6f = 1000;
    double WC7f = 200;

    // constants for wobblecognate interaction in 1/sec
    double wobble1f = 1.4e8*codon.wobblecognate_conc;
    double wobble1r = 85;
    double wobble2f = 190;
    double wobble2r = 1;
    double wobble3f = 25;
    double wobble4f = 1000;
    double wobble5f = 1000;
    double wobblediss = 1.1;
    double wobble6f = 1.6;
    double wobble7f = 200;

    // constants for nearcognate interaction in 1/sec
    double near1f = 1.4e8*codon.nearcognate_conc;
    double near1r = 85;
    double near2f = 190;
    double near2r = 80;
    double near3f = 0.4;
    double near4f = 1000;
    double near5f = 1000;
    double neardiss = 1000;
    double near6f = 60;
    double near7f = 200;

    // constants for noncognate interaction in 1/sec.
    // Non-cognates are assumed to not undergo any significant
    // interaction but to simply dissociate quickly.
    double non1f = 1.4e8*nonconc;
    double non1r = 1e5;

    // constants for translocation in 1/sec
    // 150 uM-1 s-1 = is from Fluitt et al 2007 (PMID 17897886)
    double trans1f = eEF2conc*1.5e8;
    double trans1r = 140;
    double trans2 = 250;
    double trans3 = 350;
    double trans4 = 1000;
    double trans5 = 1000;
    double trans6 = 1000;
    double trans7 = 1000;
    double trans8 = 1000;
    double trans9 = 1000;
    
    std::vector<double> ks = {non1f, near1f, wobble1f, WC1f, non1r, near1r, near2f, near2r, near3f, near4f, near5f, neardiss, near6f, near7f, trans1f, wobble1r, wobble2f, wobble2r, wobble3f, wobble4f, wobble5f, wobblediss, wobble6f, wobble7f, trans1f, WC1r, WC2f, WC2r, WC3f, WC4f, WC5f, WCdiss, WC6f, WC7f, trans1f, trans1r, trans2, trans3, trans4, trans5, trans6, trans7, trans8, trans9};
    
    std::vector<std::string> reactions_identifiers = {"non1f", "near1f", "wobble1f", "WC1f", "non1r", "near1r", "near2f", "near2r", "near3f", "near4f", "near5f", "neardiss", "near6f", "near7f", "trans1f", "wobble1r", "wobble2f", "wobble2r", "wobble3f", "wobble4f", "wobble5f", "wobblediss", "wobble6f", "wobble7f", "trans1f", "WC1r", "WC2f", "WC2r", "WC3f", "WC4f", "WC5f", "WCdiss", "WC6f", "WC7f", "trans1f", "trans1r", "trans2", "trans3", "trans4", "trans5", "trans6", "trans7", "trans8", "trans9"};

    Eigen::MatrixXi reactionMatrix[44];
    // build the vector of reactions.
    // [] x=0 -> non1f:(x'=1);
    reactionMatrix[0].resize(32, 1);
    reactionMatrix[0].fill(0);
    reactionMatrix[0](0,0) = -1;
    reactionMatrix[0](1,0) = 1;
    
    // [] x=0 -> near1f:(x'=2);
    reactionMatrix[1].resize(32, 1);
    reactionMatrix[1].fill(0);
    reactionMatrix[1](0,0) = -1;
    reactionMatrix[1](2,0) = 1;
    
    // [] x=0 -> wobble1f:(x'=9);
    reactionMatrix[2].resize(32, 1);
    reactionMatrix[2].fill(0);
    reactionMatrix[2](0,0) = -1;
    reactionMatrix[2](9,0) = 1;

    // [] x=0 -> WC1f:(x'=16);
    reactionMatrix[3].resize(32, 1);
    reactionMatrix[3].fill(0);
    reactionMatrix[3](0,0) = -1;
    reactionMatrix[3](16,0) = 1;

    // [] x=1 -> non1r:(x'=0);
    reactionMatrix[4].resize(32, 1);
    reactionMatrix[4].fill(0);
    reactionMatrix[4](1,0) = -1;
    reactionMatrix[4](0,0) = 1;
    
    // [] x=2 -> near1r:(x'=0);
    reactionMatrix[5].resize(32, 1);
    reactionMatrix[5].fill(0);
    reactionMatrix[5](2,0) = -1;
    reactionMatrix[5](0,0) = 1;
    
    // [] x=2 -> near2f:(x'=3);
    reactionMatrix[6].resize(32, 1);
    reactionMatrix[6].fill(0);
    reactionMatrix[6](2,0) = -1;
    reactionMatrix[6](3,0) = 1;
    
    // [] x=3 -> near2r:(x'=2);
    reactionMatrix[7].resize(32, 1);
    reactionMatrix[7].fill(0);
    reactionMatrix[7](3,0) = -1;
    reactionMatrix[7](2,0) = 1;
    
    // [] x=3 -> near3f:(x'=4);
    reactionMatrix[8].resize(32, 1);
    reactionMatrix[8].fill(0);
    reactionMatrix[8](3,0) = -1;
    reactionMatrix[8](4,0) = 1;
    
    // [] x=4 -> near4f:(x'=5);
    reactionMatrix[9].resize(32, 1);
    reactionMatrix[9].fill(0);
    reactionMatrix[9](4,0) = -1;
    reactionMatrix[9](5,0) = 1;
    
    // [] x=5 -> near5f:(x'=6);
    reactionMatrix[10].resize(32, 1);
    reactionMatrix[10].fill(0);
    reactionMatrix[10](5,0) = -1;
    reactionMatrix[10](6,0) = 1;
    
    // [] x=6 -> neardiss:(x'=0);
    reactionMatrix[11].resize(32, 1);
    reactionMatrix[11].fill(0);
    reactionMatrix[11](6,0) = -1;
    reactionMatrix[11](0,0) = 1;
    
    // [] x=6 -> near6f:(x'=7);
    reactionMatrix[12].resize(32, 1);
    reactionMatrix[12].fill(0);
    reactionMatrix[12](6,0) = -1;
    reactionMatrix[12](7,0) = 1;
    
    // [] x=7 -> near7f:(x'=8);
    reactionMatrix[13].resize(32, 1);
    reactionMatrix[13].fill(0);
    reactionMatrix[13](7,0) = -1;
    reactionMatrix[13](8,0) = 1;
    
    // [] x=8 -> trans1f:(x'=23);
    reactionMatrix[14].resize(32, 1);
    reactionMatrix[14].fill(0);
    reactionMatrix[14](8,0) = -1;
    reactionMatrix[14](23,0) = 1;
    
    // [] x=9 -> wobble1r:(x'=0);
    reactionMatrix[15].resize(32, 1);
    reactionMatrix[15].fill(0);
    reactionMatrix[15](9,0) = -1;
    reactionMatrix[15](0,0) = 1;
    
    // [] x=9 -> wobble2f:(x'=10);
    reactionMatrix[16].resize(32, 1);
    reactionMatrix[16].fill(0);
    reactionMatrix[16](9,0) = -1;
    reactionMatrix[16](10,0) = 1;
    
    // [] x=10 -> wobble2r:(x'=9);
    reactionMatrix[17].resize(32, 1);
    reactionMatrix[17].fill(0);
    reactionMatrix[17](10,0) = -1;
    reactionMatrix[17](9,0) = 1;

    // [] x=10 -> wobble3f:(x'=11);
    reactionMatrix[18].resize(32, 1);
    reactionMatrix[18].fill(0);
    reactionMatrix[18](10,0) = -1;
    reactionMatrix[18](11,0) = 1;

    // [] x=11 -> wobble4f:(x'=12);
    reactionMatrix[19].resize(32, 1);
    reactionMatrix[19].fill(0);
    reactionMatrix[19](11,0) = -1;
    reactionMatrix[19](12,0) = 1;

    // [] x=12 -> wobble5f:(x'=13);
    reactionMatrix[20].resize(32, 1);
    reactionMatrix[20].fill(0);
    reactionMatrix[20](12,0) = -1;
    reactionMatrix[20](13,0) = 1;
    
    // [] x=13 -> wobblediss:(x'=0);
    reactionMatrix[21].resize(32, 1);
    reactionMatrix[21].fill(0);
    reactionMatrix[21](13,0) = -1;
    reactionMatrix[21](0,0) = 1;
    
    // [] x=13 -> wobble6f:(x'=14);
    reactionMatrix[22].resize(32, 1);
    reactionMatrix[22].fill(0);
    reactionMatrix[22](13,0) = -1;
    reactionMatrix[22](14,0) = 1;
    
    // [] x=14 -> wobble7f:(x'=15);
    reactionMatrix[23].resize(32, 1);
    reactionMatrix[23].fill(0);
    reactionMatrix[23](14,0) = -1;
    reactionMatrix[23](15,0) = 1;
    
    // [] x=15 -> trans1f:(x'=23);
    reactionMatrix[24].resize(32, 1);
    reactionMatrix[24].fill(0);
    reactionMatrix[24](15,0) = -1;
    reactionMatrix[24](23,0) = 1;
    
    // [] x=16 -> WC1r:(x'=0);
    reactionMatrix[25].resize(32, 1);
    reactionMatrix[25].fill(0);
    reactionMatrix[25](16,0) = -1;
    reactionMatrix[25](0,0) = 1;
    
    // [] x=16 -> WC2f:(x'=17);
    reactionMatrix[26].resize(32, 1);
    reactionMatrix[26].fill(0);
    reactionMatrix[26](16,0) = -1;
    reactionMatrix[26](17,0) = 1;
    
    // [] x=17 -> WC2r:(x'=16);
    reactionMatrix[27].resize(32, 1);
    reactionMatrix[27].fill(0);
    reactionMatrix[27](17,0) = -1;
    reactionMatrix[27](16,0) = 1;
    
    // [] x=17 -> WC3f:(x'=18);
    reactionMatrix[28].resize(32, 1);
    reactionMatrix[28].fill(0);
    reactionMatrix[28](17,0) = -1;
    reactionMatrix[28](18,0) = 1;
    
    // [] x=18 -> WC4f:(x'=19);
    reactionMatrix[29].resize(32, 1);
    reactionMatrix[29].fill(0);
    reactionMatrix[29](18,0) = -1;
    reactionMatrix[29](19,0) = 1;
    
    // [] x=19 -> WC5f:(x'=20);
    reactionMatrix[30].resize(32, 1);
    reactionMatrix[30].fill(0);
    reactionMatrix[30](19,0) = -1;
    reactionMatrix[30](20,0) = 1;
    
    // [] x=20 -> WCdiss:(x'=0);
    reactionMatrix[31].resize(32, 1);
    reactionMatrix[31].fill(0);
    reactionMatrix[31](20,0) = -1;
    reactionMatrix[31](0,0) = 1;
    
    // [] x=20 -> WC6f:(x'=21);
    reactionMatrix[32].resize(32, 1);
    reactionMatrix[32].fill(0);
    reactionMatrix[32](20,0) = -1;
    reactionMatrix[32](21,0) = 1;
    
    // [] x=21 -> WC7f:(x'=22);
    reactionMatrix[33].resize(32, 1);
    reactionMatrix[33].fill(0);
    reactionMatrix[33](21,0) = -1;
    reactionMatrix[33](22,0) = 1;
    
    // [] x=22 -> trans1f:(x'=23);
    reactionMatrix[34].resize(32, 1);
    reactionMatrix[34].fill(0);
    reactionMatrix[34](22,0) = -1;
    reactionMatrix[34](23,0) = 1;
    
    // [] x=23 -> trans1r:(x'=22);
    reactionMatrix[35].resize(32, 1);
    reactionMatrix[35].fill(0);
    reactionMatrix[35](23,0) = -1;
    reactionMatrix[35](22,0) = 1;
    
    // [] x=23 -> trans2:(x'=24);
    reactionMatrix[36].resize(32, 1);
    reactionMatrix[36].fill(0);
    reactionMatrix[36](23,0) = -1;
    reactionMatrix[36](24,0) = 1;
    
    // [] x=24 -> trans3:(x'=25);
    reactionMatrix[37].resize(32, 1);
    reactionMatrix[37].fill(0);
    reactionMatrix[37](24,0) = -1;
    reactionMatrix[37](25,0) = 1;
    
    // [] x=25 -> trans4:(x'=26);
    reactionMatrix[38].resize(32, 1);
    reactionMatrix[38].fill(0);
    reactionMatrix[38](25,0) = -1;
    reactionMatrix[38](26,0) = 1;
    
    // [] x=26 -> trans5:(x'=27);
    reactionMatrix[39].resize(32, 1);
    reactionMatrix[39].fill(0);
    reactionMatrix[39](26,0) = -1;
    reactionMatrix[39](27,0) = 1;
    
    // [] x=27 -> trans6:(x'=28);
    reactionMatrix[40].resize(32, 1);
    reactionMatrix[40].fill(0);
    reactionMatrix[40](27,0) = -1;
    reactionMatrix[40](28,0) = 1;
    
    // [] x=28 -> trans7:(x'=29);
    reactionMatrix[41].resize(32, 1);
    reactionMatrix[41].fill(0);
    reactionMatrix[41](28,0) = -1;
    reactionMatrix[41](29,0) = 1;
    
    // [] x=29 -> trans8:(x'=30);
    reactionMatrix[42].resize(32, 1);
    reactionMatrix[42].fill(0);
    reactionMatrix[42](29,0) = -1;
    reactionMatrix[42](30,0) = 1;
    
    // [] x=30 -> trans9:(x'=31);
    reactionMatrix[43].resize(32, 1);
    reactionMatrix[43].fill(0);
    reactionMatrix[43](30,0) = -1;
    reactionMatrix[43](31,0) = 1;

    int ii = 0;
    std::vector<std::vector<std::tuple<double, int>>> r_g;
    r_g.resize(32);
    std::fill(r_g.begin(), r_g.end(), std::vector<std::tuple<double, int>>());
    // the vector reactions_graph (I know, not a good name. needs to be changed at some point.), have the following format:
    // reactions_graph[current ribisome state] = [vector of tuples(reaction propensity, ribosome state)]
    // this way, if the ribosome state is, say, 0, we check the possible reactions at reactions_graph[0].
    // if, say we select the reaction with the tuple (0.3, 16), it means that the reaction propensity is 0.3 and
    // it will make the ribosome state go to 16. This is purely for the sake of optimization.
    // the loop below populates reactions_graph automatically. It assumes that each reaction is first-degree.
    for (Eigen::MatrixXi m:reactionMatrix) {
        if (ks.at(ii) > 0){ 
            //populate the local index.
            Eigen::Index maxRow, maxCol, minRow, minCol;
            m.maxCoeff(&maxRow, &maxCol); // 1
            m.minCoeff(&minRow, &minCol); // -1
            r_g.at(minRow).push_back({ks.at(ii), maxRow});
        }
        ii++;
    }
    return r_g;
}

int Simulations::RibosomeSimulator::getState()
{
    Eigen::MatrixXi::Index state;
    current_population.col(0).maxCoeff(&state);
    return state;
}
void Simulations::RibosomeSimulator::setState(int s)
{
    current_population.fill(0);
    current_population(s,0) = 1;
}

void Simulations::RibosomeSimulator::getAlphas(std::vector<double>& as, std::vector<int>& reactions_index)
{
    as.clear();
    reactions_index.clear();
    Eigen::Index state;
    current_population.col(0).maxCoeff(&state); // get the current ribosome state
    auto alphas_and_indexes = reactions_graph[state]; //go the possible reactions of that state.
    double k;
    int index;
    for (auto element:alphas_and_indexes){
        std::tie(k, index) = element;
        as.push_back(k);
        reactions_index.push_back(index);
    }
}
