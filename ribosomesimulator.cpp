#include "ribosomesimulator.h"
#include <float.h>
#include <algorithm>
#include <eigen3/Eigen/Dense>
#include <numeric>
#include "concentrationsreader.h"

#ifdef COMIPLE_PYTHON_MODULE

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

PYBIND11_MODULE(ribosomesimulator, mod) {
  py::class_<Simulations::RibosomeSimulator>(mod, "ribosomesimulator")
      .def(py::init<>())  // constructor
      .def("loadConcentrations",
           &Simulations::RibosomeSimulator::loadConcentrations)
      .def("setCodonForSimulation",
           &Simulations::RibosomeSimulator::setCodonForSimulation)
      .def("setState", &Simulations::RibosomeSimulator::setState)
      .def("run_and_get_times",
           [](Simulations::RibosomeSimulator& rs) {
             double d = 0.0;
             double t = 0.0;
             rs.run_and_get_times(d, t);
             return std::make_tuple(d, t);
           })
      .def("setPropensities",
           &Simulations::RibosomeSimulator::setPropensities)

      .def("getPropensities", &Simulations::RibosomeSimulator::getPropensities)
      .def_readonly("dt_history", &Simulations::RibosomeSimulator::dt_history)
      .def_readonly("ribosome_state_history",
                    &Simulations::RibosomeSimulator::ribosome_state_history);
}
#endif

Simulations::RibosomeSimulator::RibosomeSimulator() {
  // set initial state to 0
  current_state = 0;
}

void Simulations::RibosomeSimulator::loadConcentrations(
    const std::string& file_name) {
  concentrations_reader.loadConcentrations(file_name);
  buildReactionsMap();
}

void Simulations::RibosomeSimulator::buildReactionsMap() {
  std::vector<csv_utils::concentration_entry> codons_concentrations;
  concentrations_reader.getContents(codons_concentrations);
  reactions_map.clear();  // make sure the map is clear.
  for (csv_utils::concentration_entry entry : codons_concentrations) {
    auto result =
        std::find(stop_codons.begin(), stop_codons.end(), entry.codon);
    if (result == end(stop_codons)) {
      // Not a stop codon. Proceed.
      double nonconc = totalconc - entry.wc_cognate_conc -
                       entry.wobblecognate_conc - entry.nearcognate_conc;
      // constants for WCcognate interaction in 1/sec
      WC1f[entry.codon] = 1.4e8 * entry.wc_cognate_conc;

      // constants for wobblecognate interaction in 1/sec
      wobble1f[entry.codon] = 1.4e8 * entry.wobblecognate_conc;

      // constants for nearcognate interaction in 1/sec
      near1f[entry.codon] = 1.4e8 * entry.nearcognate_conc;

      // constants for noncognate interaction in 1/sec.
      // Non-cognates are assumed to not undergo any significant
      // interaction but to simply dissociate quickly.
      non1f[entry.codon] = 1.4e8 * nonconc;

      reactions_map[entry.codon] = createReactionsGraph(entry);
    }
  }
  if (!simulation_codon_3_letters.empty()) {
    reactions_graph = reactions_map.at(simulation_codon_3_letters);
  }
}

void Simulations::RibosomeSimulator::setPropensities(
    std::array<double, 40> prop) {
  for (std::size_t i = 0; i < prop.size(); i++) {
    if (prop.at(i) >= 0) {
      switch (i) {
        case 0:
          non1f[simulation_codon_3_letters] = prop.at(i);
          break;
        case 1:
          near1f[simulation_codon_3_letters] = prop.at(i);
          break;
        case 2:
          wobble1f[simulation_codon_3_letters] = prop.at(i);
          break;
        case 3:
          WC1f[simulation_codon_3_letters] = prop.at(i);
          break;
        case 4:
          non1r = prop.at(i);
          break;
        case 5:
          near1r = prop.at(i);
          break;
        case 6:
          near2f = prop.at(i);
          break;
        case 7:
          near2r = prop.at(i);
          break;
        case 8:
          near3f = prop.at(i);
          break;
        case 9:
          near4f = prop.at(i);
          break;
        case 10:
          near5f = prop.at(i);
          break;
        case 11:
          neardiss = prop.at(i);
          break;
        case 12:
          near6f = prop.at(i);
          break;
        case 13:
          wobble1r = prop.at(i);
          break;
        case 14:
          wobble2f = prop.at(i);
          break;
        case 15:
          wobble2r = prop.at(i);
          break;
        case 16:
          wobble3f = prop.at(i);
          break;
        case 17:
          wobble4f = prop.at(i);
          break;
        case 18:
          wobble5f = prop.at(i);
          break;
        case 19:
          wobblediss = prop.at(i);
          break;
        case 20:
          wobble6f = prop.at(i);
          break;
        case 21:
          WC1r = prop.at(i);
          break;
        case 22:
          WC2f = prop.at(i);
          break;
        case 23:
          WC2r = prop.at(i);
          break;
        case 24:
          WC3f = prop.at(i);
          break;
        case 25:
          WC4f = prop.at(i);
          break;
        case 26:
          WC5f = prop.at(i);
          break;
        case 27:
          WCdiss = prop.at(i);
          break;
        case 28:
          WC6f = prop.at(i);
          break;
        case 29:
          dec7f = prop.at(i);
          break;
        case 30:
          trans1f= prop.at(i);
          break;
        case 31:
          trans1r = prop.at(i);
          break;
        case 32:
          trans2 = prop.at(i);
          break;
        case 33:
          trans3 = prop.at(i);
          break;
        case 34:
          trans4 = prop.at(i);
          break;
        case 35:
          trans5 = prop.at(i);
          break;
        case 36:
          trans6 = prop.at(i);
          break;
        case 37:
          trans7 = prop.at(i);
          break;
        case 38:
          trans8 = prop.at(i);
          break;
        case 39:
          trans9 = prop.at(i);
          break;
      }
    }
  }
}

std::map<std::string, double>
Simulations::RibosomeSimulator::getPropensities() {
  std::map<std::string, double> result;
  std::vector<double> ks = {non1f[simulation_codon_3_letters],
                            near1f[simulation_codon_3_letters],
                            wobble1f[simulation_codon_3_letters],
                            WC1f[simulation_codon_3_letters],
                            non1r,
                            near1r,
                            near2f,
                            near2r,
                            near3f,
                            near4f,
                            near5f,
                            neardiss,
                            near6f,
                            wobble1r,
                            wobble2f,
                            wobble2r,
                            wobble3f,
                            wobble4f,
                            wobble5f,
                            wobblediss,
                            wobble6f,
                            WC1r,
                            WC2f,
                            WC2r,
                            WC3f,
                            WC4f,
                            WC5f,
                            WCdiss,
                            WC6f,
                            dec7f,
                            trans1f,
                            trans1r,
                            trans2,
                            trans3,
                            trans4,
                            trans5,
                            trans6,
                            trans7,
                            trans8,
                            trans9};

//   std::vector<std::string> reactions_identifiers = {
//       "non1f",    "near1f",   "wobble1f", "WC1f",       "non1r",    "near1r",
//       "near2f",   "near2r",   "near3f",   "near4f",     "near5f",   "neardiss",
//       "near6f",   "near7f",   "trans1f",  "wobble1r",   "wobble2f", "wobble2r",
//       "wobble3f", "wobble4f", "wobble5f", "wobblediss", "wobble6f", "wobble7f",
//       "trans1f",  "WC1r",     "WC2f",     "WC2r",       "WC3f",     "WC4f",
//       "WC5f",     "WCdiss",   "WC6f",     "WC7f",       "trans1f",  "trans1r",
//       "trans2",   "trans3",   "trans4",   "trans5",     "trans6",   "trans7",
//       "trans8",   "trans9"};
  for (std::size_t i = 0; i < ks.size(); i++) {
    result[reactions_identifiers[i]] = ks[i];
  }
  return result;
}

void Simulations::RibosomeSimulator::setCodonForSimulation(
    const std::string& codon) {
  simulation_codon_3_letters = codon;
  reactions_graph = reactions_map.at(codon);
}

void Simulations::RibosomeSimulator::run_and_get_times(
    double& decoding_time, double& translocation_time) {
  dt_history.clear();
  ribosome_state_history.clear();

  // initialize the random generator
  std::random_device
      rd;  // Will be used to obtain a seed for the random number engine
  std::mt19937 gen(rd());  // Standard mersenne_twister_engine seeded with rd()
  std::uniform_real_distribution<> dis(0, 1);

  double r1 = 0, r2 = 0;
  double tau = 0, clock = 0.0;
  std::vector<double> alphas;
  std::vector<int> next_state;
  while (true) {
    // update history
    dt_history.push_back(tau);
    ribosome_state_history.push_back(getState());
    // randomly generate parameter for calculating dt
    r1 = dis(gen) + DBL_MIN;  // adding minumum double value in order to avoid
                              // division by zero and infinities.
    // randomly generate parameter for selecting reaction
    r2 = dis(gen) + DBL_MIN;  // adding minumum double value in order to avoid
                              // division by zero and infinities.
    // calculate an
    getAlphas(alphas, next_state);
    if (alphas.empty()) {
      translocation_time = 0;
      decoding_time = 0;
      // no available reactions, get times and quit.
      bool is_translocating = true;
      for (int i = static_cast<int>(ribosome_state_history.size() - 1); i >= 0;
           i--) {
        if (is_translocating) {
          translocation_time += dt_history[static_cast<std::size_t>(i)];
          if (ribosome_state_history[static_cast<std::size_t>(i)] < 23) {
            is_translocating = false;
          }
        } else {
          decoding_time += dt_history[static_cast<std::size_t>(i)];
        }
      }
      return;
    }
    double a0 = std::accumulate(alphas.begin(), alphas.end(), 0.0);
    // select next reaction to execute
    double cumsum = 0;
    int selected_alpha_vector_index = -1;
    // TODO(Heday): vectorization of this loop would increase performance
    do {
      selected_alpha_vector_index++;
      cumsum += alphas[static_cast<std::size_t>(selected_alpha_vector_index)];
    } while (cumsum < a0 * r2);
    // Apply reaction
    setState(next_state[static_cast<std::size_t>(selected_alpha_vector_index)]);
    // Update time
    tau = (1.0 / a0) * log(1.0 / r1);
    clock += tau;
  }
}

std::vector<std::vector<std::tuple<std::reference_wrapper<double>, int>>>
Simulations::RibosomeSimulator::createReactionsGraph(
    const csv_utils::concentration_entry& codon) {
  std::array<std::reference_wrapper<double>, 40> ks = {{non1f[codon.codon],
                                                        near1f[codon.codon],
                                                        wobble1f[codon.codon],
                                                        WC1f[codon.codon],
                                                        non1r,
                                                        near1r,
                                                        near2f,
                                                        near2r,
                                                        near3f,
                                                        near4f,
                                                        near5f,
                                                        neardiss,
                                                        near6f,
                                                        wobble1r,
                                                        wobble2f,
                                                        wobble2r,
                                                        wobble3f,
                                                        wobble4f,
                                                        wobble5f,
                                                        wobblediss,
                                                        wobble6f,
                                                        WC1r,
                                                        WC2f,
                                                        WC2r,
                                                        WC3f,
                                                        WC4f,
                                                        WC5f,
                                                        WCdiss,
                                                        WC6f,
                                                        dec7f,
                                                        trans1f,
                                                        trans1r,
                                                        trans2,
                                                        trans3,
                                                        trans4,
                                                        trans5,
                                                        trans6,
                                                        trans7,
                                                        trans8,
                                                        trans9}};

  Eigen::MatrixXi reactionMatrix[40];
  // build the vector of reactions.
  // [] x=0 -> non1f:(x'=1);
  reactionMatrix[0].resize(32, 1);
  reactionMatrix[0].fill(0);
  reactionMatrix[0](0, 0) = -1;
  reactionMatrix[0](1, 0) = 1;

  // [] x=0 -> near1f:(x'=2);
  reactionMatrix[1].resize(32, 1);
  reactionMatrix[1].fill(0);
  reactionMatrix[1](0, 0) = -1;
  reactionMatrix[1](2, 0) = 1;

  // [] x=0 -> wobble1f:(x'=9);
  reactionMatrix[2].resize(32, 1);
  reactionMatrix[2].fill(0);
  reactionMatrix[2](0, 0) = -1;
  reactionMatrix[2](9, 0) = 1;

  // [] x=0 -> WC1f:(x'=16);
  reactionMatrix[3].resize(32, 1);
  reactionMatrix[3].fill(0);
  reactionMatrix[3](0, 0) = -1;
  reactionMatrix[3](16, 0) = 1;

  // [] x=1 -> non1r:(x'=0);
  reactionMatrix[4].resize(32, 1);
  reactionMatrix[4].fill(0);
  reactionMatrix[4](1, 0) = -1;
  reactionMatrix[4](0, 0) = 1;

  // [] x=2 -> near1r:(x'=0);
  reactionMatrix[5].resize(32, 1);
  reactionMatrix[5].fill(0);
  reactionMatrix[5](2, 0) = -1;
  reactionMatrix[5](0, 0) = 1;

  // [] x=2 -> near2f:(x'=3);
  reactionMatrix[6].resize(32, 1);
  reactionMatrix[6].fill(0);
  reactionMatrix[6](2, 0) = -1;
  reactionMatrix[6](3, 0) = 1;

  // [] x=3 -> near2r:(x'=2);
  reactionMatrix[7].resize(32, 1);
  reactionMatrix[7].fill(0);
  reactionMatrix[7](3, 0) = -1;
  reactionMatrix[7](2, 0) = 1;

  // [] x=3 -> near3f:(x'=4);
  reactionMatrix[8].resize(32, 1);
  reactionMatrix[8].fill(0);
  reactionMatrix[8](3, 0) = -1;
  reactionMatrix[8](4, 0) = 1;

  // [] x=4 -> near4f:(x'=5);
  reactionMatrix[9].resize(32, 1);
  reactionMatrix[9].fill(0);
  reactionMatrix[9](4, 0) = -1;
  reactionMatrix[9](5, 0) = 1;

  // [] x=5 -> near5f:(x'=6);
  reactionMatrix[10].resize(32, 1);
  reactionMatrix[10].fill(0);
  reactionMatrix[10](5, 0) = -1;
  reactionMatrix[10](6, 0) = 1;

  // [] x=6 -> neardiss:(x'=0);
  reactionMatrix[11].resize(32, 1);
  reactionMatrix[11].fill(0);
  reactionMatrix[11](6, 0) = -1;
  reactionMatrix[11](0, 0) = 1;

//   // [] x=6 -> near6f:(x'=7);
//   reactionMatrix[12].resize(32, 1);
//   reactionMatrix[12].fill(0);
//   reactionMatrix[12](6, 0) = -1;
//   reactionMatrix[12](7, 0) = 1;

  //   // [] x=6 -> dec7f:(x'=21);
  reactionMatrix[12].resize(32, 1);
  reactionMatrix[12].fill(0);
  reactionMatrix[12](6, 0) = -1;
  reactionMatrix[12](21, 0) = 1;

  // [] x=7 -> near7f:(x'=8);
//   reactionMatrix[13].resize(32, 1);
//   reactionMatrix[13].fill(0);
//   reactionMatrix[13](7, 0) = -1;
//   reactionMatrix[13](8, 0) = 1;

  // [] x=8 -> trans1f:(x'=23);
//   reactionMatrix[14].resize(32, 1);
//   reactionMatrix[14].fill(0);
//   reactionMatrix[14](8, 0) = -1;
//   reactionMatrix[14](23, 0) = 1;

  // [] x=9 -> wobble1r:(x'=0);
  reactionMatrix[13].resize(32, 1);
  reactionMatrix[13].fill(0);
  reactionMatrix[13](9, 0) = -1;
  reactionMatrix[13](0, 0) = 1;

  // [] x=9 -> wobble2f:(x'=10);
  reactionMatrix[14].resize(32, 1);
  reactionMatrix[14].fill(0);
  reactionMatrix[14](9, 0) = -1;
  reactionMatrix[14](10, 0) = 1;

  // [] x=10 -> wobble2r:(x'=9);
  reactionMatrix[15].resize(32, 1);
  reactionMatrix[15].fill(0);
  reactionMatrix[15](10, 0) = -1;
  reactionMatrix[15](9, 0) = 1;

  // [] x=10 -> wobble3f:(x'=11);
  reactionMatrix[16].resize(32, 1);
  reactionMatrix[16].fill(0);
  reactionMatrix[16](10, 0) = -1;
  reactionMatrix[16](11, 0) = 1;

  // [] x=11 -> wobble4f:(x'=12);
  reactionMatrix[17].resize(32, 1);
  reactionMatrix[17].fill(0);
  reactionMatrix[17](11, 0) = -1;
  reactionMatrix[17](12, 0) = 1;

  // [] x=12 -> wobble5f:(x'=13);
  reactionMatrix[18].resize(32, 1);
  reactionMatrix[18].fill(0);
  reactionMatrix[18](12, 0) = -1;
  reactionMatrix[18](13, 0) = 1;

  // [] x=13 -> wobblediss:(x'=0);
  reactionMatrix[19].resize(32, 1);
  reactionMatrix[19].fill(0);
  reactionMatrix[19](13, 0) = -1;
  reactionMatrix[19](0, 0) = 1;

//   // [] x=13 -> wobble6f:(x'=14);
//   reactionMatrix[22].resize(32, 1);
//   reactionMatrix[22].fill(0);
//   reactionMatrix[22](13, 0) = -1;
//   reactionMatrix[22](14, 0) = 1;

// [] x=13 -> dec7f:(x'=21);
  reactionMatrix[20].resize(32, 1);
  reactionMatrix[20].fill(0);
  reactionMatrix[20](13, 0) = -1;
  reactionMatrix[20](21, 0) = 1;

//   // [] x=14 -> wobble7f:(x'=15);
//   reactionMatrix[23].resize(32, 1);
//   reactionMatrix[23].fill(0);
//   reactionMatrix[23](14, 0) = -1;
//   reactionMatrix[23](15, 0) = 1;

//   // [] x=15 -> trans1f:(x'=23);
//   reactionMatrix[24].resize(32, 1);
//   reactionMatrix[24].fill(0);
//   reactionMatrix[24](15, 0) = -1;
//   reactionMatrix[24](23, 0) = 1;

  // [] x=16 -> WC1r:(x'=0);
  reactionMatrix[21].resize(32, 1);
  reactionMatrix[21].fill(0);
  reactionMatrix[21](16, 0) = -1;
  reactionMatrix[21](0, 0) = 1;

  // [] x=16 -> WC2f:(x'=17);
  reactionMatrix[22].resize(32, 1);
  reactionMatrix[22].fill(0);
  reactionMatrix[22](16, 0) = -1;
  reactionMatrix[22](17, 0) = 1;

  // [] x=17 -> WC2r:(x'=16);
  reactionMatrix[23].resize(32, 1);
  reactionMatrix[23].fill(0);
  reactionMatrix[23](17, 0) = -1;
  reactionMatrix[23](16, 0) = 1;

  // [] x=17 -> WC3f:(x'=18);
  reactionMatrix[24].resize(32, 1);
  reactionMatrix[24].fill(0);
  reactionMatrix[24](17, 0) = -1;
  reactionMatrix[24](18, 0) = 1;

  // [] x=18 -> WC4f:(x'=19);
  reactionMatrix[25].resize(32, 1);
  reactionMatrix[25].fill(0);
  reactionMatrix[25](18, 0) = -1;
  reactionMatrix[25](19, 0) = 1;

  // [] x=19 -> WC5f:(x'=20);
  reactionMatrix[26].resize(32, 1);
  reactionMatrix[26].fill(0);
  reactionMatrix[26](19, 0) = -1;
  reactionMatrix[26](20, 0) = 1;

  // [] x=20 -> WCdiss:(x'=0);
  reactionMatrix[27].resize(32, 1);
  reactionMatrix[27].fill(0);
  reactionMatrix[27](20, 0) = -1;
  reactionMatrix[27](0, 0) = 1;

  // [] x=20 -> WC6f:(x'=21);
  reactionMatrix[28].resize(32, 1);
  reactionMatrix[28].fill(0);
  reactionMatrix[28](20, 0) = -1;
  reactionMatrix[28](21, 0) = 1;

  // [] x=21 -> WC7f:(x'=22);
  reactionMatrix[29].resize(32, 1);
  reactionMatrix[29].fill(0);
  reactionMatrix[29](21, 0) = -1;
  reactionMatrix[29](22, 0) = 1;

  // [] x=22 -> trans1f:(x'=23);
  reactionMatrix[30].resize(32, 1);
  reactionMatrix[30].fill(0);
  reactionMatrix[30](22, 0) = -1;
  reactionMatrix[30](23, 0) = 1;

  // [] x=23 -> trans1r:(x'=22);
  reactionMatrix[31].resize(32, 1);
  reactionMatrix[31].fill(0);
  reactionMatrix[31](23, 0) = -1;
  reactionMatrix[31](22, 0) = 1;

  // [] x=23 -> trans2:(x'=24);
  reactionMatrix[32].resize(32, 1);
  reactionMatrix[32].fill(0);
  reactionMatrix[32](23, 0) = -1;
  reactionMatrix[32](24, 0) = 1;

  // [] x=24 -> trans3:(x'=25);
  reactionMatrix[33].resize(32, 1);
  reactionMatrix[33].fill(0);
  reactionMatrix[33](24, 0) = -1;
  reactionMatrix[33](25, 0) = 1;

  // [] x=25 -> trans4:(x'=26);
  reactionMatrix[34].resize(32, 1);
  reactionMatrix[34].fill(0);
  reactionMatrix[34](25, 0) = -1;
  reactionMatrix[34](26, 0) = 1;

  // [] x=26 -> trans5:(x'=27);
  reactionMatrix[35].resize(32, 1);
  reactionMatrix[35].fill(0);
  reactionMatrix[35](26, 0) = -1;
  reactionMatrix[35](27, 0) = 1;

  // [] x=27 -> trans6:(x'=28);
  reactionMatrix[36].resize(32, 1);
  reactionMatrix[36].fill(0);
  reactionMatrix[36](27, 0) = -1;
  reactionMatrix[36](28, 0) = 1;

  // [] x=28 -> trans7:(x'=29);
  reactionMatrix[37].resize(32, 1);
  reactionMatrix[37].fill(0);
  reactionMatrix[37](28, 0) = -1;
  reactionMatrix[37](29, 0) = 1;

  // [] x=29 -> trans8:(x'=30);
  reactionMatrix[38].resize(32, 1);
  reactionMatrix[38].fill(0);
  reactionMatrix[38](29, 0) = -1;
  reactionMatrix[38](30, 0) = 1;

  // [] x=30 -> trans9:(x'=31);
  reactionMatrix[39].resize(32, 1);
  reactionMatrix[39].fill(0);
  reactionMatrix[39](30, 0) = -1;
  reactionMatrix[39](31, 0) = 1;

  int ii = 0;
  std::vector<std::vector<std::tuple<std::reference_wrapper<double>, int>>> r_g;
  r_g.resize(32);
  std::fill(r_g.begin(), r_g.end(),
            std::vector<std::tuple<std::reference_wrapper<double>, int>>());
  // the vector reactions_graph (I know, not a good name. needs to be changed at
  // some point.), have the following format: reactions_graph[current ribisome
  // state] = [vector of tuples(reaction propensity, ribosome state)] this way,
  // if the ribosome state is, say, 0, we check the possible reactions at
  // reactions_graph[0]. if, say we select the reaction with the tuple (0.3,
  // 16), it means that the reaction propensity is 0.3 and it will make the
  // ribosome state go to 16. This is purely for the sake of optimization. the
  // loop below populates reactions_graph automatically. It assumes that each
  // reaction is first-degree.
  for (const Eigen::MatrixXi& m : reactionMatrix) {
    if (ks.at(static_cast<std::size_t>(ii)) > 0) {
      // populate the local index.
      Eigen::Index maxRow, maxCol, minRow, minCol;
      m.maxCoeff(&maxRow, &maxCol);  // 1
      m.minCoeff(&minRow, &minCol);  // -1
      r_g.at(static_cast<std::size_t>(minRow))
          .push_back({ks.at(static_cast<std::size_t>(ii)), maxRow});
    }
    ii++;
  }
  return r_g;
}

int Simulations::RibosomeSimulator::getState() { return current_state; }
void Simulations::RibosomeSimulator::setState(int s) { current_state = s; }

void Simulations::RibosomeSimulator::getAlphas(
    std::vector<double>& as, std::vector<int>& reactions_index) {
  as.clear();
  reactions_index.clear();
  auto alphas_and_indexes = reactions_graph[static_cast<std::size_t>(
      current_state)];  // go the possible
                        // reactions of that
                        // state.
  double k;
  int index;
  for (auto element : alphas_and_indexes) {
    std::tie(k, index) = element;
    as.push_back(k);
    reactions_index.push_back(index);
  }
}

void Simulations::RibosomeSimulator::getDecodingAlphas(
    std::vector<double>& as, std::vector<int>& reactions_index) {
  as.clear();
  reactions_index.clear();
  auto alphas_and_indexes = reactions_graph[static_cast<std::size_t>(
      current_state)];  // go the possible
                        // reactions of that
                        // state.
  double k;
  int index;
  for (auto element : alphas_and_indexes) {
    std::tie(k, index) = element;
    if (index < 23) {
      as.push_back(k);
      reactions_index.push_back(index);
    }
  }
}
