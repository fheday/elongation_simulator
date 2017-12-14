#ifndef SIMULATIONS_GILLESPIE_H
#define SIMULATIONS_GILLESPIE_H

#include <eigen3/Eigen/Dense>
#include <list>
#include <memory>
#include "reactionsset.h"

namespace Simulations {

class Gillespie {
 public:
  Gillespie();
  Gillespie(int, const Eigen::MatrixXi&, const ReactionsSet&);
  virtual void run();
  void enableHistory(bool);
  void setReactionsSet(const ReactionsSet&);
  void setInitialPopulation(const Eigen::MatrixXi&);
  void setIterationLimit(int);
  void setTimeLimit(double);
  virtual double getReactionTime(double&, double&, std::string&);
  virtual ~Gillespie();
  int iteration_limit = -1;
  double time_limit = -1;
  Eigen::MatrixXi initial_population;
  Eigen::MatrixXi current_population;
  ReactionsSet reactions;
  std::vector<double> dt_history;
  std::vector<Eigen::MatrixXi> population_history;
  double simulation_total_time;
};
}  // namespace Simulations

#endif  // SIMULATIONS_GILLESPIE_H
