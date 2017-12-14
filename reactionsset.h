#ifndef REACTIONSSET_H
#define REACTIONSSET_H

#include <eigen3/Eigen/Dense>
#include <tuple>
#include <vector>

namespace Simulations {

class ReactionsSet {
 public:
  ReactionsSet();
  void addReaction(Eigen::MatrixXi, double);
  void addReaction(Eigen::MatrixXi, double, std::string);
  Eigen::MatrixXi getReaction(int);
  void getAlphas(const Eigen::MatrixXi&, std::vector<double>&,
                 std::vector<int>&);
  ~ReactionsSet();
  std::vector<Eigen::MatrixXi> reactions_vector;
  std::vector<double> ks;
  std::vector<std::string> descriptions;

 private:
  // the elements of this vector are: k (reaction propensity),
  // vector of indexes (col, row) of the reactants.
  std::vector<std::tuple<double, std::vector<std::tuple<int, int>>>>
      k_pop_index;
};
}  // namespace Simulations
#endif  // REACTIONSSET_H
