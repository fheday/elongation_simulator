#ifndef MRNAELEMENT_H
#define MRNAELEMENT_H

/*
 * @file  mrnaelement.h
 * 
 * @brief class to represent codons in a mRNA
 *
 * @author Fabio Hedayioglu
 * Contact: fheday@gmail.com
 *
 */


#include <array>
#include <map>
#include <string>
#include <vector>

namespace Simulations {

class mRNAElement {
 public:
  std::string codon;
  int index;
  mRNAElement();
  virtual ~mRNAElement();
  void setAvailable(bool);
  void setOccupied(bool);
  bool isAvailable();
  bool isOccupied();
  virtual void executeReaction(int) {}
  virtual int getState() { return -1; }
  virtual void setState(int) {}
  void setNextCodon(mRNAElement *);
  void setPreviousCodon(mRNAElement *);
  virtual void updateAlphas() {}
  void addReactionToHistory(int state, double dt);
  std::pair<std::vector<int>, std::vector<double>> getHistory();

  virtual void setPropensities(std::map<std::string, double>) {}
  virtual void setNoNonCognate(bool) {}

  virtual std::map<std::string, double> getPropensities() {
    return std::map<std::string, double>();
  }

  std::vector<double> alphas;
  std::vector<int> reactions_index;

 protected:
  bool is_available = true;  // true if the position can be used.
  bool is_occupied = false;  // true if there is  ribosome in the position. As
                             // the ribosome moves, it sets the next
                             // 'isAvailable' to false, and the 10th previous
                             // 'isAvailable' to true. When terminates, sets
                             // last 10 'isAvailable' to true.
  mRNAElement *next_mRNA_element, *previous_mRNA_element;
  std::vector<int>
      state_history;  // here we store the history of the codon states
  std::vector<double>
      dt_history;  // here we store the history of reactions times.
};
}  // namespace Simulations

#endif  // MRNAELEMENT_H
