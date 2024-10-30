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


#include <cstddef>
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
  [[nodiscard]] bool isAvailable() const;
  [[nodiscard]] bool isOccupied() const ;
  virtual void executeReaction(int) {}
  virtual int getState() { return -1; }
  virtual void setState(int) {}
  void setNextCodon(mRNAElement *);
  void setPreviousCodon(mRNAElement *);
  virtual void updateAlphas() {}
  void addReactionToHistory(std::size_t state, float dt);
  std::pair<std::vector<std::size_t>, std::vector<float>> getHistory();

  virtual void setPropensities(std::map<std::string, float>&) {}
  virtual void setNoNonCognate(bool) {}

  virtual std::map<std::string, float> getPropensities() {
    return {};
  }

  std::vector<float> alphas;
  std::vector<int> reactions_index;

 protected:
  bool is_available = true;  // true if the position can be used.
  bool is_occupied = false;  // true if there is  ribosome in the position. As
                             // the ribosome moves, it sets the next
                             // 'isAvailable' to false, and the 10th previous
                             // 'isAvailable' to true. When terminates, sets
                             // last 10 'isAvailable' to true.
  mRNAElement *next_mRNA_element, *previous_mRNA_element;
  std::vector<std::size_t>
      state_history;  // here we store the history of the codon states
  std::vector<float>
      dt_history;  // here we store the history of reactions times.
};
}  // namespace Simulations

#endif  // MRNAELEMENT_H
