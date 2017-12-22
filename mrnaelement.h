#ifndef MRNAELEMENT_H
#define MRNAELEMENT_H

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
  virtual void getAlphas(std::vector<double> &, std::vector<int> &) {}
  virtual void executeReaction(int) {}
  virtual int getState() { return -1; }
  virtual void setState(int) {}
  void setNextCodon(mRNAElement *);
  void setPreviousCodon(mRNAElement *);
  virtual void updateAlphas() {}

 protected:
  std::vector<double> alphas;
  std::vector<int> reactions_index;
  bool is_available = true;  // true if the position can be used.
  bool is_occupied = false;  // true if there is  ribosome in the position. As
                             // the ribosome moves, it sets the next
                             // 'isAvailable' to false, and the 10th previous
                             // 'isAvailable' to true. When terminates, sets
                             // last 10 'isAvailable' to true.
  mRNAElement *next_mRNA_element, *previous_mRNA_element;
};
}  // namespace Simulations

#endif  // MRNAELEMENT_H
