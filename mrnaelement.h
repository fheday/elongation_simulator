#ifndef MRNAELEMENT_H
#define MRNAELEMENT_H

#include <vector>
#include <string>

namespace Simulations{
    
    class mRNAElement
    {
    public:
        virtual ~mRNAElement() {};
        bool isAvailable = true; // true if the position can be used.
        bool isOccupied = false; // true if there is  ribosome in the position. As the ribosome moves, it sets the next 'isAvailable' to false, and the 10th previous 'isAvailable' to true. When terminates, sets last 10 'isAvailable' to true.
        virtual void updateAlphas(bool = true) {};
        std::string codon;
        int index;
        virtual void getAlphas(std::vector<double>&, std::vector<int>&) {};
        virtual void executeReaction(int) {};
        virtual int getState() {return -1;};
        virtual void setState(int) {};
        std::vector<double> alphas;
        std::vector<int> reactions_index;
    };
}


#endif // MRNAELEMENT_H
