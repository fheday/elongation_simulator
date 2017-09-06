#ifndef MRNAELEMENT_H
#define MRNAELEMENT_H

#include <eigen3/Eigen/Dense>

namespace Simulations{
    
    class mRNAElement
    {
    public:
        virtual ~mRNAElement() {};
        bool isAvailable = true; // true if the position can be used.
        bool isOccupied = false; // true if there is  ribosome in the position. As the ribosome moves, it sets the next 'isAvailable' to false, and the 10th previous 'isAvailable' to true. When terminates, sets last 10 'isAvailable' to true.
        std::string codon;
        int index;
        virtual void getAlphas(Eigen::VectorXd& as, Eigen::VectorXi& reactions_index) {};
        virtual void executeReaction(int) {};
        virtual int getState() {return -1;};
        virtual void setState(int) {};
    };
}


#endif // MRNAELEMENT_H
