#ifndef REGGEON_H
#define REGGEON_H

#include <string>
#include <vector>

class Reggeon
{
    private :
        double J, dJdt ;
        std::vector< std::vector<double> > wf ;
        std::string name ;
        int index ;
        void copy(const Reggeon &rhs);
    public :
        // Default Reggeon constructor
        Reggeon();
        // Reggeon constructor
        Reggeon(const double j, const double djdt, const std::vector< std::vector<double> > &wwf, 
                const std::string n, const int i);
        // Copy constructor
        Reggeon(const Reggeon &reg);
        double getJ() const ;                                           // Getter for js
        double getdJdt() const;                                         // Getter for dJdt
        std::vector< std::vector<double> > getWf() const;               // Getter for wf
        std::string getName() const ;                                   // Getter for name
        int getIndex() const ;                                          // Getter for index
        void setJ(const double j);                                      // Setter for js
        void setdJdt(const double djdt);                                // Setter for dJdt
        void setWf(const std::vector< std::vector<double> > &wwf);      // Setter for wf
        void setName(const std::string n);                              // Setter for name
        void setIndex(const int i) ;                                    // Setter for index
        Reggeon& operator= (const Reggeon &rhs);                        // Assignment operator
        ~Reggeon();                                                     // Destructor

};
#endif