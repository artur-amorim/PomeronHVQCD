#include "Reggeon.h"

// Default Reggeon constructor
Reggeon::Reggeon():
    J(0), dJdt(0),
    wf({}), name(""),
    index(0){}

// Reggeon constructor
Reggeon::Reggeon(const double j, const double djdt, const std::vector< std::vector<double> > &wwf, 
                const std::string n, const int i):
                J(j), dJdt(djdt), wf(wwf), name(n), index(i) {}

// Copy Reggeon constructor
Reggeon::Reggeon(const Reggeon &reg):
    J(reg.J), dJdt(reg.dJdt), wf(reg.wf),
    name(reg.name), index(reg.index) {}

double Reggeon::getJ() const
{
    // Getter for J
    return this->J ;
}

double Reggeon::getdJdt() const
{
    // Getter for dJdt
    return this->dJdt ;
}

std::vector< std::vector<double> > Reggeon::getWf() const
{
    // Getter for wf
    return this->wf ;
}

std::string Reggeon::getName() const
{
    // Getter for name
    return this->name ;
}

int Reggeon::getIndex() const
{
    // Getter for index
    return this->index ;
}

void Reggeon::setJ(const double j)
{
    // Setter for J
    this->J = j ;
}

void Reggeon::setdJdt(const double djdt)
{
    // Setter for dJdt
    this->dJdt = djdt ;
}

void Reggeon::setWf(const std::vector< std::vector<double> > &wwf)
{
    // Getter for wf
    this->wf = wwf ;
}

void Reggeon::setName(const std::string n)
{
    // Setter for name
    this->name = n;
}

void Reggeon::setIndex(const int i)
{
    // Setter for index
    this->index = i ;
}

void Reggeon::copy(const Reggeon &rhs)
{
    //Copies Reggeon data
    J = rhs.getJ();
    dJdt = rhs.getdJdt() ;
    wf = rhs.getWf();
    name = rhs.getName();
    index = rhs.getIndex();
}

Reggeon& Reggeon::operator= (const Reggeon &rhs)
{
    if (this == &rhs) return *this;
    copy(rhs);
    return *this;
}

Reggeon::~Reggeon(){}