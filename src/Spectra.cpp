#include "Spectra.h"

Spectra::Spectra(const double tval): t(tval) {}

Spectra::Spectra(const double &tval, const std::vector<Reggeon> &regs): t(tval), reggeons(regs){}

Spectra::Spectra(const Spectra &spec):
    t(spec.t), reggeons(spec.reggeons){}

void Spectra::setT(const double tval)
{
    // Setter of ts
    t = tval ;
}

void Spectra::setReggeons(const std::vector<Reggeon> &regs)
{
    // Setter of kernels
    reggeons = regs ;
}

double Spectra::getT() const
{
    // Getter of ts
    return t;
}

std::vector<Reggeon> Spectra::getReggeons() const
{
    // Getter of Kernels
    return reggeons;
}

void Spectra::copy(const Spectra& rhs)
{
    //Copies Reggeon data
    t = rhs.t;
    reggeons = rhs.reggeons;
}

Spectra& Spectra::operator= (const Spectra& rhs)
{
    if (this == &rhs) return *this;
    copy(rhs);
    return *this;
}

bool Spectra::operator< (const Spectra &rhs)
{
    return this->t < rhs.t ;
}

bool Spectra::operator> (const Spectra &rhs)
{
    return this->t > rhs.t ;
}

Spectra::~Spectra(){}