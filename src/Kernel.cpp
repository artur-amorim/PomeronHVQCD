#include <functional>
#include <thread>
#include "HolographicVQCD.h"
#include "Kernel.h"
#include "Reggeon.h"
#include "schrodinger/schrodinger.h"
#include "methods/vectorOperators.hpp"
#include "methods/interpolation/Spline_Interp.hpp"

void Kernel::copy(const Kernel &rhs)
{
    //Copies Kernel data
    name = rhs.name;
    nTrajectories = rhs.nTrajectories;
    trajectories = rhs.trajectories;
    kernelPars = rhs.kernelPars;
}

Kernel::Kernel():
    name(""), nTrajectories(0),
    trajectories({}), kernelPars({})
{}


Kernel::Kernel(const Kernel &kernel):
    name(kernel.name), nTrajectories(kernel.nTrajectories),
    trajectories(kernel.trajectories), kernelPars(kernel.kernelPars)
{}

Kernel::Kernel(const std::string nname, const int ntraj, const std::vector<double> &pars):
    name(nname), nTrajectories(ntraj),
    trajectories({}), kernelPars(pars)
    {
        trajectories.resize(ntraj);
    }

std::string Kernel::Name() const
{
    // Returns name of the Kernel
    return this->name ;
}

int Kernel::NTrajectories() const
{
    // Returns the number of Regge trajectories
    return this->nTrajectories;
}

Spline_Interp<double> Kernel::Trajectory(const int n) const
{
    return trajectories[n];
}

std::vector<double> Kernel::KernelPars() const
{
    // Returns the reggeons in the Kernel
    return kernelPars ;
}

void Kernel::computeReggeTrajectories(const std::vector<double> &pars)
{
    // Update the parameters of the kernel
    if(pars.size() != 0) setKernelPars(pars);
    // Define the js and ts vectors
    const double jmax = 2.1;
    const double h = 0.025;
    int n_js = jmax / h ;
    std::vector<double> js(n_js, 0.0);
    std::vector<std::vector<double> > ts(nTrajectories,std::vector<double>(n_js,0.0));
    // Get the number of threads and number of js to compute per thread
    int nthreads = std::thread::hardware_concurrency();
    int n_js_thread = n_js / nthreads;
    int r = n_js % nthreads; // Reminder of n_js divided by nthreads
    // Function that defines each thread
    std::function<void(int)> f = [&js, &ts, this, h, n_js_thread] (const int th_id)
    {
        unsigned int imin = th_id * n_js_thread;
        unsigned int imax = (th_id+1) * n_js_thread;
        for(unsigned int i = imin; i < imax; i++)
        {
            const double j = i * h;
            js[i] = j;
            // Compute the potential
            std::vector<double> VSch = this->computePotential(j);
            // Compute all the tn(J)s
            List  spectrum = computeSpectrum(ihqcd().getZ(), VSch, this->NTrajectories(), "cheb");
            for(int k = 0; k < this->NTrajectories(); k++) ts[k][i] = spectrum.Es[k];
        }
        return;
    };
    std::thread * ths = new std::thread[nthreads];
    for(int i = 0; i < nthreads; i++) ths[i] = std::thread(f, i);
    for(int i = 0; i < nthreads; i++) ths[i].join();
    delete[] ths;
    f = [&js, &ts, this, h, n_js_thread, nthreads] (const int th_id)
    {
        unsigned int i = n_js_thread * nthreads + th_id;
        const double j = i * h ;
        js[i] = j;
        // Compute the potential
        std::vector<double> VSch = this->computePotential(j);
        // Compute all the tn(J)s
        List  spectrum = computeSpectrum(ihqcd().getZ(), VSch, this->NTrajectories(), "cheb");
        for(int k = 0; k < this->NTrajectories(); k++) ts[k][i] = spectrum.Es[k]; 
        return;
    };
    ths = new std::thread[r];
    for(int i = 0; i < r; i++) ths[i] = std::thread(f, i);
    for(int i = 0; i < r; i++) ths[i].join();
    delete[] ths;
    // Now creat the Spline_Interp objects
    for(int i = 0; i < nTrajectories; i++)
    {
        Spline_Interp<double> trajectory(ts[i], js);
        trajectories[i] = trajectory;
    }
}

void Kernel::setKernelPars(const std::vector<double> &pars)
{    
    // Sets the kernel parameters in the Kernel
    this->kernelPars = pars ;
}

std::vector<double> Kernel::computePotential(const double J) const
{
    std::cout << "Inside computePotential from Kernel class." << std::endl;
    std::vector<double> VSch(0);
    return VSch;
};

Kernel& Kernel::operator= (const Kernel &rhs)
{
    if (this == &rhs) return *this;
    copy(rhs);
    return *this;
};

Kernel::~Kernel(){}


std::vector<Reggeon> computeReggeons(const Kernel &kernel, const double t, const int n)
{
    /* 
        Finds the n Reggeons corresponding to a given t
        First we find J such that tn(J) = t 
    */

   std::vector<Reggeon> reggeons(n);
   for(int i = 0; i < n; i++)
   {
        // Get Regge trajectory
        Spline_Interp<double> reggeTraj = kernel.Trajectory(i);
        // Find jn
        const double jn = reggeTraj.interp(t);
        // Find djn/dt
        const double dJdt = reggeTraj.der1(t);
        std::vector<double> VSch = kernel.computePotential(jn);
        // Compute the wavefunction
        const List spectrum = computeSpectrum(ihqcd().getZ(), VSch, i+1);
        std::vector< std::vector<double> > wf = spectrum.wfs[i];
        // We need to normalize the wavefunction
        Spline_Interp<double> wf2(wf[0], wf[1] * wf[1], 0, 0);
        double c = wf2.integrate(wf[0].back());
        wf[1] = (1/std::sqrt(c)) * (wf[1][2] / std::fabs(wf[1][2])) * wf[1] ;
        Reggeon reggeon(jn, dJdt, wf, kernel.Name() + "." + std::to_string(i+1), i+1);
        reggeons[i] = reggeon;
   }
   return reggeons ;
}