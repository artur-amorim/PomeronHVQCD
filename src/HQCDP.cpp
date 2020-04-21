#include <algorithm>
#include <stdexcept>
#include "HQCDP.h"
#include "methods/optimization/NelderMead.hpp"

HQCDP::HQCDP(const bool rrsslog, const double spconst):
            rsslog(rrsslog), SPconstraint(spconst),
            processes({}), useTVals({}),
            kernels({}), spectrum({}),
            gns({})
            {}

void HQCDP::setGNs(const std::vector<double> &gs)
{
    gns = gs;
}

std::vector<Kernel * > HQCDP::getKernels() const
{
    // Get Kernels
    return kernels ;
}

std::vector<Spectra> HQCDP::getSpectrum() const
{
    return spectrum;
}

std::vector<double> HQCDP::getUseTVals() const
{
    // Get the tvals to be used
    return useTVals ;
}

void HQCDP::computeNeededTVals()
{
    const int nprocs = processes.size() ;
    // Iterate over all processes
    for(int i = 0; i < nprocs ; i++)
    {
        std::vector<double> proc_ts = processes[i]->getNeededTVals();
        if(proc_ts.size() == 0)
        {
            std::cout << "Process with no t values has been found" << std::endl;
        }
        // Get all the tvals in proc_ts
        for(int j = 0 ; j < proc_ts.size(); j++)
        {
            useTVals.push_back(proc_ts[j]) ;
        }
    }
    // Now we std::sort the values of t and remove repeated values
    std::sort(useTVals.begin(), useTVals.end()) ;
    // Erase any duplicates
    useTVals.erase(std::unique(useTVals.begin(), useTVals.end() ), useTVals.end() );
}

void HQCDP::addProcessObservable(ProcessObservable &proc)
{
    // Add a process to be fitted or predicted
    // For consistency change NMC, alpha and rsslog data members of proc
    // to be the same as the HQCDP object.
    ProcessObservable * newProcess = &proc;
    if(proc.getRSSLOG() != rsslog) newProcess->setRSSLOG(rsslog);
    // Append the process observable proc
    processes.push_back(newProcess) ;
}

void HQCDP::addKernel(Kernel &f)
{
    // Add a kernel that will be used in later fits or predictions
    Kernel * newKernel = &f;
    kernels.push_back(newKernel) ;
}

void HQCDP::computeSpectrum(const std::vector< std::vector<double> > &kernelPars)
{
    /*
        Compute the spectrum given kernelPars. kernelPars is a std::vector containing vectors with the parameters
        of the corresponding kernels presentend in the kernels attribute of the HQCDP object.
        The program starts by checking that the number of such vectors is the same as the number of kernels.
    */
   // If the number of vectors is the same update the parameters of each kernel.
   if (kernels.size() == 0) throw std::runtime_error("Add one kernel first");
   if(kernelPars.size() != 0)
   {
       for(int ker_idx = 0; ker_idx < kernels.size(); ker_idx+=1)
       {
           std::vector<double> pars = kernelPars[ker_idx];
           kernels[ker_idx]->setKernelPars(pars);
       }
   }
   // Compute all the needed tvals if not computed
   if(useTVals.size() == 0) computeNeededTVals();
   // for each value of t we are going to compute the Reggeons of each Kernel object in kernels
   // Clean previous spectrum
   spectrum.clear();
   for(int i = 0; i < useTVals.size(); i++)
   {
       const double t = useTVals[i];
       std::cout << "Computing Reggeons for t = " << t << std::endl;
       std::vector<Reggeon> reggeons;
       // Given t, compute the Reggeons for each kernel
       for(int j = 0; j < kernels.size(); j++)
       {
            Kernel * ker = kernels[j];
            const int n_regs = ker->NTrajectories();
            ker->computeReggeTrajectories();
            std::vector<Reggeon> ker_reggeons = computeReggeons(*ker, t, n_regs);
            for(int k = 0; k < ker_reggeons.size(); k++) reggeons.push_back(ker_reggeons[k]);
       }
       // Now we need to get the indices of the reggeons right
       for(int i = 0; i < reggeons.size(); i++) reggeons[i].setIndex(i+1);
       Spectra spec(t, reggeons);
       spectrum.push_back(spec);
   }
}

int HQCDP::NumberOfDegreesOfFreedom()
{
    /* 
        Computes the number of degrees of freedom in the fit.
        It is equal to the number of fitting points minus the number of parameters.
        Then we iterate over the ProcessObservable vector,and for each ProcessObservable we compute
        the number of data points. Then we iterate over the kernels and compute the number
        of kernel parameters. Finally we add to the number of kernel parameters the number of gns.
    */
   // Iterate over the processes vector
   int nPoints = 0;
   for(int i = 0; i < processes.size(); i++) nPoints += processes[i]->getDataPts()[0].size();
   // Iterate over the kernel vector;
   int nfitPars = 0;
   for(int i = 0; i < kernels.size(); i++) nfitPars += kernels[i]->KernelPars().size();
   // Add to nfitPars the number of gns
   nfitPars += gns.size();
   // # of degrees of freedom = Number of fitting Points - Number of fitting parameters
   return nPoints - nfitPars;
}

double HQCDP::chi2()
{
    /*
        This method computes the total chi2 of the object.
        It iterates over the list of processes and for each process it computes
        the associated chi2.
        The total chi2 will be the sum of the individual chi2s.
    */
    double chi2 = 0.0;
    for(int i = 0; i < processes.size(); i++)
    {
        ProcessObservable * proc = processes[i];
        std::vector<std::vector<double> > points = proc->expKinematics();
        std::vector<kinStruct> izs = proc->getIzs(points, spectrum);
        std::vector<kinStruct> izsbar = proc->getIzsBar(points, spectrum, gns);
        chi2 += proc->rss(izs, izsbar, points);
    }
    return chi2;
}

double HQCDP::operator()(const std::vector<double> &x)
{
    /* 
        The operator() is used as the fitting function.
        It uses the values of x to update the kernel parameters
        as well the gn parameters.
        Then it updates the spectrum since we have new kernel parameters.
        Afther that we also update the GNs values.
        After that we compute the new chi2 value and it is returned
    */
    // Kernel pars correspond to x(0), x(1), x(2), x(3), x(4)
    std::vector<double> gluon_pars = {x[0], x[1], x[2], x[3], x[4]};
    std::vector<std::vector<double> > kernels_pars = {gluon_pars};
    // GNs correspond to x(5), ..., x(n_Kernel_regs+4)
    const int n_Kernel_regs = kernels[0]->NTrajectories();
    std::vector<double> gpars(n_Kernel_regs, 0);
    for(int i = 0; i < n_Kernel_regs; i++) gpars[i] = x[5+i];
    // Update the spectrum
    computeSpectrum(kernels_pars);
    // Update the GNs
    setGNs(gpars);
    // Compute chi2
    double CHI2 = chi2();
    std::cout << "invls: " << x[0] << " a: " << x[1] << " b: " << x[2] << " c: " << x[3] << " d: " << x[4];
    std::cout << "g1: " << x[5] << " g2: " << x[6] << " g3: " << x[7] <<  " g4: " << x[8] << std::endl;
    std::cout << "chi2: " << CHI2 << std::endl;
    // Return chi2
    return CHI2;
}

void HQCDP::fit(const std::vector<double> &xguess, const double delta)
{
    /*
        This function takes an initial xguess and starts to fit
        the kernel parameters as well the gns to the data available
        in the processes vector. The function we want to minimize
        is the chi2 one and to find the optimum set of parameters
        we will use the Nelder-Mead algorithm.
        After the optimization it computes the chi2 with the
        optimum value of parameters and prints the optimum
        parameters as well the value of chi2
    */
    int ndof  = this->NumberOfDegreesOfFreedom();
    // Optimize using Nelder-Mead algorithm
    std::function<double(std::vector<double>)> f = [this] (const std::vector<double> x) {return this->operator()(x);};
    std::vector<double> xopt = optimFunction(xguess, f, delta);
    // Compute chi2
    // Kernel pars correspond to xopt(0), xopt(1), xopt(2), xopt(3), xopt(4)
    std::vector<double> gluon_pars = {xopt[0], xopt[1], xopt[2], xopt[3], xopt[4]};
    std::vector<std::vector<double> > kernels_pars = {gluon_pars};
    // GNs correspond to x(5), ..., x(n_Kernel_regs+4)
    const int n_Kernel_regs = kernels[0]->NTrajectories();
    std::vector<double> gpars(n_Kernel_regs, 0);
    for(int i = 0; i < n_Kernel_regs; i++) gpars[i] = xopt[5+i];
    // Update the spectrum
    computeSpectrum(kernels_pars);
    // Update the GNs
    setGNs(gpars);
    // Compute chi2
    double CHI2 = chi2();
   // Print the results
   std::cout << "Best chi2 found for:" << std::endl;
   std::cout << "invls:\t" << gluon_pars[0] << std::endl;
   std::cout << "a:\t"     << gluon_pars[1] << std::endl;
   std::cout << "b:\t"     << gluon_pars[2] << std::endl;
   std::cout << "c:\t"     << gluon_pars[3] << std::endl;
   std::cout << "d:\t"     << gluon_pars[4] << std::endl;
   // Print the GNs
   for(int i = 0; i < gpars.size(); i++) std::cout << "g" << std::to_string(i+1) << '\t' << gpars[i] << std::endl;
   std::cout << "chi2:\t"  << CHI2 << std::endl;
   std::cout << "Number of degrees of freedom (Ndof):\t" << ndof << std::endl;
   std::cout << "chi2/Ndof:\t" << CHI2/ndof << std::endl;

}