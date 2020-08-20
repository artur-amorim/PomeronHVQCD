#include <algorithm>
#include <stdexcept>
#include "HQCDP.h"

HQCDP::HQCDP():
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
        if(proc_ts.size() == 0) std::cout << "Process with no t values has been found" << std::endl;
        // Get all the tvals in proc_ts
        for(int j = 0 ; j < proc_ts.size(); j++) useTVals.push_back(proc_ts[j]) ;
    }
    // Now we std::sort the values of t and remove repeated values
    std::sort(useTVals.begin(), useTVals.end()) ;
    // Erase any duplicates
    useTVals.erase(std::unique(useTVals.begin(), useTVals.end() ), useTVals.end() );
}

void HQCDP::addProcessObservable(Process &proc)
{
    // Add a process to be fitted or predicted
    // For consistency change NMC, alpha and rsslog data members of proc
    // to be the same as the HQCDP object.
    Process * newProcess = &proc;
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
   if(kernelPars.size() == kernels.size())
   {
       for(int ker_idx = 0; ker_idx < kernels.size(); ker_idx+=1) kernels[ker_idx]->computeReggeTrajectories(kernelPars[ker_idx]);
   }
   else
   {
       // Computes with default values
       for(int ker_idx = 0; ker_idx < kernels.size(); ker_idx+=1) kernels[ker_idx]->computeReggeTrajectories();
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
            std::cout << "Computing Reggeons of " << ker->Name() << " kernel." << std::endl;
            const int n_regs = ker->NTrajectories();
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
        Then we iterate over the Process * vector,and for each Process we compute
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
        Process * proc = processes[i];
        std::vector<std::vector<double> > points = proc->expKinematics();
        std::vector<kinStruct> izs = proc->getIzs(points, spectrum);
        std::vector<kinStruct> izsbar = proc->getIzsBar(points, spectrum, gns);
        chi2 += proc->chi2(izs, izsbar, points);
    }
    return chi2;
}