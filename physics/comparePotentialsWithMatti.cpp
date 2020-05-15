#include <iostream>
#include <string>
#include <fstream>
#include <boost/algorithm/string.hpp>
#include "HolographicVQCD.h"
#include "methods/interpolation/Poly_Interp.hpp"

using namespace std;

int main(int argc, char ** argv)
{
    double sc, ksc, wsc, W0, w0, kU1, wU1;
    double VgIR, WIR, kIR, wIR, W1, k1, w1;
    double xf = 2.0/3, tau0, Za, ca;
    if (argc < 18)
    {
        sc = 3.0; ksc = 3.0; wsc = 1.56; W0 = 2.5; w0 = 1.26; kU1 = 11./9; wU1 = 0.0;
        VgIR = 2.05; WIR = 0.9; kIR = 1.8; wIR = 5.0; W1 = 0.0; k1 = -0.23;
        w1 = 0.0; tau0 = 1.; Za = 133; ca = 0.26;
    }
    else
    {
        sc = stod(argv[1]); ksc = stod(argv[2]); wsc = stod(argv[3]); W0 = stod(argv[4]); w0 = stod(argv[5]);
        kU1 = stod(argv[6]);
        wU1 = stod(argv[7]); VgIR = stod(argv[8]); WIR = stod(argv[9]); kIR = stod(argv[10]); wIR = stod(argv[11]);
        W1 = stod(argv[12]); k1 = stod(argv[13]); w1 = stod(argv[14]); tau0 = stod(argv[15]);
        Za = stod(argv[16]); ca = stod(argv[17]);
    }
    
    cout << "Computing potentials with parameter values" << endl;
    cout << "sc: " << sc << " ksc: " << ksc << " wsc: " << wsc << " W0: " << W0 << " w0: " << w0 << " kU1: " << kU1;
    cout << " wU1: " << wU1 << " VgIR: " << VgIR << " WIR: " << WIR << " kIR: " << kIR << " wIR: " << wIR << " W1: " << W1;
    cout << " k1: " << k1 << " w1: " << w1 << " tau0: " << tau0 << " Za: " << Za << " ca: " << ca << endl;

    HVQCD hvqcd(sc, ksc, wsc, W0, w0, kU1, wU1, VgIR, WIR, kIR, wIR, W1, k1, w1, xf, tau0, Za, ca);
    hvqcd.solve(-80, 20);

    // Get u and compute the pseudoscalar potential
    vector<double> VPSM = computePseudoScalarMesonPotential(hvqcd), u = hvqcd.u();
    Poly_Interp<double> VPSM_func(u, VPSM, 4);

    
    ofstream potenial_values;
    ifstream matti_potenial_values;
    potenial_values.open("PseudoscalarPotential.txt");
    matti_potenial_values.open("pseudoscalar_potential.dat");

    potenial_values << "u\tMatti\tMine\tError(%)" << endl;

    string line;
    vector<string> result;
    double matti_uval, matti_pot_val, my_val, error;
    while ( getline (matti_potenial_values,line) )
    {
        boost::split(result, line, boost::is_any_of("\t") );
        matti_uval = stod(result[0]);
        matti_pot_val = stod(result[1]);
        my_val = VPSM_func.interp(matti_uval);
        error = 100 * fabs((matti_pot_val-my_val)/matti_pot_val); 
        potenial_values << matti_uval << '\t' << matti_pot_val << '\t' << my_val << '\t' << error << endl;
    }
    //for(int i = 0; i < u.size(); i++) potenial_values << u[i] << '\t' << VPSM[i] << endl;

    matti_potenial_values.close();
    potenial_values.close();

    return 0;
}