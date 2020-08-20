#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include "SigmaProtonProton.h"

using namespace std;

int main(int argc, char ** argv)
{
    string data_path;
    if(argc < 2)
    {
        data_path = "expdata/SigmaProtonProton/SigmaProtonProton_data_W_lt_10000.txt";
    }
    else
    {
        data_path = argv[1];
    }
    SigmaProtonProton sigma_pp(data_path);

    // We now start to test the SigmaProtonProton class
    // Test methods related to data loaded
    const double mb_to_GEVMINUS2 = 10 / 3.894 ;
    const double mrho_U_mrho_GeV = 4.3669;
    const double mb_to_UMINUS2 = mb_to_GEVMINUS2 / std::pow(mrho_U_mrho_GeV,2);
    vector<double> sigma_exp_val = sigma_pp.expVal();
    vector<double> sigma_exp_err = sigma_pp.expErr();
    vector<vector<double> > sigma_kinematics = sigma_pp.expKinematics();
    cout << "W(GeV)\tWPlus\tWMinus\tsigma(mb)\tError" << endl;
    for(int i = 0; i < sigma_exp_val.size(); i++) 
    {
        cout << sigma_kinematics[0][i] / mrho_U_mrho_GeV  << '\t' << sigma_kinematics[1][i] / mrho_U_mrho_GeV  << '\t' << sigma_kinematics[2][i] / mrho_U_mrho_GeV << '\t';
        cout << sigma_exp_val[i] / mb_to_UMINUS2 << '\t' << sigma_exp_err[i] / mb_to_UMINUS2 << endl;
    }

    // Create 4 reggeons just to test
    Reggeon reg1(1.4, 0, {{}}, "gluon", 1), reg2(1.08, 0, {{}}, "gluon", 2), reg3(0.45, 0, {{}}, "meson", 3), reg4(0.45, 0, {{}}, "meson", 4);
    // Test IzN and IzNBar methods
    cout << "Testing IzN" << endl;
    cout << "IzN1\tIzN2\tIzN3\tIzN4" << endl;
    cout << sigma_pp.IzN({}, reg1) << '\t' << sigma_pp.IzN({}, reg2) << '\t' << sigma_pp.IzN({}, reg3) << '\t' <<  sigma_pp.IzN({}, reg4) << endl;
    cout << "Testing IzNBar" << endl;
    vector<double> gs = {1, 2, 3, 4};
    vector<double> Ws_vals = sigma_kinematics[0];
    cout << "W(GeV)\tIzNBar1\tIzNBar2\tIzNBar3\tIzNBar4" << endl;
    for(int i = 0; i < Ws_vals.size(); i++)
    {
       double  W_GeV = Ws_vals[i] / mrho_U_mrho_GeV;
       cout << W_GeV << '\t' << sigma_pp.IzNBar({W_GeV}, reg1, gs) << '\t' << sigma_pp.IzNBar({W_GeV}, reg2, gs) << '\t' << sigma_pp.IzNBar({W_GeV}, reg3,  gs);
       cout << '\t' << sigma_pp.IzNBar({W_GeV}, reg4, gs) << endl; 
    }

    return 0;
}
