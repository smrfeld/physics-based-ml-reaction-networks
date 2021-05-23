#include <iostream>
#include <vector>
#include <thread>
#include <random>

#include <gilsp>

using namespace gilsp;
using namespace std;

// ***************
// MARK: - Make reactions IP3R
// ***************

struct Params_DeYoung_Keizer {
    
    // Total ca2+ in terms of cytosolic vol
    double c0 = 2.0; // uM

    // ER vol / cytosol vol ratio
    double c1 = 0.185;

    // Max ca2+ channel flux
    double v1 = 6; // 1/s

    // Ca2+ leak flux constant
    double v2 = 0.11; // 1/s

    // Max ca2+ uptake
    double v3 = 0.9; // 1/(uM * s)

    // Activation constant for ATP-ca2+ pump
    double k3 = 0.1; // uM

    // IP3
    double a1 = 400; // 1/(uM * s)

    // Ca2+ (inhibiton)
    double a2 = 0.2; // 1/(uM * s)

    // IP3
    double a3 = 400; // 1/(uM * s)

    // Ca2+ (inhibition)
    double a4 = 0.2; // 1/(uM * s)

    // Ca2+ (activation)
    double a5 = 20; // 1/(uM * s)

    // IP3
    double d1 = 0.13; // uM

    // Ca2+ (inhibition)
    double d2 = 1.049; // uM

    // IP3
    double d3 = 943.4 * pow(10,-3); // uM

    // Ca2+ (inhibition)
    double d4 = 144.5 * pow(10,-3); // uM

    // Ca2+ (activation)
    double d5 = 82.34 * pow(10,-3); // uM    
};
