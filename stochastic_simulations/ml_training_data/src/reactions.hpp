#include "params.hpp"

// ***************
// MARK: - Rates & params
// ***************

struct Rates_IP3R {
    double a1, a2, a3, a4, a5, b1, b2, b3, b4, b5;
    
    Rates_IP3R() {};
    
    Rates_IP3R(Params_DeYoung_Keizer params) {
        a1 = params.a1;
        a2 = params.a2;
        a3 = params.a3;
        a4 = params.a4;
        a5 = params.a5;
        
        // IP3
        b1 = params.d1 * params.a1;

        // Ca2+ (inhibition)
        b2 = params.d2 * params.a2;

        // IP3
        b3 = params.d3 * params.a3;

        // Ca2+ (inhibition)
        b4 = params.d4 * params.a4;

        // Ca2+ (activation)
        b5 = params.d5 * params.a5;
    }
};

struct Vols {
    double vol_cytosol_in_l, vol_er_in_l, vol_total_in_l;
    
    Vols() {};
    
    Vols(Params_DeYoung_Keizer params, int vol_exp) {
        vol_cytosol_in_l = pow(10,-vol_exp);      // L
        vol_er_in_l = params.c1 * vol_cytosol_in_l;        // L
        vol_total_in_l = vol_cytosol_in_l + vol_er_in_l;
        
        std::cout << "The total volume is: " << vol_total_in_l << " L" << std::endl;
        std::cout << "The cytosol volume is: " << vol_cytosol_in_l << " L" << std::endl;
        std::cout << "The ER volume is: " << vol_er_in_l << " L" << std::endl;
    }
};

struct Rates_Transport {
    double ktransport_fwd;
    double ktransport_bkwd;
    
    Rates_Transport() {};
    
    Rates_Transport(Params_DeYoung_Keizer params, double conc_ip3r) {
        ktransport_fwd = params.v1 * pow(conc_ip3r,-3); // 1 / (uM^3 * s)
        ktransport_bkwd = params.c1 * params.v1 * pow(conc_ip3r,-3); // 1 / (uM^3 * s)
    }
};

struct Times {
    double dt_st_every, t_max, dt_diff_eq;
    
    Times() {};
    
    Times(double dt_st_every, double t_max, double dt_diff_eq) {
        this->dt_st_every = dt_st_every;
        this->t_max = t_max;
        this->dt_diff_eq = dt_diff_eq;
    }
};

struct NoInit {
    double vol_cytosol_in_l;
    int no_ip3r, no_ca2er;
    double conc_ip3r_uM, conc_ca2er_uM;
    
    // Init dist
    double mean_conc_ca2i_uM, std_conc_ca2i_uM, std_conc_ip3_uM;
    
    // Number space
    int get_mean_ca2i_number_space() const {
        return conc_in_um_to_no(mean_conc_ca2i_uM, vol_cytosol_in_l);
    }
    
    double get_std_ca2i_number_space() const {
        return conc_in_um_to_no(std_conc_ca2i_uM, vol_cytosol_in_l);
    }
    
    double get_std_ip3_number_space() const {
        return conc_in_um_to_no(std_conc_ip3_uM, vol_cytosol_in_l);
    }
            
    NoInit() {};
    
    NoInit(int no_ip3r, double mean_conc_ca2i_uM, double std_conc_ca2i_uM, double std_conc_ip3_uM, Vols vols, Params_DeYoung_Keizer params_dk) {
        this->std_conc_ca2i_uM = std_conc_ca2i_uM;
        this->std_conc_ip3_uM = std_conc_ip3_uM;
        this->vol_cytosol_in_l = vols.vol_cytosol_in_l;
        this->mean_conc_ca2i_uM = mean_conc_ca2i_uM; // uM

        // IP3Rs
        this->no_ip3r = no_ip3r;
        conc_ip3r_uM = no_to_conc_in_um(no_ip3r, vols.vol_cytosol_in_l);

        // Ca2ER
        conc_ca2er_uM = (params_dk.c0 - mean_conc_ca2i_uM) / params_dk.c1; // uM
        no_ca2er = conc_in_um_to_no(conc_ca2er_uM, vols.vol_er_in_l);
        
        std::cout << "The mean number of Ca2 in cytoplasm is: " << get_mean_ca2i_number_space() << " and in ER is: " << no_ca2er << std::endl;
        
        /*
        conc_ip3r_uM = no_to_conc_in_um(no_ip3r, vols.vol_cytosol_in_l);
        std::cout << "The no of IP3R is: " << no_ip3r << " concentration: " << conc_ip3r_uM << " uM" << std::endl;
        
        this->mean_conc_ca2i_uM = mean_conc_ca2i_uM; // uM
        mean_no_ca2i = conc_in_um_to_no(conc_ca2i_uM, vols.vol_cytosol_in_l);
            
        conc_ca2er_uM = (params_dk.c0 - conc_ca2i_uM) / params_dk.c1; // uM
        no_ca2er = conc_in_um_to_no(conc_ca2er_uM, vols.vol_er_in_l);

        std::cout << "The mean number of Ca2 in cytoplasm is: " << mean_no_ca2i << " and in ER is: " << no_ca2er << std::endl;
         */
    }
    
    int sample_ca2i(std::mt19937 e2) const {
        int mean_no_ca2i = get_mean_ca2i_number_space();
        double std_no_ca2i = get_std_ca2i_number_space();
        
        std::normal_distribution<double> dist_ca2i(mean_no_ca2i, std_no_ca2i);
        return std::round(dist_ca2i(e2));
    }
    
    int sample_ip3(std::mt19937 e2, int mean_no_ip3) const {
        double std_no_ip3 = get_std_ca2i_number_space();
        
        std::normal_distribution<double> dist_ip3(mean_no_ip3, std_no_ip3);
        return std::round(dist_ip3(e2));
    }
};

struct Params {
    Times times;
    Rates_Transport rates_transport;
    Vols vols;
    Rates_IP3R rates_ip3r;
    NoInit no_init;
    
    Params(Times times, Vols vols, Rates_Transport rates_transport, Rates_IP3R rates_ip3r, NoInit no_init) {
        this->times = times;
        this->vols = vols;
        this->rates_transport = rates_transport;
        this->rates_ip3r = rates_ip3r;
        this->no_init = no_init;
    }
};

// ***************
// MARK: - Make reactions IP3R
// ***************

vector<Rxn> make_reactions_ip3r(Rates_IP3R rates, double vol_cytosol_in_l) {
    
    vector<Rxn> rxn_list;
    string r1,r2,p1,p2,name;
    
    // Panel B of Fig 1
    for (auto k=0; k<=1; k++) {
        
        // Bottom left to bottom right
        r1 = "s0" + to_string(k) + "0";
        r2 = "ca2i";
        p1 = "s0" + to_string(k) + "1";
        name = r1 + "_and_" + r2 + "_to_" + p1;
        rxn_list.push_back(Rxn(name, rates.a4 / (pow(10,-6) * AVOGADRO * vol_cytosol_in_l), {r1, r2}, {p1}));
        
        // Bottom right to bottom left
        r1 = "s0" + to_string(k) + "1";
        p1 = "s0" + to_string(k) + "0";
        p2 = "ca2i";
        name = r1 + "_to_" + p1 + "_and_" + p2;
        rxn_list.push_back(Rxn(name, rates.b4, {r1}, {p1,p2}));
        
        // Bottom left to top right
        r1 = "s0" + to_string(k) + "1";
        r2 = "ip3";
        p1 = "s1" + to_string(k) + "1";
        name = r1 + "_and_" + r2 + "_to_" + p1;
        rxn_list.push_back(Rxn(name, rates.a3 / (pow(10,-6) * AVOGADRO * vol_cytosol_in_l), {r1, r2}, {p1}));
        
        // Top right to bottom left
        r1 = "s1" + to_string(k) + "1";
        p1 = "s0" + to_string(k) + "1";
        p2 = "ip3";
        name = r1 + "_to_" + p1 + "_and_" + p2;
        rxn_list.push_back(Rxn(name, rates.b3, {r1}, {p1,p2}));
        
        // Top left to top right
        r1 = "s1" + to_string(k) + "0";
        r2 = "ca2i";
        p1 = "s1" + to_string(k) + "1";
        name = r1 + "_and_" + r2 + "_to_" + p1;
        rxn_list.push_back(Rxn(name, rates.a2 / (pow(10,-6) * AVOGADRO * vol_cytosol_in_l), {r1, r2}, {p1}));
        
        // Top right to top left
        r1 = "s1" + to_string(k) + "1";
        p1 = "s1" + to_string(k) + "0";
        p2 = "ca2i";
        name = r1 + "_to_" + p1 + "_and_" + p2;
        rxn_list.push_back(Rxn(name, rates.b2, {r1}, {p1,p2}));
        
        // Bottom left to top left
        r1 = "s0" + to_string(k) + "0";
        r2 = "ip3";
        p1 = "s1" + to_string(k) + "0";
        name = r1 + "_and_" + r2 + "_to_" + p1;
        rxn_list.push_back(Rxn(name, rates.a1 / (pow(10,-6) * AVOGADRO * vol_cytosol_in_l), {r1, r2}, {p1}));
        
        // Top left to bottom left
        r1 = "s1" + to_string(k) + "0";
        p1 = "s0" + to_string(k) + "0";
        p2 = "ip3";
        name = r1 + "_to_" + p1 + "_and_" + p2;
        rxn_list.push_back(Rxn(name, rates.b1, {r1}, {p1,p2}));
    }
    
    // Panel C of Figure 1
    for (int i=0; i<=1; i++) {
        for (int j=0; j<=1; j++) {
            
            r1 = "s" + to_string(i) + "0" + to_string(j);
            r2 = "ca2i";
            p1 = "s" + to_string(i) + "1" + to_string(j);
            name = r1 + "_and_" + r2 + "_to_" + p1;
            rxn_list.push_back(Rxn(name, rates.a5 / (pow(10,-6) * AVOGADRO * vol_cytosol_in_l), {r1, r2}, {p1}));
            
            r1 = "s" + to_string(i) + "1" + to_string(j);
            p1 = "s" + to_string(i) + "0" + to_string(j);
            p2 = "ca2i";
            name = r1 + "_to_" + p1 + "_and_" + p2;
            rxn_list.push_back(Rxn(name, rates.b5, {r1}, {p1,p2}));
        }
    }
    
    return rxn_list;
}

// ***************
// MARK: - Reactions transport
// ***************

vector<Rxn> make_reactions_transport(Rates_Transport rates_transport, Vols vols) {
    
    vector<Rxn> rxn_list;
    
    double gamma_bkwd = 6 * rates_transport.ktransport_bkwd / pow(pow(10,-6) * AVOGADRO * vols.vol_cytosol_in_l,3);
    double gamma_fwd = 6 * rates_transport.ktransport_fwd / pow(pow(10,-6) * AVOGADRO * vols.vol_cytosol_in_l,3);

    rxn_list.push_back(Rxn("transport_fwd", gamma_fwd, {"s110","s110","s110","ca2er"}, {"s110","s110","s110","ca2i"}));
    rxn_list.push_back(Rxn("transport_bkwd", gamma_bkwd, {"s110","s110","s110","ca2i"}, {"s110","s110","s110","ca2er"}));

    return rxn_list;
}
