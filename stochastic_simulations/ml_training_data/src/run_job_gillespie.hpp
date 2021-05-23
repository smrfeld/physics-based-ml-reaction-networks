#include "reactions.hpp"

map<string, double> initialize(Params params, int no_ca2i, int no_ip3) {
    
    // INITIALIZE
    
    // Initial counts of particles
    Counts counts_init = Counts();
    counts_init.set_count("s000", params.no_init.no_ip3r);
    counts_init.set_count("s100", 0);
    counts_init.set_count("s010", 0);
    counts_init.set_count("s001", 0);
    counts_init.set_count("s110", 0);
    counts_init.set_count("s101", 0);
    counts_init.set_count("s011", 0);
    counts_init.set_count("s111", 0);
    counts_init.set_count("ip3", no_ip3);
    counts_init.set_count("ca2i", no_ca2i);
    counts_init.set_count("ca2er", params.no_init.no_ca2er);

    // Make reactions
    vector<Rxn> rxn_list_init = make_reactions_ip3r(params.rates_ip3r, params.vols.vol_cytosol_in_l);

    // Run
    double dt_st_every_init = 0.1;
    double t_max_init = 10.0;
    Gillespie g;
    std::cout << "Running init...." << std::endl;
    OptionsGillespie options;
    options.conserved_species = {"ca2i", "ip3", "ca2er"};
    auto counts_hist_init = g.run(rxn_list_init, counts_init, dt_st_every_init, t_max_init, options);
    std::cout << "Finished running init...." << std::endl;

    // Average
    return counts_hist_init.get_average_counts(0.4 * counts_hist_init.get_t_max());
}

double calculate_no_ca2i_influx_current(int no_ca2i_curr, int no_ca2er_curr, const Params &params, const Params_DeYoung_Keizer &params_dk) {
    double conc_sqrd = pow(no_to_conc_in_um(no_ca2i_curr, params.vols.vol_cytosol_in_l),2);
    double j1, j2;
    j2 = conc_in_um_to_no_as_double(params_dk.v3 * conc_sqrd / (conc_sqrd + pow(params_dk.k3,2)), params.vols.vol_cytosol_in_l);
    j1 = conc_in_um_to_no_as_double(params_dk.c1 * params_dk.v2 * (no_to_conc_in_um(no_ca2er_curr, params.vols.vol_er_in_l) - no_to_conc_in_um(no_ca2i_curr, params.vols.vol_cytosol_in_l)), params.vols.vol_cytosol_in_l);
    double no_ca2i_influx = params.times.dt_diff_eq * (j1 - j2);
    
    return no_ca2i_influx;
}

vector<Rxn> make_rxns(Params params) {
    auto rxn_list_ip3r = make_reactions_ip3r(params.rates_ip3r, params.vols.vol_cytosol_in_l);
    auto rxn_list_transport = make_reactions_transport(params.rates_transport, params.vols);

    vector<Rxn> rxn_list;
    rxn_list.reserve( rxn_list_ip3r.size()); // preallocate memory
    rxn_list.insert( rxn_list.end(), rxn_list_ip3r.begin(), rxn_list_ip3r.end() );
    rxn_list.insert( rxn_list.end(), rxn_list_transport.begin(), rxn_list_transport.end() );
    return rxn_list;
}

void run_with_diff_eq(string fname_write, Params params, Params_DeYoung_Keizer params_dk, int no_ca2i, int no_ip3, bool verbose) {
    
    // Init
    auto counts_init_avg = initialize(params, no_ca2i, no_ip3);
    
    Gillespie g;
    
    // Initial counts of particles
    Counts counts = Counts();
    for (auto const &pr: counts_init_avg) {
        counts.set_count(pr.first, round(pr.second));
    }
    
    // Make reactions
    auto rxn_list = make_rxns(params);
    
    // Setup count storage
    CountsHist counts_hist;
    counts_hist.store_counts(0, counts);
    double t_st_next = params.times.dt_st_every;
    
    // Keep track of the Ca current
    // Calculate the initial current
    double no_ca2i_influx = calculate_no_ca2i_influx_current(counts.get_count("ca2i"), counts.get_count("ca2er"), params, params_dk);

    // Next time to store do diff eq
    double t_diff_eq_next = params.times.dt_diff_eq;
    
    // Run
    std::cout << "Running actual" << std::endl;
    double t_curr = 0.0;
    while (t_curr <= params.times.t_max) {
    
        // Pick next rxn
        auto pr = g.choose_next_rxn(rxn_list, counts);
                
        // Check if exists
        if (pr.first == nullptr) {
            break;
        }
        
        // Store
        while (t_st_next < t_curr + pr.second) {
            
            if (verbose) {
                std::cout << t_st_next << " / " << params.times.t_max << " ca2i: " << counts.get_count("ca2i") << std::endl;
            }
            
            // Store
            counts_hist.store_counts(t_st_next, counts);
                        
            // Advance
            t_st_next += params.times.dt_st_every;
        }
        
        while (t_diff_eq_next < t_curr + pr.second) {
            /*
            if (verbose) {
                std::cout << "Taking a diff eq step: " << t_diff_eq_next << " until next rxn: " << t_curr + pr.second << std::endl;
            }
             */
            
            // Influx/efflux
            if (no_ca2i_influx >= 1) {
                // Influx
                int no_ca2i_influx_int = min(int(no_ca2i_influx), counts.get_count("ca2er"));
                counts.increment_count("ca2i", no_ca2i_influx_int, true);
                counts.increment_count("ca2er", -1 * no_ca2i_influx_int, true);
                no_ca2i_influx -= double(no_ca2i_influx_int);
            } else if (no_ca2i_influx <= -1) {
                // Efflux
                int no_ca2i_efflux_int = min(int(-1 * no_ca2i_influx), counts.get_count("ca2i"));
                counts.increment_count("ca2i", -1 * no_ca2i_efflux_int, true);
                counts.increment_count("ca2er", no_ca2i_efflux_int, true);
                no_ca2i_influx += double(no_ca2i_efflux_int);
            }

            // Add for next
            no_ca2i_influx += calculate_no_ca2i_influx_current(counts.get_count("ca2i"), counts.get_count("ca2er"), params, params_dk);

            // Advance
            t_diff_eq_next += params.times.dt_diff_eq;
        }

        // Advance time
        t_curr += pr.second;
        
        // Do the reaction
        g.do_rxn(pr.first, counts);
    }
    std::cout << "Finished running actual" << std::endl;

    // Write it
    counts_hist.write_count_hist_all_species_single_file(fname_write);
}
