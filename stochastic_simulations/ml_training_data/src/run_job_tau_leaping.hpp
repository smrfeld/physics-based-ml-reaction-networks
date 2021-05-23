#include "reactions.hpp"

map<string, double> initialize(Params params, int no_ca2i, int no_ip3, double tau_fixed) {
    
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
    TauLeaping g;
    cout << "Running init...." << endl;
    OptionsTauLeaping options;
    options.verbose = false;
    options.conserved_species = {"ca2i", "ip3", "ca2er"};
    options.tau_fixed = tau_fixed;
    options.with_fixed_tau = true;
    auto counts_hist_init = g.run(rxn_list_init, counts_init, dt_st_every_init, t_max_init, options);
    cout << "Finished running init...." << endl;

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

void run_with_diff_eq(string fname_write, Params params, Params_DeYoung_Keizer params_dk, int no_ca2i, int no_ip3, double tau_fixed, bool verbose) {
    
    // Init
    auto counts_init_avg = initialize(params, no_ca2i, no_ip3, tau_fixed);
    
    TauLeaping g;
    
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
    cout << "Running actual" << endl;
    double t_curr = 0.0;
    while (t_curr <= params.times.t_max) {
    
        bool with_fixed_tau = true;
        auto pr0 = g.calculate_no_times_reaction_occurs_in_tau(rxn_list, counts, with_fixed_tau, tau_fixed);
        if (!pr0.first) {
            // Stop, no rxns left
            break;
        }
        double tau = pr0.second;
                
        // Store current counts
        while ((t_st_next < t_curr + tau) && (t_st_next <= params.times.t_max)) {
            
            if (verbose) {
                cout << t_st_next << " / " << params.times.t_max << endl;
            }
            
            // Store
            counts_hist.store_counts(t_st_next, counts);
            
            // Write?
            //...
            
            // Advance
            t_st_next += params.times.dt_st_every;
        }

        while (t_diff_eq_next < t_curr + tau) {
            /*
            if (verbose) {
                cout << "Taking a diff eq step: " << t_diff_eq_next << " until next rxn: " << t_curr + pr.second << endl;
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
        t_curr += tau;

        // Reached the end of time?
        if (t_curr > params.times.t_max) {
            break;
        }
        
        // Update states = do rxns
        g.update_states_in_tau(rxn_list, counts);
    }
    cout << "Finished running actual" << endl;
    
    // Write it
    counts_hist.write_count_hist_all_species_single_file(fname_write);
}
