#include "run_job_tau_leaping.hpp"

// ***************
// MARK: - Main
// ***************

void launch_jobs(vector<double> concs_ip3, int fname_idx_start, int no_fname_idx, std::string data_dir, const Params &params, const Params_DeYoung_Keizer &params_dk);

int main() {
        
    // Random seed
    srand (time(NULL));
    // srand(10);
    
    // ***************
    // MARK: - De Young Keizer rates
    // ***************
  
    Params_DeYoung_Keizer params_dk;
    
    // ***************
    // MARK: - Volumes
    // ***************
        
    int vol_exp = 14;
    Vols vols = Vols(params_dk, vol_exp);
        
    // ***************
    // MARK: - Init
    // ***************
    
    vector<int> no_ip3rs({100,1000});
    for (auto no_ip3r: no_ip3rs) {
        
        // int no_ip3r = 3000;
        double mean_conc_ca2i_uM = 0.25;
        double std_conc_ca2i_uM = 0.001;
        double std_conc_ip3_uM = 0.001;
        NoInit no_init(no_ip3r, mean_conc_ca2i_uM, std_conc_ca2i_uM, std_conc_ip3_uM, vols, params_dk);
            
        // ***************
        // MARK: - Rates
        // ***************
        
        Rates_Transport rates_transport(params_dk, no_init.conc_ip3r_uM);
        Rates_IP3R rates_ip3r(params_dk);
        
        // ***************
        // MARK: - Times
        // ***************
        
        double dt_st_every = 0.1;
        double t_max = 50.0; // 100.0;
        double dt_diff_eq = 0.001;
        Times times(dt_st_every, t_max, dt_diff_eq);
        
        // ***************
        // MARK: - Combine into params object
        // ***************
        
        Params params(times, vols, rates_transport, rates_ip3r, no_init);
        
        // ***************
        // MARK: - Writing directory
        // ***************
        
        // Clear dir for writing
        string data_dir = "../data_tau_leaping/";
        if (string(1,data_dir.back()) != "/") {
            data_dir += "/";
        }

        // Also put the ip3r count
        data_dir += "vol_exp_" + pad_str(vol_exp, 2) + "/ip3r_" + pad_str(no_init.no_ip3r,5) + "/";
        system(("mkdir -p " + data_dir).c_str());

        // Clear the directory...
        if (false) {
            system(("rm -rf " + data_dir + "*").c_str());
        }
        
        // ***************
        // MARK: - Concentrations of IP3 to do, no samples for each
        // ***************
        
        // Concentrations of ip3 to do
        vector<double> concs_ip3({
            // 0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0
            1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0
        });
        
        // No runs for each
        int fname_idx_start = 0;
        int no_fname_idx = 100;
        
        // ***************
        // MARK: - Run
        // ***************
        
        launch_jobs(concs_ip3, fname_idx_start, no_fname_idx, data_dir, params, params_dk);
    }
    
    return 0;
}

// ***************
// MARK: - Launch jobs
// ***************

void launch_jobs(vector<double> concs_ip3, int fname_idx_start, int no_fname_idx, std::string data_dir, const Params &params, const Params_DeYoung_Keizer &params_dk) {
    
    // Random
    random_device rd;
    mt19937 e2(rd());
    
    // Setup all combinations
    vector<pair<int, int>> jobs;
    for (auto idx_ip3=0; idx_ip3<concs_ip3.size(); idx_ip3++) {
        for (auto fname_idx=fname_idx_start; fname_idx<fname_idx_start+no_fname_idx; fname_idx++) {
            jobs.push_back(make_pair(idx_ip3, fname_idx));
        }
    }
    
    // Number of jobs to run at a time
    int no_threads = 8; // 8
    int idx_next_job = 0;
    vector<thread> threads;
    while (idx_next_job < jobs.size()) {
        int no_jobs_left = jobs.size()-idx_next_job;
        int no_threads_launch = min(no_threads, no_jobs_left);
        
        cout << "Launching jobs: " << idx_next_job << " to " << idx_next_job+no_threads_launch << " of " << jobs.size() << endl;
        
        // Launch the next no_threads
        for (auto idx_thread=0; idx_thread<no_threads_launch; idx_thread++) {
            auto pr = jobs[idx_next_job + idx_thread];
            double conc_ip3 = concs_ip3.at(pr.first);
            int mean_no_ip3 = conc_in_um_to_no(conc_ip3, params.vols.vol_cytosol_in_l);
            int fname_idx = pr.second;

            cout << "   IP3 conc: " << conc_ip3 << " fname_idx: " << fname_idx << endl;
            
            // Sample the actual no particles
            int no_ca2i = params.no_init.sample_ca2i(e2);
            int no_ip3 = params.no_init.sample_ip3(e2,mean_no_ip3);
            
            // Dir name to write to
            string dir_name_write = data_dir + "ip3_" + prec_p_str(conc_ip3, 3) + "/";
            system(("mkdir -p " + dir_name_write).c_str());
            string fname_write = dir_name_write + pad_str(fname_idx, 4) + ".txt";
            
            // Launch job
            bool verbose = false;
            double tau_fixed = 0.0001;
            threads.push_back(thread(run_with_diff_eq, fname_write, params, params_dk, no_ca2i, no_ip3, tau_fixed, verbose));
        }

        // Wait for the threads to join
        for (auto idx_thread=0; idx_thread<no_threads_launch; idx_thread++) {
            threads[idx_thread].join();
        }
        
        // Clear threads for next time
        threads.clear();
        
        // Advance
        idx_next_job += no_threads_launch;
    }
}
