/*
 * ==============================================================================================
 * Project P20 (Team T05): Code Red Worm Propagation Simulation Engine
 * Authors:
 *   - Bhavashruthi M (Roll No: 24011103008)
 *   - D M Subashri   (Roll No: 24011103010)
 *
 * Subject: High Performance Computing & Big Data (Semester 5)
 * Berkeley HPC Motif: Monte Carlo Methods (Dwarf #4)
 * Parallel Programming Framework: OpenMP (Shared-Memory Multi-Threading)
 *
 * Reference Publication:
 *   Cliff Changchun Zou, Weibo Gong, Don Towsley,
 *   "Code Red Worm Propagation Modeling and Analysis",
 *   ACM Conference on Computer and Communications Security (CCS), 2002.
 *
 * This implementation models 5 distinct scenarios using parallel stochastic Monte Carlo
 * simulation to quantify saturation dynamics and delayed patching response effects:
 *   - Scenario 1: Baseline Epidemic (Classical Simple Epidemic Model, no congestion, no patching)
 *   - Scenario 2: Congestion Slowdown (Network router packet-flooding throttling)
 *   - Scenario 3: Human Countermeasures (Immediate dynamic immunization/patching from t=0)
 *   - Scenario 4: Two-Factor Worm Model (Both congestion slowdown and dynamic countermeasures)
 *   - Scenario 5: Delayed Patching Response (Countermeasures delayed until t = T_delay)
 * ==============================================================================================
 */

#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <fstream>
#include <string>
#include <algorithm>
#include <iomanip>
#include <cstdint>
#include <omp.h>

using namespace std;

// -----------------------------------------------------------------------------
// HOST EPIDEMIOLOGICAL STATES
// -----------------------------------------------------------------------------
enum HostState : uint8_t {
    SUSCEPTIBLE = 0, // Vulnerable host not yet reached by worm
    INFECTIOUS  = 1, // Actively scanning and infecting others
    REMOVED     = 2  // Immunized, patched, or quarantined (cannot infect / be infected)
};

// -----------------------------------------------------------------------------
// SIMULATION CONFIGURATION
// -----------------------------------------------------------------------------
struct Config {
    int N = 360000;              // Total vulnerable IIS web servers worldwide (Zou et al.)
    int initial_infected = 10;   // Initial bootstrap seed machines
    int max_steps = 2000;        // Discrete simulation time steps

    // Two-Factor Model parameters from Zou et al. (CCS 2002) Section 5 & 6
    double D0 = 10.0;            // Base scanning delay time
    double k1 = 150.0;           // Congestion delay mean coefficient
    double k2 = 70.0;            // Congestion delay variance coefficient
    double n  = 2.0;             // Power sensitivity exponent
    double a  = 0.5;             // Immunization removal fraction C(t) = a * J(t)
    double gamma = 0.012;        // Continuous disinfection rate dR/dt = gamma * I(t) (Section 6, Eq. 14)
    double alpha = 0.003;        // Pre-immunization rate dQ/dt = alpha * S(t) * (J(t)/N) (Section 6, Eq. 15)

    bool enable_congestion   = false; // Factor 2: Throttles probing rate via congestion
    bool enable_immunization = false; // Factor 1: Dynamic human countermeasures

    int immunization_start = 0;       // Delayed patching response start step (T_delay)
    uint64_t seed = 12345ULL;         // Base PRNG seed for reproducibility
};

// -----------------------------------------------------------------------------
// METRICS PER TIME STEP
// -----------------------------------------------------------------------------
struct StepStats {
    int susceptible;
    int infectious;
    int removed;
    int cumulative_infected;
};

// -----------------------------------------------------------------------------
// RESULT OF A SINGLE MONTE CARLO TRAJECTORY
// -----------------------------------------------------------------------------
struct TrialResult {
    int final_infected = 0;
    int peak_infectious = 0;
    int peak_time = 0;
    int final_removed = 0;
    int duration = 0;
    vector<StepStats> history;
};

// -----------------------------------------------------------------------------
// STATISTICAL SUMMARY ACROSS MONTE CARLO TRIALS
// -----------------------------------------------------------------------------
struct ScenarioSummary {
    string name;
    double mean_final_infected = 0.0;
    double mean_peak_infectious = 0.0;
    double mean_peak_time = 0.0;
    double mean_final_removed = 0.0;
    double mean_duration = 0.0;
    double elapsed_wall_time = 0.0;
    double saturation_percentage = 0.0;
};

// -----------------------------------------------------------------------------
// FAST THREAD-SAFE SEED GENERATOR (SplitMix64)
// Prevents false sharing and cross-thread PRNG correlations
// -----------------------------------------------------------------------------
static inline uint64_t splitmix64(uint64_t x) {
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}

// -----------------------------------------------------------------------------
// DYNAMIC INFECTION DELAY SAMPLING (Zou et al. CCS 2002 Equations 11-13)
// -----------------------------------------------------------------------------
static inline int sample_scanning_delay(
    const Config& cfg,
    int cumulative_infected,
    mt19937_64& gen
) {
    if (!cfg.enable_congestion) {
        return max(1, static_cast<int>(ceil(cfg.D0)));
    }

    double p = static_cast<double>(cumulative_infected) / static_cast<double>(cfg.N);
    double p_pow = pow(p, cfg.n);

    double mean = cfg.k1 * p_pow;
    double variance = cfg.k2 * p_pow;
    double stddev = sqrt(variance > 0.0 ? variance : 0.0);

    normal_distribution<double> dist(mean, stddev > 1e-6 ? stddev : 1e-6);
    double X = dist(gen);
    double Y = (X > 0.0) ? ceil(X) : 0.0;

    return max(1, static_cast<int>(ceil(cfg.D0 + Y)));
}

// -----------------------------------------------------------------------------
// DYNAMIC IMMUNIZATION & DISINFECTION (Zou et al. CCS 2002 Section 6)
// 1. dR/dt = gamma * I(t) : continuous cleaning of active infectious hosts (Eq. 14)
// 2. dQ/dt = alpha * S(t) * (J(t)/N) : pre-immunization of susceptible hosts (Eq. 15)
// -----------------------------------------------------------------------------
static inline void apply_immunization_step(
    vector<HostState>& hosts,
    vector<int>& infectious_list,
    int cumulative_infected,
    int& current_susceptible,
    int& current_infectious,
    int& current_removed,
    const Config& cfg,
    mt19937_64& gen
) {
    // 1. Continuous Disinfection of Infectious Hosts (Section 6, Eq. 14: dR/dt = gamma * I)
    if (cfg.gamma > 0.0 && current_infectious > 0) {
        int num_to_disinfect = static_cast<int>(round(cfg.gamma * current_infectious));
        if (num_to_disinfect > current_infectious) num_to_disinfect = current_infectious;

        for (int k = 0; k < num_to_disinfect && !infectious_list.empty(); k++) {
            uniform_int_distribution<size_t> idx_dist(0, infectious_list.size() - 1);
            size_t idx = idx_dist(gen);
            int target = infectious_list[idx];

            if (hosts[target] == INFECTIOUS) {
                hosts[target] = REMOVED;
                current_infectious--;
                current_removed++;
            }
            // Fast O(1) unordered removal from active list
            infectious_list[idx] = infectious_list.back();
            infectious_list.pop_back();
        }
    }

    // 2. Pre-immunization of Susceptible Hosts (Section 6, Eq. 15: dQ/dt = alpha * S * (J/N))
    if (cfg.alpha > 0.0 && current_susceptible > 0 && cumulative_infected > 0) {
        double ratio = static_cast<double>(cumulative_infected) / static_cast<double>(cfg.N);
        int num_to_prepatch = static_cast<int>(round(cfg.alpha * current_susceptible * ratio));
        if (num_to_prepatch > current_susceptible) num_to_prepatch = current_susceptible;

        uniform_int_distribution<int> host_dist(0, cfg.N - 1);
        int patched = 0;
        int attempts = 0;
        int max_attempts = num_to_prepatch * 20 + 50;

        while (patched < num_to_prepatch && attempts < max_attempts && current_susceptible > 0) {
            attempts++;
            int target = host_dist(gen);
            if (hosts[target] == SUSCEPTIBLE) {
                hosts[target] = REMOVED;
                current_susceptible--;
                current_removed++;
                patched++;
            }
        }
    }
}

// -----------------------------------------------------------------------------
// SIMULATE A SINGLE STOCHASTIC MONTE CARLO TRAJECTORY
// -----------------------------------------------------------------------------
TrialResult simulate_single_trajectory(const Config& cfg) {
    vector<HostState> hosts(cfg.N, SUSCEPTIBLE);
    vector<int> next_attempt_time(cfg.N, -1);
    vector<int> infectious_list;
    infectious_list.reserve(cfg.N);

    mt19937_64 gen(cfg.seed);
    uniform_int_distribution<int> target_dist(0, cfg.N - 1);

    int cumulative_infected = cfg.initial_infected;
    int current_susceptible = cfg.N - cfg.initial_infected;
    int current_infectious  = cfg.initial_infected;
    int current_removed     = 0;

    // Bootstrap initial infections
    for (int i = 0; i < cfg.initial_infected; i++) {
        hosts[i] = INFECTIOUS;
        infectious_list.push_back(i);
        int delay = sample_scanning_delay(cfg, cumulative_infected, gen);
        next_attempt_time[i] = delay;
    }

    TrialResult result;
    result.history.reserve(cfg.max_steps);

    // Discrete time step loop
    for (int t = 0; t < cfg.max_steps; t++) {

        // --- Phase 1: Human Countermeasures ---
        if (cfg.enable_immunization && t >= cfg.immunization_start) {
            apply_immunization_step(
                hosts, infectious_list, cumulative_infected,
                current_susceptible, current_infectious, current_removed,
                cfg, gen
            );
        }

        // --- Phase 2: Identify Active Hosts Ready to Probe ---
        vector<int> ready_hosts;
        vector<int> active_infectious;
        active_infectious.reserve(infectious_list.size());

        for (int host : infectious_list) {
            if (hosts[host] == INFECTIOUS) {
                active_infectious.push_back(host);
                if (next_attempt_time[host] <= t) {
                    ready_hosts.push_back(host);
                }
            }
        }
        infectious_list.swap(active_infectious);
        current_infectious = static_cast<int>(infectious_list.size());

        // Outbreak terminates if no infectious hosts remain
        if (current_infectious == 0) {
            result.history.push_back({
                current_susceptible, current_infectious, current_removed, cumulative_infected
            });
            break;
        }

        // --- Phase 3 & 4: Uniform Random Scanning & Infection ---
        vector<int> newly_infected;
        for (size_t i = 0; i < ready_hosts.size(); i++) {
            int target = target_dist(gen);
            if (hosts[target] == SUSCEPTIBLE) {
                hosts[target] = INFECTIOUS;
                cumulative_infected++;
                current_susceptible--;
                current_infectious++;
                newly_infected.push_back(target);
                infectious_list.push_back(target);
            }
        }

        // --- Phase 5: Schedule Next Scanning Timers ---
        for (int host : ready_hosts) {
            if (hosts[host] == INFECTIOUS) {
                int delay = sample_scanning_delay(cfg, cumulative_infected, gen);
                next_attempt_time[host] = t + delay;
            }
        }
        for (int host : newly_infected) {
            int delay = sample_scanning_delay(cfg, cumulative_infected, gen);
            next_attempt_time[host] = t + delay;
        }

        // --- Phase 6: Record Metrics ---
        result.history.push_back({
            current_susceptible, current_infectious, current_removed, cumulative_infected
        });

        if (current_infectious > result.peak_infectious) {
            result.peak_infectious = current_infectious;
            result.peak_time = t;
        }
        result.final_infected = cumulative_infected;
        result.final_removed  = current_removed;
        result.duration       = t + 1;
    }

    if (result.history.empty()) {
        result.history.push_back({
            current_susceptible, current_infectious, current_removed, cumulative_infected
        });
        result.final_infected = cumulative_infected;
        result.final_removed  = current_removed;
        result.duration       = 1;
    }

    return result;
}

// -----------------------------------------------------------------------------
// RUN PARALLEL MONTE CARLO EXPERIMENT FOR A SCENARIO
// -----------------------------------------------------------------------------
ScenarioSummary run_scenario_experiment(
    Config cfg,
    const string& scenario_name,
    const string& output_filename,
    int num_trials
) {
    cout << "  [Running] " << scenario_name << " (" << num_trials << " Monte Carlo trials)..." << flush;

    vector<TrialResult> trials(num_trials);
    double t_start = omp_get_wtime();

    // Parallelize independent stochastic Monte Carlo trials across OpenMP threads
    #pragma omp parallel for schedule(dynamic)
    for (int run = 0; run < num_trials; run++) {
        Config run_cfg = cfg;
        run_cfg.seed = splitmix64(cfg.seed + static_cast<uint64_t>(run * 104729ULL));
        trials[run] = simulate_single_trajectory(run_cfg);
    }

    double t_end = omp_get_wtime();
    double wall_time = t_end - t_start;

    // Aggregate statistics across Monte Carlo runs
    ScenarioSummary summary;
    summary.name = scenario_name;
    summary.elapsed_wall_time = wall_time;

    for (const auto& r : trials) {
        summary.mean_final_infected  += r.final_infected;
        summary.mean_peak_infectious += r.peak_infectious;
        summary.mean_peak_time       += r.peak_time;
        summary.mean_final_removed   += r.final_removed;
        summary.mean_duration        += r.duration;
    }

    summary.mean_final_infected  /= num_trials;
    summary.mean_peak_infectious /= num_trials;
    summary.mean_peak_time       /= num_trials;
    summary.mean_final_removed   /= num_trials;
    summary.mean_duration        /= num_trials;
    summary.saturation_percentage = (summary.mean_final_infected / static_cast<double>(cfg.N)) * 100.0;

    // Determine maximum time steps recorded
    int max_steps = 0;
    for (const auto& r : trials) {
        max_steps = max(max_steps, static_cast<int>(r.history.size()));
    }

    // Export averaged trajectory to CSV
    ofstream out(output_filename);
    if (out.is_open()) {
        out << "time,avg_susceptible,avg_infectious,avg_removed,avg_cumulative_infected\n";
        for (int t = 0; t < max_steps; t++) {
            double s_sum = 0.0, i_sum = 0.0, r_sum = 0.0, c_sum = 0.0;
            for (const auto& r : trials) {
                const auto& s = (t < static_cast<int>(r.history.size())) ? r.history[t] : r.history.back();
                s_sum += s.susceptible;
                i_sum += s.infectious;
                r_sum += s.removed;
                c_sum += s.cumulative_infected;
            }
            out << t << ","
                << fixed << setprecision(2)
                << (s_sum / num_trials) << ","
                << (i_sum / num_trials) << ","
                << (r_sum / num_trials) << ","
                << (c_sum / num_trials) << "\n";
        }
        out.close();
    }

    cout << " Done in " << fixed << setprecision(3) << wall_time << "s | Peak I: "
         << static_cast<int>(round(summary.mean_peak_infectious)) << " ("
         << fixed << setprecision(1) << (summary.mean_peak_infectious / cfg.N * 100.0) << "%) | Saturation: "
         << fixed << setprecision(1) << summary.saturation_percentage << "%\n";

    return summary;
}

// -----------------------------------------------------------------------------
// MAIN ENTRY POINT
// -----------------------------------------------------------------------------
int main(int argc, char** argv) {
    int num_trials = 100;
    if (argc > 1) {
        num_trials = atoi(argv[1]);
        if (num_trials <= 0) num_trials = 100;
    }

    omp_set_dynamic(0);
    int num_threads = omp_get_max_threads();

    cout << "=================================================================================\n";
    cout << "  HPC Project P20 (Team T05): Code Red Worm Propagation Simulation Engine        \n";
    cout << "  Authors: Bhavashruthi M (24011103008) & D M Subashri (24011103010)             \n";
    cout << "  Berkeley HPC Dwarf: Monte Carlo (#4) | Parallel Framework: OpenMP             \n";
    cout << "  Zou et al. (ACM CCS 2002) Two-Factor Worm Model | 5 Comparative Scenarios      \n";
    cout << "=================================================================================\n\n";

    cout << "Execution Environment:\n";
    cout << "  Total Vulnerable Hosts (N) : 360,000 Microsoft IIS servers\n";
    cout << "  Initial Infected Seeds (I0): 10\n";
    cout << "  Monte Carlo Trials / Run   : " << num_trials << "\n";
    cout << "  OpenMP Threads Available   : " << num_threads << "\n\n";

    vector<ScenarioSummary> summaries;
    Config base_cfg;

    // -------------------------------------------------------------------------
    // SCENARIO 1: BASELINE (Classical Simple Epidemic Model)
    // No congestion slowdown, no human countermeasures.
    // -------------------------------------------------------------------------
    Config s1_cfg = base_cfg;
    s1_cfg.enable_congestion   = false;
    s1_cfg.enable_immunization = false;
    s1_cfg.seed = 11111ULL;
    summaries.push_back(run_scenario_experiment(
        s1_cfg, "Scenario 1: Baseline Epidemic", "scenario1_baseline.csv", num_trials
    ));

    // -------------------------------------------------------------------------
    // SCENARIO 2: CONGESTION SLOWDOWN ONLY
    // Network packet-flooding and router queue drops throttle scanning rates.
    // -------------------------------------------------------------------------
    Config s2_cfg = base_cfg;
    s2_cfg.enable_congestion   = true;
    s2_cfg.enable_immunization = false;
    s2_cfg.seed = 22222ULL;
    summaries.push_back(run_scenario_experiment(
        s2_cfg, "Scenario 2: Congestion Slowdown", "scenario2_congestion.csv", num_trials
    ));

    // -------------------------------------------------------------------------
    // SCENARIO 3: HUMAN COUNTERMEASURES ONLY
    // Dynamic immunization & patching from t = 0 without router congestion.
    // -------------------------------------------------------------------------
    Config s3_cfg = base_cfg;
    s3_cfg.enable_congestion   = false;
    s3_cfg.enable_immunization = true;
    s3_cfg.immunization_start  = 0;
    s3_cfg.seed = 33333ULL;
    summaries.push_back(run_scenario_experiment(
        s3_cfg, "Scenario 3: Human Countermeasures", "scenario3_countermeasures.csv", num_trials
    ));

    // -------------------------------------------------------------------------
    // SCENARIO 4: TWO-FACTOR WORM MODEL
    // Both congestion slowdown and human countermeasures from t = 0.
    // -------------------------------------------------------------------------
    Config s4_cfg = base_cfg;
    s4_cfg.enable_congestion   = true;
    s4_cfg.enable_immunization = true;
    s4_cfg.immunization_start  = 0;
    s4_cfg.seed = 44444ULL;
    summaries.push_back(run_scenario_experiment(
        s4_cfg, "Scenario 4: Two-Factor Model", "scenario4_both_factors.csv", num_trials
    ));

    // -------------------------------------------------------------------------
    // SCENARIO 5: DELAYED PATCHING RESPONSE
    // Both factors active, but countermeasures start only after delay (t = 250).
    // Quantifies how delayed reaction alters saturation and peak infection.
    // -------------------------------------------------------------------------
    Config s5_cfg = base_cfg;
    s5_cfg.enable_congestion   = true;
    s5_cfg.enable_immunization = true;
    s5_cfg.immunization_start  = 250; // Critical response delay (T_delay = 250 steps)
    s5_cfg.seed = 55555ULL;
    summaries.push_back(run_scenario_experiment(
        s5_cfg, "Scenario 5: Delayed Patching Response", "scenario5_delayed_patching.csv", num_trials
    ));

    // -------------------------------------------------------------------------
    // EXPORT SUMMARY TABLE & DISPLAY RESULTS
    // -------------------------------------------------------------------------
    cout << "\n=================================================================================\n";
    cout << "                           5-SCENARIO SUMMARY RESULTS                            \n";
    cout << "=================================================================================\n";
    cout << left << setw(38) << "Scenario Name"
         << right << setw(14) << "Peak Infected"
         << setw(12) << "Peak Time"
         << setw(16) << "Final Infected"
         << setw(14) << "Saturation"
         << setw(12) << "Wall Time" << "\n";
    cout << string(106, '-') << "\n";

    ofstream sum_out("simulation_summary.csv");
    if (sum_out.is_open()) {
        sum_out << "scenario,mean_peak_infectious,mean_peak_time,mean_final_infected,saturation_percent,mean_final_removed,wall_time_sec\n";
    }

    for (const auto& s : summaries) {
        cout << left  << setw(38) << s.name
             << right << setw(14) << static_cast<int>(round(s.mean_peak_infectious))
             << setw(12) << static_cast<int>(round(s.mean_peak_time))
             << setw(16) << static_cast<int>(round(s.mean_final_infected))
             << setw(13) << fixed << setprecision(1) << s.saturation_percentage << "%"
             << setw(11) << fixed << setprecision(3) << s.elapsed_wall_time << "s\n";

        if (sum_out.is_open()) {
            sum_out << "\"" << s.name << "\","
                    << s.mean_peak_infectious << ","
                    << s.mean_peak_time << ","
                    << s.mean_final_infected << ","
                    << s.saturation_percentage << ","
                    << s.mean_final_removed << ","
                    << s.elapsed_wall_time << "\n";
        }
    }
    if (sum_out.is_open()) sum_out.close();

    cout << "\nAll 5 scenario datasets successfully generated and exported to CSV files.\n";
    cout << "Open dashboard.html in your browser to inspect the interactive multi-scenario charts.\n";
    cout << "=================================================================================\n";

    return 0;
}

