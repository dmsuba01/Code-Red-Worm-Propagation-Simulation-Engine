
#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <fstream>
#include <string>
#include <algorithm>
#include <cstdint>
#include <omp.h>

using namespace std;

// ---------------------------------------------
// HOST STATES
// ---------------------------------------------

enum HostState {
    SUSCEPTIBLE,
    INFECTIOUS,
    REMOVED
};

// ---------------------------------------------
// SIMULATION CONFIGURATION
// ---------------------------------------------

struct Config {
    int N = 360000;
    int initial_infected = 10;
    int max_steps = 2000;

    // Parameters from the paper
    double D0 = 10.0;
    double k1 = 150.0;
    double k2 = 70.0;
    double n = 2.0;
    double a = 0.5;

    bool enable_congestion = false;
    bool enable_immunization = false;

    // Time at which immunization starts
    int immunization_start = 0;

    uint64_t seed = 12345;
};

// ---------------------------------------------
// STATISTICS FOR EACH TIME STEP
// ---------------------------------------------

struct StepStats {
    int susceptible;
    int infectious;
    int removed;
    int cumulative_infected;
};

// ---------------------------------------------
// RESULT OF ONE SIMULATION
// ---------------------------------------------

struct Result {
    int final_infected = 0;
    int peak_infectious = 0;
    int peak_time = 0;
    int final_removed = 0;
    int duration = 0;

    vector<StepStats> history;
};

// ---------------------------------------------
// AVERAGE RESULTS FOR ONE SCENARIO
// ---------------------------------------------

struct ScenarioSummary {
    double mean_final_infected = 0;
    double mean_peak_infectious = 0;
    double mean_peak_time = 0;
    double mean_final_removed = 0;
    double mean_duration = 0;
};

// ---------------------------------------------
// GENERATE REPRODUCIBLE RANDOM SEEDS
// ---------------------------------------------

uint64_t make_seed(uint64_t x) {

    x += 0x9e3779b97f4a7c15ULL;

    x = (x ^ (x >> 30)) *
        0xbf58476d1ce4e5b9ULL;

    x = (x ^ (x >> 27)) *
        0x94d049bb133111ebULL;

    return x ^ (x >> 31);
}

// ---------------------------------------------
// SAMPLE RANDOMIZED SCANNING DELAY
// ---------------------------------------------

int sample_scanning_delay(
    const Config& cfg,
    int cumulative_infected,
    mt19937_64& gen
) {

    // Without congestion, use the baseline delay.
    if (!cfg.enable_congestion) {

        return max(
            1,
            static_cast<int>(ceil(cfg.D0))
        );
    }

    // Proportion of cumulatively infected hosts.
    double p =
        static_cast<double>(cumulative_infected) / cfg.N;

    double p_power = pow(p, cfg.n);

    // Normal distribution parameters.
    double mean = cfg.k1 * p_power;

    double variance = cfg.k2 * p_power;

    double standard_deviation = sqrt(variance);

    normal_distribution<double> distribution(
        mean,
        standard_deviation
    );

    // Generate random variable X(t).
    double X = distribution(gen);

    // Calculate Y(t).
    double Y = (X > 0.0) ? ceil(X) : 0.0;

    // Calculate D(t) = D0 + Y(t).
    return max(
        1,
        static_cast<int>(ceil(cfg.D0 + Y))
    );
}

// ---------------------------------------------
// DYNAMIC HUMAN IMMUNIZATION
// ---------------------------------------------

void apply_immunization(
    vector<HostState>& hosts,
    int cumulative_infected,
    const Config& cfg,
    mt19937_64& gen
) {

    // Target number of removed hosts.
    int target_removed = static_cast<int>(
        floor(cfg.a * cumulative_infected)
    );

    int current_removed = 0;

    vector<int> eligible_hosts;

    // Identify hosts that have not been removed.
    for (int i = 0; i < cfg.N; i++) {

        if (hosts[i] == REMOVED) {
            current_removed++;
        }
        else {
            eligible_hosts.push_back(i);
        }
    }

    int additional_to_remove =
        target_removed - current_removed;

    if (additional_to_remove <= 0) {
        return;
    }

    // Randomly select hosts for immunization.
    shuffle(
        eligible_hosts.begin(),
        eligible_hosts.end(),
        gen
    );

    additional_to_remove = min(
        additional_to_remove,
        static_cast<int>(eligible_hosts.size())
    );

    for (int i = 0; i < additional_to_remove; i++) {

        hosts[eligible_hosts[i]] = REMOVED;
    }
}

// ---------------------------------------------
// RUN ONE STOCHASTIC SIMULATION
// ---------------------------------------------

Result simulate(const Config& cfg) {

    // Initialize all hosts as susceptible.
    vector<HostState> hosts(
        cfg.N,
        SUSCEPTIBLE
    );

    // Each host has its own next infection-attempt time.
    // -1 means no attempt is scheduled.
    vector<int> next_attempt_time(
        cfg.N,
        -1
    );

    // Independent random generator for this trial.
    mt19937_64 gen(cfg.seed);

    uniform_int_distribution<int> target_distribution(
        0,
        cfg.N - 1
    );

    // -----------------------------------------
    // INITIAL INFECTIONS
    // -----------------------------------------

    int cumulative_infected =
        cfg.initial_infected;

    for (int i = 0; i < cfg.initial_infected; i++) {

        hosts[i] = INFECTIOUS;

        // Assign an individual scanning timer.
        int delay = sample_scanning_delay(
            cfg,
            cumulative_infected,
            gen
        );

        next_attempt_time[i] = delay;
    }

    Result result;

    // -----------------------------------------
    // MAIN SIMULATION LOOP
    // -----------------------------------------

    for (int t = 0; t < cfg.max_steps; t++) {

        // =====================================
        // PHASE 1: HUMAN COUNTERMEASURES
        // =====================================

        if (cfg.enable_immunization &&
            t >= cfg.immunization_start) {

            apply_immunization(
                hosts,
                cumulative_infected,
                cfg,
                gen
            );
        }

        // =====================================
        // PHASE 2: FIND INFECTIOUS HOSTS
        // =====================================

        vector<int> infectious_hosts;
        vector<int> ready_hosts;

        for (int i = 0; i < cfg.N; i++) {

            if (hosts[i] == INFECTIOUS) {

                infectious_hosts.push_back(i);

                // Check whether its scanning timer expired.
                if (next_attempt_time[i] <= t) {

                    ready_hosts.push_back(i);
                }
            }
        }

        int current_infectious =
            static_cast<int>(infectious_hosts.size());

        // Stop if there are no infectious hosts.
        if (current_infectious == 0) {
            break;
        }

        // =====================================
        // PHASE 3: GENERATE INFECTION ATTEMPTS
        // =====================================

        vector<int> targets(
            ready_hosts.size(),
            -1
        );

        for (
            int i = 0;
            i < static_cast<int>(ready_hosts.size());
            i++
        ) {

            // Select a random target IP address.
            targets[i] = target_distribution(gen);
        }

        // =====================================
        // PHASE 4: COMMIT INFECTIONS
        // =====================================

        vector<int> newly_infected;

        for (
            int i = 0;
            i < static_cast<int>(ready_hosts.size());
            i++
        ) {

            int target = targets[i];

            // Infect only susceptible hosts.
            if (hosts[target] == SUSCEPTIBLE) {

                hosts[target] = INFECTIOUS;

                cumulative_infected++;

                newly_infected.push_back(target);
            }
        }

        // =====================================
        // PHASE 5: SCHEDULE NEXT ATTEMPTS
        // =====================================

        // Schedule the next attempt for hosts
        // that have just completed an attempt.
        for (int host : ready_hosts) {

            int delay = sample_scanning_delay(
                cfg,
                cumulative_infected,
                gen
            );

            next_attempt_time[host] = t + delay;
        }

        // Schedule the first attempt for newly
        // infected hosts.
        for (int host : newly_infected) {

            int delay = sample_scanning_delay(
                cfg,
                cumulative_infected,
                gen
            );

            next_attempt_time[host] = t + delay;
        }

        // =====================================
        // PHASE 6: COLLECT STATISTICS
        // =====================================

        int susceptible = 0;
        int infectious = 0;
        int removed = 0;

        for (int i = 0; i < cfg.N; i++) {

            if (hosts[i] == SUSCEPTIBLE) {

                susceptible++;
            }
            else if (hosts[i] == INFECTIOUS) {

                infectious++;
            }
            else {

                removed++;
            }
        }

        result.history.push_back({
            susceptible,
            infectious,
            removed,
            cumulative_infected
        });

        // Track peak infectious population.
        if (infectious > result.peak_infectious) {

            result.peak_infectious = infectious;

            result.peak_time = t;
        }

        result.final_infected =
            cumulative_infected;

        result.final_removed =
            removed;

        result.duration = t + 1;
    }

    // -----------------------------------------
    // HANDLE EMPTY SIMULATION HISTORY
    // -----------------------------------------

    if (result.history.empty()) {

        int susceptible = 0;
        int infectious = 0;
        int removed = 0;

        for (HostState state : hosts) {

            if (state == SUSCEPTIBLE) {
                susceptible++;
            }
            else if (state == INFECTIOUS) {
                infectious++;
            }
            else {
                removed++;
            }
        }

        result.history.push_back({
            susceptible,
            infectious,
            removed,
            cumulative_infected
        });

        result.final_infected =
            cumulative_infected;

        result.final_removed =
            removed;
    }

    return result;
}

// ---------------------------------------------
// RUN 100 INDEPENDENT TRIALS IN PARALLEL
// ---------------------------------------------

ScenarioSummary run_scenario(
    Config cfg,
    const string& scenario_name,
    const string& filename,
    int num_runs
) {

    cout << "\nRunning "
         << scenario_name
         << " (" << num_runs << " trials)...\n";

    vector<Result> results(num_runs);

    double start_time = omp_get_wtime();

    // Parallelize independent simulation trials.
    #pragma omp parallel for schedule(dynamic)
    for (int run = 0; run < num_runs; run++) {

        Config run_cfg = cfg;

        // Assign an independent reproducible seed.
        run_cfg.seed = make_seed(
            cfg.seed + static_cast<uint64_t>(run)
        );

        results[run] = simulate(run_cfg);
    }

    double end_time = omp_get_wtime();

    // -----------------------------------------
    // CALCULATE AVERAGE RESULTS
    // -----------------------------------------

    ScenarioSummary summary;

    for (const Result& r : results) {

        summary.mean_final_infected +=
            r.final_infected;

        summary.mean_peak_infectious +=
            r.peak_infectious;

        summary.mean_peak_time +=
            r.peak_time;

        summary.mean_final_removed +=
            r.final_removed;

        summary.mean_duration +=
            r.duration;
    }

    summary.mean_final_infected /= num_runs;

    summary.mean_peak_infectious /= num_runs;

    summary.mean_peak_time /= num_runs;

    summary.mean_final_removed /= num_runs;

    summary.mean_duration /= num_runs;

    // -----------------------------------------
    // FIND LONGEST SIMULATION
    // -----------------------------------------

    int max_duration = 0;

    for (const Result& r : results) {

        max_duration = max(
            max_duration,
            static_cast<int>(r.history.size())
        );
    }

    // -----------------------------------------
    // EXPORT AVERAGE RESULTS TO CSV
    // -----------------------------------------

    ofstream out(filename);

    out << "time,avg_susceptible,avg_infectious,"
           "avg_removed,avg_cumulative_infected\n";

    for (int t = 0; t < max_duration; t++) {

        double susceptible_sum = 0;
        double infectious_sum = 0;
        double removed_sum = 0;
        double cumulative_sum = 0;

        for (const Result& r : results) {

            const vector<StepStats>& history =
                r.history;

            // Carry forward the final state of a
            // simulation that has already terminated.
            const StepStats& s =
                (t < static_cast<int>(history.size()))
                ? history[t]
                : history.back();

            susceptible_sum += s.susceptible;

            infectious_sum += s.infectious;

            removed_sum += s.removed;

            cumulative_sum += s.cumulative_infected;
        }

        out << t << ","
            << susceptible_sum / num_runs << ","
            << infectious_sum / num_runs << ","
            << removed_sum / num_runs << ","
            << cumulative_sum / num_runs << "\n";
    }

    out.close();

    // -----------------------------------------
    // DISPLAY SCENARIO RESULTS
    // -----------------------------------------

    cout << "Execution time: "
         << end_time - start_time
         << " seconds\n";

    cout << "Mean final infected: "
         << summary.mean_final_infected << "\n";

    cout << "Mean peak infectious: "
         << summary.mean_peak_infectious << "\n";

    cout << "Mean peak time: "
         << summary.mean_peak_time << "\n";

    cout << "Mean final removed: "
         << summary.mean_final_removed << "\n";

    cout << "Mean duration: "
         << summary.mean_duration << "\n";

    cout << "CSV saved: "
         << filename << "\n";

    return summary;
}

// ---------------------------------------------
// MAIN FUNCTION
// ---------------------------------------------

int main() {

    omp_set_dynamic(0);

    Config cfg;

    const int NUM_RUNS = 100;

    cout << "========================================\n";
    cout << " STOCHASTIC WORM PROPAGATION SIMULATION\n";
    cout << " OpenMP Parallel Implementation\n";
    cout << "========================================\n";

    cout << "Total hosts: "
         << cfg.N << "\n";

    cout << "Initially infected: "
         << cfg.initial_infected << "\n";

    cout << "Trials per scenario: "
         << NUM_RUNS << "\n";

    cout << "OpenMP threads: "
         << omp_get_max_threads() << "\n";

    // =========================================
    // SCENARIO 1: BASELINE
    // No congestion, no human countermeasures.
    // =========================================

    cfg.enable_congestion = false;

    cfg.enable_immunization = false;

    cfg.immunization_start = 0;

    cfg.seed = 12345;

    ScenarioSummary s1 = run_scenario(
        cfg,
        "Scenario 1 - Baseline",
        "scenario1_baseline.csv",
        NUM_RUNS
    );

    // =========================================
    // SCENARIO 2: CONGESTION
    // Congestion only.
    // =========================================

    cfg.enable_congestion = true;

    cfg.enable_immunization = false;

    cfg.seed = 23456;

    ScenarioSummary s2 = run_scenario(
        cfg,
        "Scenario 2 - Congestion",
        "scenario2_congestion.csv",
        NUM_RUNS
    );

    // =========================================
    // SCENARIO 3: HUMAN COUNTERMEASURES
    // Immunization only.
    // =========================================

    cfg.enable_congestion = false;

    cfg.enable_immunization = true;

    cfg.immunization_start = 0;

    cfg.seed = 34567;

    ScenarioSummary s3 = run_scenario(
        cfg,
        "Scenario 3 - Human Countermeasures",
        "scenario3_countermeasures.csv",
        NUM_RUNS
    );

    // =========================================
    // SCENARIO 4: BOTH FACTORS
    // Congestion and immunization.
    // =========================================

    cfg.enable_congestion = true;

    cfg.enable_immunization = true;

    cfg.immunization_start = 0;

    cfg.seed = 45678;

    ScenarioSummary s4 = run_scenario(
        cfg,
        "Scenario 4 - Both Factors",
        "scenario4_both_factors.csv",
        NUM_RUNS
    );

    // =========================================
    // SCENARIO 5: DELAYED PATCHING
    // Both factors, immunization starts at t=100.
    // This is our project extension.
    // =========================================

    cfg.enable_congestion = true;

    cfg.enable_immunization = true;

    cfg.immunization_start = 100;

    cfg.seed = 56789;

    ScenarioSummary s5 = run_scenario(
        cfg,
        "Scenario 5 - Delayed Patching",
        "scenario5_delayed_patching.csv",
        NUM_RUNS
    );

    // =========================================
    // FINAL COMPARISON
    // =========================================

    cout << "\n========================================\n";
    cout << "           FINAL COMPARISON\n";
    cout << "========================================\n";

    cout << "\nScenario                  Mean Final Infected\n";

    cout << "1. Baseline               "
         << s1.mean_final_infected << "\n";

    cout << "2. Congestion             "
         << s2.mean_final_infected << "\n";

    cout << "3. Countermeasures        "
         << s3.mean_final_infected << "\n";

    cout << "4. Both factors           "
         << s4.mean_final_infected << "\n";

    cout << "5. Delayed patching       "
         << s5.mean_final_infected << "\n";

    cout << "\nAll five experiments completed successfully.\n";

    return 0;
}
