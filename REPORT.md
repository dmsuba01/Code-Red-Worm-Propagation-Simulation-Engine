# High Performance Computing Project Report
## P20 | Code Red Worm Propagation Modeling and Analysis
**Team Registration:** T05  
**Students:**  
* Bhavashruthi M (Roll No: 24011103008)  
* D M Subashri (Roll No: 24011103010)  

**Berkeley HPC Dwarf / Motif:** Monte Carlo Methods  
**Parallel Programming Framework:** OpenMP (Shared-Memory Multi-threading)  
**Target Domain:** Cybersecurity / Internet Infrastructure Resilience  
**Reference Paper:** Cliff C. Zou, Weibo Gong, Don Towsley, *"Code Red Worm Propagation Modeling and Analysis"*, ACM Conference on Computer and Communications Security (CCS), 2002.

---

## Executive Summary

On July 19, 2001, the **Code Red v2** worm infected over **359,000 Microsoft IIS web servers** across the global Internet in less than 14 hours. Standard biological epidemic models (such as the classic Kermack-McKendrick epidemic model) failed to predict the real-world trajectory of Code Red because they assumed static infection rates and ignored human intervention and network congestion. In this project, we implement and evaluate the **Two-Factor Worm Model** proposed by Zou, Gong, and Towsley using **parallel stochastic Monte Carlo simulations accelerated by OpenMP**. 

Our primary scientific objectives are:
1. To quantify how internet-scale random scanning saturates a vulnerable population.
2. To determine the critical tipping point where delayed human patching transforms an outbreak from a manageable incident into irreversible global infrastructure saturation.
3. To evaluate the parallel scaling characteristics and computational bottlenecks using the **Roofline Performance Model**.

---

## 1. Application: Cybersecurity & Internet Infrastructure

### 1.1 The Threat: Fast-Spreading Random Scanning Worms
A computer worm is an autonomous, self-propagating program that exploits software vulnerabilities to spread across computer networks without human interaction. Unlike email viruses, internet worms require no user action (such as opening an attachment).

**Code Red v2 Details:**
* **Vulnerability:** Buffer overflow in Microsoft IIS `idq.dll` index server ISAPI extension.
* **Scanning Mechanism:** Each infected server spawns 100 worker threads. Each thread generates pseudo-random 32-bit IPv4 addresses uniformly over $\Omega = 2^{32} \approx 4.29 \times 10^9$ possible addresses and sends an HTTP `GET /default.ida?NNNN...` probe on TCP port 80.
* **Global Impact:** In less than 14 hours, it saturated 359,000+ servers worldwide and generated severe packet storms that degraded global BGP routing tables.

```
+-----------------------------------------------------------------------------------+
|                            INTERNET IPV4 ADDRESS SPACE                            |
|                            (Omega = 2^32 = 4.29 Billion IPs)                      |
|                                                                                   |
|    [Infected Host A]                                                              |
|          |----(Random Scan #1: 198.51.100.42)-----> [Linux Server - Ignored]      |
|          |----(Random Scan #2: 203.0.113.88)------> [Empty IP - Timeout]          |
|          |----(Random Scan #3: 192.0.2.15)--------> [Vulnerable IIS] -> INFECTED! |
+-----------------------------------------------------------------------------------+
```

### 1.2 The Core Cybersecurity Research Problem
How long do national security teams and network operators have before the infection becomes uncontainable? 
We model the **Patch Response Delay ($T$)**:
* If system administrators patch vulnerable servers immediately ($T = 0\text{ to }2$ hours), the worm dies out with less than 8% total population compromise.
* If response is delayed beyond a critical tipping point ($T \ge 6\text{ to }8$ hours), infection breaches 50%, causing runaway internet-wide saturation (>90%) and triggering severe infrastructure collapse.

---

## 2. Model: From Classical SIR to the Zou-Gong-Towsley Two-Factor Model

### 2.1 The Classic Epidemic Model & Why It Fails
In classical epidemiology, the spread of an infectious agent is governed by the logistic differential equation:
$$\frac{dI(t)}{dt} = \beta \cdot I(t) \cdot (N - I(t))$$
Where:
* $N$ = Total vulnerable population (estimated at $\approx 360,000$ for Code Red).
* $I(t)$ = Number of infected hosts at time $t$.
* $S(t) = N - I(t)$ = Number of susceptible hosts.
* $\beta$ = Pairwise infection transmission rate: $\beta = \frac{\eta}{\Omega}$.
* $\eta$ = Average scan rate per host ($\approx 358\text{ scans/minute}$).
* $\Omega$ = Total address space ($2^{32} = 4,294,967,296$).

**Failure of the Simple Model:** The classic model predicts that the worm continues accelerating until nearly 100% of hosts are infected. Real-world telemetry collected by CAIDA (Cooperative Association for Internet Data Analysis) showed that Code Red's propagation slowed down significantly after infecting ~180,000 hosts.

### 2.2 The Two-Factor Worm Model Formulation
Zou, Gong, and Towsley (2002) introduced two real-world operational factors that fundamentally alter the epidemic trajectory:

#### Factor 1: Dynamic Human Countermeasures with Patch Delay ($T$)
System administrators do not respond instantaneously. There is a detection and reaction delay $T$:
* For $t < T$: Countermeasure rate $\gamma(t) = 0$.
* For $t \ge T$: Emergency advisories are published. Susceptible hosts are patched proactively ($S \rightarrow R$), and infected hosts are cleaned, patched, or quarantined ($I \rightarrow R$) at removal rate $\gamma$.

#### Factor 2: Network & Host Congestion Slowdown
As hundreds of thousands of infected machines simultaneously emit hundreds of scans per second, the aggregate traffic exceeds the switching capacity of core routers and local gateways. Packet drops, connection timeouts, and server crashes throttle the effective scan rate:
$$\beta(t) = \beta_0 \cdot \frac{1}{1 + \alpha \cdot I(t)}$$
Where $\alpha$ is the network congestion sensitivity coefficient and $\beta_0 = \frac{\eta}{\Omega}$.

#### Discrete Stochastic Monte Carlo Transition Equations
In our Monte Carlo simulation, we discretize time into steps $\Delta t = 1.0\text{ minute}$. At each step:
1. **Effective Infection Probability:**
   $$\lambda_{inf}(t) = \beta(t) \cdot I(t) \cdot \Delta t$$
   $$P_{inf}(t) = 1 - e^{-\lambda_{inf}(t)}$$
2. **Stochastic State Transitions:**
   $$\Delta I_{new} \sim \text{Binomial}(S(t), P_{inf}(t)) \approx \mathcal{N}(\mu_{inf}, \sigma_{inf}^2)$$
   $$\Delta R_S \sim \text{Binomial}(S(t), P_{patch, S}(t))$$
   $$\Delta R_I \sim \text{Binomial}(I(t), P_{patch, I}(t))$$
3. **State Conservation Updates:**
   $$S(t + \Delta t) = S(t) - \Delta I_{new} - \Delta R_S$$
   $$I(t + \Delta t) = I(t) + \Delta I_{new} - \Delta R_I$$
   $$R(t + \Delta t) = R(t) + \Delta R_S + \Delta R_I$$

---

## 3. Parallelism: The Monte Carlo Motif & OpenMP Multi-Threading

### 3.1 HPC Dwarf / Motif: Monte Carlo Methods
Monte Carlo methods rely on repeated random sampling to compute numerical results for complex stochastic physical systems. In epidemic propagation, microscopic randomness (which IP address is chosen, whether a packet is dropped, which admin patches first) creates divergent outbreak paths, especially in the early seed phase ($I_0 = 10$).

To compute expected value trajectories, variance, and 95% confidence intervals, we execute thousands of independent Monte Carlo trials.

### 3.2 OpenMP Parallelization Strategy
We employ **trial-level domain decomposition** using OpenMP shared-memory multi-threading:

```c
#pragma omp parallel num_threads(cfg->num_threads) \
    reduction(+:total_peak_infected, total_time_to_peak, total_flops, total_bytes)
{
    int tid = omp_get_thread_num();
    // Unique 64-bit seed per thread prevents false sharing and PRNG correlation
    uint64_t thread_seed = 0x853c49e6748fea9bULL ^ ((uint64_t)tid * 0x9e3779b97f4a7c15ULL);

    #pragma omp for schedule(dynamic, 1)
    for (int trial = 0; trial < num_trials; trial++) {
        uint64_t trial_seed = thread_seed + (uint64_t)trial * 1337ULL;
        SimMetrics res = run_single_stochastic_trial(cfg, trial_seed, NULL, total_steps);
        // Accumulate reductions
        total_peak_infected += res.peak_infected;
        total_flops         += res.total_flops;
        total_bytes         += res.total_bytes_accessed;
    }
}
```

#### Parallel Optimization Highlights:
1. **Thread-Safe PRNG (XorShift64):** Replaces POSIX `rand()` or `drand48()` (which have internal global mutex locks) with zero-lock 64-bit XorShift registers, eliminating thread synchronization bottlenecks.
2. **Dynamic Work Scheduling (`schedule(dynamic, 1)`):** Because some trials terminate early (if the worm is eradicated early by rapid patching) while others simulate the full 24-hour horizon, dynamic scheduling provides load balancing across threads.
3. **OpenMP Reductions:** Prevents race conditions when aggregating global statistics across multi-core architectures without explicit mutual exclusion (`#pragma omp critical`).

---

## 4. Bottleneck & Roofline Model Analysis

### 4.1 Theoretical Roofline Model Formulation
The **Roofline Model** (Williams, Waterman, Patterson, 2009) bounds the maximum attainable floating-point performance ($P$) of an application based on hardware parameters:
$$P = \min\left(P_{peak}, \; B_{peak} \times I\right)$$
Where:
* $P_{peak}$ = Theoretical Peak Double-Precision Compute Performance ($\text{GFLOP/s}$).
* $B_{peak}$ = Peak DRAM Memory Bandwidth ($\text{GB/s}$).
* $I$ = Operational Intensity (Arithmetic Intensity), defined as:
  $$I = \frac{\text{Total Arithmetic Floating-Point Operations (FLOPs)}}{\text{Total Memory Traffic to/from DRAM (Bytes)}}$$

The **Ridge Point (Knee)** represents the minimum operational intensity required to achieve peak computational performance:
$$I_{knee} = \frac{P_{peak}}{B_{peak}}$$

### 4.2 Hardware Specifications & Measured Parameters

| Parameter | Theoretical / Measured Value |
| :--- | :--- |
| **Processor Architecture** | 14-Core Intel Core Processor |
| **Peak DP Floating-Point Performance ($P_{peak}$)** | **780.0 GFLOP/s** |
| **Peak DRAM Bandwidth ($B_{peak}$)** | **76.8 GB/s** |
| **Roofline Ridge Point ($I_{knee}$)** | **10.15 FLOPs/Byte** |
| **Simulation Arithmetic Operations per Step** | ~39 FLOPs |
| **Memory Traffic per Step** | 24 Bytes ($S, I, R$ read & write) |
| **Measured Operational Intensity ($I_{sim}$)** | **0.95 to 1.62 FLOPs/Byte** |

### 4.3 Bottleneck Classification: Memory-Bound
Because our measured operational intensity ($I_{sim} \approx 1.62\text{ FLOPs/Byte}$) is substantially smaller than the architectural ridge point ($I_{knee} = 10.15\text{ FLOPs/Byte}$):
$$I_{sim} \ll I_{knee}$$

$$\therefore \textbf{The application is strictly MEMORY-BANDWIDTH & CACHE-LATENCY BOUND.}$$

#### Architectural Explanation:
1. In Monte Carlo stochastic simulations, the CPU spends the vast majority of time fetching agent state vectors from memory, computing pseudo-random state indices, and writing updated states back to RAM.
2. The number of floating-point operations per memory access is low (evaluating 2-3 exponential/trigonometric expressions and basic multiplications per state update).
3. The arithmetic execution units (ALUs / FPUs / AVX-512) sit idle while awaiting memory operands. As a result, the application's performance is limited by the DRAM bus and cache hierarchy rather than ALU throughput.

---

## 5. Experimental Results & Discussion

### 5.1 Experiment 1: Patch Latency Sensitivity Results

| Patch Delay $T$ (Hours) | Peak Infected ($I_{max}$) | Peak Saturation (%) | Time to Peak (Hours) | Final Patched Population (%) |
| :---: | :---: | :---: | :---: | :---: |
| **$T = 0\text{h}$** (Instant) | 6,283 | **1.7%** | 10.48h | 99.8% |
| **$T = 1\text{h}$** | 14,441 | **4.0%** | 10.01h | 99.5% |
| **$T = 2\text{h}$** | 28,292 | **7.9%** | 9.10h | 99.1% |
| **$T = 4\text{h}$** | 61,699 | **17.1%** | 8.02h | 97.8% |
| **$T = 6\text{h}$** | 106,378 | **29.5%** | 7.62h | 95.2% |
| **$T = 8\text{h}$ (Tipping Point)**| **185,209** | **51.4%** | **8.00h** | **88.6%** |
| **$T = 12\text{h}$** | 298,824 | **83.0%** | 12.00h | 64.2% |
| **$T = 16\text{h}$** | 339,751 | **94.4%** | 16.00h | 38.1% |
| **$T = 24\text{h}$** (No Patch) | 357,836 | **99.4%** | 23.98h | 0.6% |

**Key Finding:** There is a distinct non-linear tipping point at **$T \approx 6 - 8$ hours**. If security teams deploy patches within 4 hours, infection is limited to under 18%. Beyond 8 hours, the worm passes the 50% inflection point and infects >83% of the world's servers.

### 5.2 Experiment 2: Multi-Core Parallel Scaling & Speedup (15,000 Monte Carlo Trials)

| Thread Count ($p$) | Wall-Clock Time (s) | Measured Speedup ($S_p = T_1 / T_p$) | Ideal Speedup | Parallel Efficiency ($E_p = S_p / p$) | Attained Throughput |
| :---: | :---: | :---: | :---: | :---: | :---: |
| **1 Thread** | 3.889 s | **1.00x** | 1.00x | **100.0%** | 3,857 trials/s |
| **2 Threads** | 2.003 s | **1.94x** | 2.00x | **97.1%** | 7,489 trials/s |
| **4 Threads** | 1.210 s | **3.22x** | 4.00x | **80.4%** | 12,397 trials/s |
| **8 Threads** | 0.837 s | **4.65x** | 8.00x | **58.1%** | 17,921 trials/s |
| **14 Threads** | 0.892 s | **4.36x** | 14.00x | **31.2%** | 16,816 trials/s |

**Scaling Discussion:**
1. **Near-linear speedup up to 4 threads ($3.22\times$ speedup at 80.4% efficiency):** Demonstrates excellent Monte Carlo domain independence without lock contention.
2. **Sub-linear scaling beyond 8 threads:** As thread count approaches the physical core limit (Intel hybrid performance/efficient architecture), memory bus contention and L3 cache sharing begin to throttle throughput, consistent with our Roofline memory-bound diagnosis.

### 5.3 Experiment 3: 5-Scenario Comparative Analysis (Zou et al. CCS 2002 Framework)

To systematically isolate the individual and combined effects of network congestion and dynamic countermeasures, we evaluated 5 discrete scenarios matching the Zou et al. (CCS 2002) model specification across a population of $N = 360,000$ vulnerable IIS servers with continuous disinfection ($\frac{dR}{dt} = \gamma I(t)$, Section 6, Eq. 14):

| Scenario | Network Congestion | Human Countermeasures | Patch Delay ($T$) | Peak Active ($I_{max}$) | Peak Time ($t$) | Cumulative Final ($J$) | Saturation (%) | Total Patched ($R$) |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **1. Baseline Epidemic** | Disabled ($D_0=10$) | Disabled ($\gamma=0$) | N/A | **360,000** | $t = 439$ | 360,000 | **100.0%** | 0 |
| **2. Congestion Slowdown** | Enabled ($D_0 + Y$) | Disabled ($\gamma=0$) | N/A | **359,999** | $t = 2000$ | 359,999 | **100.0%** | 0 |
| **3. Human Countermeasures** | Disabled ($D_0=10$) | Enabled ($\gamma=0.012$) | $t = 0$ | **219,512** | $t = 147$ | 346,834 | **96.3%** | 360,000 |
| **4. Two-Factor Worm Model** | Enabled ($D_0 + Y$) | Enabled ($\gamma=0.012$) | $t = 0$ | **105,336** | $t = 155$ | 244,080 | **67.8%** | 357,762 |
| **5. Delayed Patching Response** | Enabled ($D_0 + Y$) | Enabled ($\gamma=0.012$) | $t = 250$ | **291,556** | $t = 249$ | 318,298 | **88.4%** | 359,643 |

#### Key Analytical Takeaways (Matching Zou et al. Figure 7):
1. **Active Infectious Waves Decay to Zero (Figure 7 Alignment):** In Scenarios 3, 4, and 5, active infectious hosts $I(t)$ do not remain permanently elevated. After reaching peak infection, continuous human disinfection ($\frac{dR}{dt} = \gamma I(t)$) outpaces the dwindling pool of susceptible hosts, causing $I(t)$ to decline smoothly back toward zero while cumulative infected $J(t) = I(t) + R(t)$ plateaus.
2. **Network Congestion as a Temporal Dampener:** Comparing Scenario 1 vs Scenario 2 shows that router ARP cache exhaustion and queue drops do not prevent complete saturation on their own, but they delay peak infection by $4.6\times$ (from $t = 439$ to $t = 2000$), buying precious time for defensive intervention.
3. **Synergy in the Two-Factor Model (Scenario 4):** Combining router congestion with prompt human countermeasures flattens peak active infections by **70.7%** (from 360,000 down to 105,336) and limits total worldwide compromise to **67.8%** ($244,080$ servers), closely replicating the observed July 2001 Code Red trajectory.
4. **Delayed Patching Response Impact (Scenario 5):** When human countermeasures are delayed until $t = 250$ steps, the worm exploits the uninhibited window to compromise an additional **+74,218 servers** (+20.6% higher saturation), proving the existence of a sharp threshold for rapid incident response.

---

## 6. Conclusion
By pairing the **Zou-Gong-Towsley Two-Factor Worm Model** with an **OpenMP-parallelized Monte Carlo simulator**, we accurately recreated the July 2001 Code Red outbreak dynamics and quantified the critical human response window. Computational profiling via the **Roofline Model** confirmed that stochastic network simulations are memory-bound ($I \approx 1.6\text{ FLOP/B} \ll 10.15\text{ FLOP/B}$), underscoring that memory hierarchy optimizations (cache blocking, register PRNG states) are far more impactful than compute vectorization for this class of HPC problems.
