# Code Red Worm Propagation Simulation Engine (5-Scenario Analysis)

[![HPC Motif](https://img.shields.io/badge/HPC%20Motif-Monte%20Carlo%20(Dwarf%20%234)-blue.svg)](#)
[![Parallel Framework](https://img.shields.io/badge/Parallelism-C%2B%2B%20%7C%20OpenMP-green.svg)](#)
[![Domain](https://img.shields.io/badge/Domain-Cybersecurity-red.svg)](#)

Academic High Performance Computing (HPC) research project modeling internet-scale stochastic worm propagation, network congestion throttling, and delayed patching response dynamics using parallel Monte Carlo simulations in C++ with OpenMP.

---

## 👥 Team Information (Team T05 / Project P20)
* **Bhavashruthi M** (Roll No: `24011103008`)
* **D M Subashri** (Roll No: `24011103010`)
* **Reference Publication:** Cliff C. Zou, Weibo Gong, Don Towsley, *"Code Red Worm Propagation Modeling and Analysis"*, ACM Conference on Computer and Communications Security (CCS), 2002.

---

## 🔬 The 5 Simulation Scenarios

| Scenario | Network Congestion | Human Countermeasures | Patch Start Time | Core Observation |
| :--- | :---: | :---: | :---: | :--- |
| **1. Baseline Epidemic** | Disabled ($D(t)=D_0$) | Disabled ($\gamma=0$) | N/A | Classical simple epidemic: unchecked exponential spread saturating 100% of hosts by $t=439$. |
| **2. Congestion Slowdown** | Enabled ($D(t)=D_0+Y(t)$) | Disabled ($\gamma=0$) | N/A | High probe volume causes ARP exhaustion and router drops, stretching saturation time by $4.6\times$ ($t=2000$). |
| **3. Human Countermeasures** | Disabled ($D(t)=D_0$) | Enabled ($\gamma=0.012$) | $t=0$ (Immediate) | Active immunization & cleaning removes active hosts, with $I(t)$ peaking at 219,512 and decaying to 0. |
| **4. Two-Factor Model** | Enabled ($D(t)=D_0+Y(t)$) | Enabled ($\gamma=0.012$) | $t=0$ (Immediate) | Both factors combined: matches Figure 7 of Zou et al. (CCS 2002), capping infections at 67.8% ($244,080$) with $I(t)$ decaying to 0. |
| **5. Delayed Patching** | Enabled ($D(t)=D_0+Y(t)$) | Enabled ($\gamma=0.012$) | $t=250$ (Delayed) | Modeling delay response: unhindered initial growth pushes infection past the tipping point to 88.4% ($318,298$). |

---

## 🚀 Compilation and Execution

### Using GCC / G++ (Linux / HPC Lab):
```bash
# Compile with OpenMP and O3 optimization
g++ -O3 -Wall -std=c++17 -fopenmp worm_simulation.cpp -o worm_simulation

# Run 100 Monte Carlo trials per scenario
./worm_simulation 100
```

### Using Makefile:
```bash
make
make run
```

---

## 📊 Interactive 5-Scenario Web Dashboard

Open [`dashboard.html`](dashboard.html) directly in any modern web browser (no web server required).

### Key Dashboard Features:
1. **Interactive Multi-Scenario S-Curves ($J(t)$):** Real-time toggling of any of the 5 scenarios with distinct color-coded trajectories.
2. **Active Infectious Dynamics ($I(t)$):** Visualizing infection wave heights, peak times, and eradication phases.
3. **Time-Step Inspector:** Interactive slider to scrub across $t = 0\dots 2000$ with live crosshair readings for all scenarios.
4. **Delayed Patching Analysis Panel:** Direct comparative evaluation of Scenario 4 (Prompt, $t=0$) versus Scenario 5 (Delayed, $t=250$).
5. **HPC & OpenMP Scaling Metrics:** Complete summary table including peak infection counts, time-to-peak, total immunized hosts, and wall-clock execution times.

---

## 📁 Repository Deliverables
* [`worm_simulation.cpp`](worm_simulation.cpp): Pure C++ OpenMP engine implementing all 5 scenarios with thread-safe SplitMix64 PRNG.
* [`dashboard.html`](dashboard.html): Standalone HTML5/Canvas interactive visualization dashboard.
* [`REPORT.md`](REPORT.md): Formal HPC academic research paper report.
* [`PRESENTATION_SCRIPT.md`](PRESENTATION_SCRIPT.md): 3-minute oral presentation defense script with anticipated professor viva questions.
* [`scenario*.csv`](scenario1_baseline.csv): Full averaged time-series datasets for all 5 scenarios.
* [`simulation_summary.csv`](simulation_summary.csv): Aggregate metrics table across all 5 scenarios.

