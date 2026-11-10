# 3-Minute Presentation Script (180 Seconds)
## Project P20: Code Red Worm Propagation Modeling and Analysis
**Team Registration:** T05  
**Speakers:** Bhavashruthi M (Speaker 1) & D M Subashri (Speaker 2)  
**Accompanying Slide:** `ONE_SLIDE.html` (4 Quadrants)

---

### Timing & Transition Overview
* **0:00 - 0:45 (45s): Application & Cybersecurity Problem** — *Speaker 1 (Bhavashruthi)*
* **0:45 - 1:15 (30s): Mathematical & Epidemic Model** — *Speaker 1 (Bhavashruthi)*
* **1:15 - 1:55 (40s): Parallelism & OpenMP Monte Carlo Scaling** — *Speaker 2 (Subashri)*
* **1:55 - 2:40 (45s): Computational Bottleneck & Roofline Analysis** — *Speaker 2 (Subashri)*
* **2:40 - 3:00 (20s): Conclusion & Key Takeaway** — *Speaker 1 & 2 Joint*

---

### Second-by-Second Presentation Script

#### [0:00 - 0:45] Section 1: Application (Cybersecurity & Patch Latency)
**Speaker 1 (Bhavashruthi):**
> *"Good morning professors and everyone. I am Bhavashruthi, and together with Subashri, we present Project P20: **Code Red Worm Propagation Modeling and Analysis**.
> 
> [Point to Quadrant 1 on Slide]
> On July 19, 2001, the Code Red worm shocked the world by infecting over **359,000 Microsoft IIS web servers in less than 14 hours**. It did this by generating completely random 32-bit IPv4 addresses and probing port 80 across the globe.
> 
> In cybersecurity, the critical question is: **How long do network operators have before internet saturation becomes irreversible?**
> 
> In our research, we modeled the **Human Patch Delay $T$**. As shown in our propagation curve:
> * If administrators deploy patches within **2 hours**, peak infection is contained under **8%**.
> * But if patching is delayed past **8 hours**, infection explodes past **51% peak and 88% cumulative saturation**, overwhelming the global internet."*

---

#### [0:45 - 1:15] Section 2: Model (Stochastic Two-Factor Epidemic)
**Speaker 1 (Bhavashruthi):**
> *"Now, how do we model this mathematically?
> 
> [Point to Quadrant 2 on Slide]
> Traditional epidemic models like the classical SIR model fail because they assume a constant transmission rate $\beta$. 
> 
> Following the landmark paper by Zou, Gong, and Towsley, we implemented the **Two-Factor Worm Model**:
> 1. **Factor 1 is Delayed Human Response ($T$):** Hosts are only patched after reaction delay $T$ at removal rate $\gamma$.
> 2. **Factor 2 is Network Congestion:** As hundreds of thousands of infected servers scan simultaneously, router queues overflow and drop packets. We model this dynamic throttling as $\beta(t) = \frac{\beta_0}{1 + \alpha I(t)}$.
> 
> Because real worm outbreaks have discrete randomness—especially in the early seed phase—we use a **stochastic Monte Carlo agent model** rather than deterministic continuous equations."*

---

#### [1:15 - 1:55] Section 3: Parallelism (OpenMP & Monte Carlo Motif)
**Speaker 2 (Subashri):**
> *"Thank you, Bhavashruthi. I am Subashri, and I will walk you through the High Performance Computing implementation.
> 
> [Point to Quadrant 3 on Slide]
> Our project maps directly to **Monte Carlo Methods**, which is **Dwarf Number 4** of the 7 Berkeley HPC Motifs. 
> 
> To obtain statistically rigorous expected values and confidence intervals, we simulate **15,000 independent stochastic trajectories** using **OpenMP shared-memory multi-threading** on a 14-core Intel processor.
> 
> Key parallel optimizations:
> 1. **Zero Lock Contention:** We assigned an independent 64-bit XorShift pseudo-random number generator state to each thread, completely eliminating synchronization overhead.
> 2. **Dynamic Work Scheduling:** We used `#pragma omp for schedule(dynamic)` to prevent thread idling when trajectories terminate early.
> 3. **Scaling Results:** We achieved a **$3.22\times$ speedup on 4 threads with 80.4% efficiency**, and a peak speedup of **$4.65\times$ on 8 threads**, reducing runtime from 3.89 seconds down to 0.84 seconds."*

---

#### [1:55 - 2:40] Section 4: Bottleneck & Roofline Model Analysis
**Speaker 2 (Subashri):**
> *"Finally, we evaluated where our application hits hardware limits using Patterson's **Roofline Model**.
> 
> [Point to Quadrant 4 on Slide]
> Our processor has a theoretical peak double-precision compute performance of **780 GFLOP/s** and a DRAM bandwidth of **76.8 GB/s**. This defines an architectural **Ridge Point of $10.16\text{ FLOPs per Byte}$**.
> 
> We profiled the simulation:
> * Total arithmetic operations per step: **39 FLOPs**.
> * Memory traffic per step: **24 Bytes** for state reads and writes.
> * Measured Operational Intensity: **$I = 1.82\text{ FLOPs per Byte}$**.
> 
> Notice that $1.82$ is far to the left of the $10.16$ ridge point!
> Therefore, **our simulation is strictly MEMORY-BOUND**. 
> 
> The CPU floating-point units are starved while waiting for memory traffic. The crucial HPC insight is that optimizing cache blocking and data locality yields far greater speedup than vectorizing with AVX-512 instructions."*

---

#### [2:40 - 3:00] Conclusion & Wrap-Up
**Speaker 1 (Bhavashruthi):**
> *"In conclusion, our parallel Monte Carlo simulation confirms that internet worm defense has a strict deadline: patch within 4 hours, or face irreversible saturation."*

**Speaker 2 (Subashri):**
> *"And computationally, Roofline profiling proves that stochastic agent simulations are fundamentally memory-bandwidth bound. Thank you, and we are ready for your questions!"*

---

## Professor Viva / Defense Questions & Bulletproof Answers

### Q1: "Why is your simulation Memory-Bound and not Compute-Bound?"
* **Answer:** *"The operational intensity $I$ is only $1.82\text{ FLOPs/Byte}$, whereas the CPU ridge point is $10.16\text{ FLOPs/Byte}$. In Monte Carlo agent simulations, each step performs only 39 basic arithmetic operations (evaluating transition probabilities and generating random numbers), but accesses 24 bytes of memory (reading and updating host state variables). Because the ratio of math to memory traffic is very low, the execution pipelines spend most cycles stalled on memory latency and DRAM bandwidth rather than maxing out the ALUs."*

### Q2: "Why choose the Monte Carlo Motif instead of solving the differential equation directly with Runge-Kutta?"
* **Answer:** *"Differential equations provide an aggregate macroscopic mean, but they assume continuous homogeneous mixing. Real-world worm propagation is stochastic: an infected host generates random IP addresses, packets are randomly dropped, and in the early outbreak phase when only 10 hosts are infected, random fluctuations determine whether the worm explodes or dies out. Monte Carlo sampling over thousands of discrete stochastic paths gives us true variance, probability distributions, and worst-case percentiles that ODEs cannot capture."*

### Q3: "Why did your parallel efficiency drop from 97% on 2 threads to 58% on 8 threads?"
* **Answer:** *"There are two main reasons: First, Amdahl's Law serial fraction: thread creation, reduction aggregation, and thread pool orchestration take a non-zero constant time. Second, as we increase thread count on a memory-bound application, all threads compete for the shared L3 cache and the dual-channel DDR memory bus, creating memory bus saturation."*

### Q4: "How does your random number generation avoid false sharing and race conditions in OpenMP?"
* **Answer:** *"We strictly avoided standard C `rand()` or `drand48()`, which use internal shared global state and mutual exclusion locks that destroy parallel speedup. Instead, we implemented a register-based 64-bit XorShift PRNG where each thread is seeded with an independent unique seed derived from `omp_get_thread_num()`. Each thread mutates its own private register state with zero synchronization locks."*

### Q5: "What is the cybersecurity takeaway regarding the patch delay $T$?"
* **Answer:** *"The response curve is highly non-linear. Between $T = 0$ and $T = 4$ hours, peak infection stays below 18%. But between 6 and 8 hours, the worm passes the critical epidemic inflection point, and the peak infection jumps to 51.4%, with final saturation exceeding 88%. This gives security operations centers (SOCs) a concrete maximum response threshold ($T < 4\text{ hours}$) before network degradation becomes catastrophic."*
