Turn-penalty-optimal Eulerian collection routes - MATLAB / GNU Octave code
===========================================================================

Companion code for: S. Naik, "Selecting among equally short collection routes:
turn-penalty-optimal Eulerian circuits via transition systems, with
capacity-feasible trip splitting".

No third-party toolbox is needed.

How to run
  MATLAB R2019b or later : cd to this folder, then   R = run_turnaware_study;
                           (exact solver: intlinprog, Optimization Toolbox)
  GNU Octave 8 or later  : same command (exact solver: built-in glpk).
                           For figures without a display:
                           xvfb-run -a octave --eval "R = run_turnaware_study;"
  R = run_turnaware_study(false) skips the random-grid benchmark (Part 10).

Outputs
  results/run_log.txt      complete console report (every number in the paper)
  results/summary.txt      one-page summary
  results/table_*.csv      streets, repairs, optimal circuit, sample, trips,
                           capacity sensitivity, benchmark
  results/results.mat      all results (struct R)
  figs/fig01..fig12        figures, PDF + PNG

Parts of the driver
   1  Ward network, degrees, odd vertices
   2  Closed parity repair: all 15 matchings, subset-DP minimum T-join
   3  Turn model, transition systems, local lower bound
   4  Hierholzer baseline, merge-and-switch heuristic, exact 0-1 programme,
      validation of the optimum
   5  Uniform random Eulerian circuits (rejection sampling), count estimate
   6  Joint analysis of all 15 repairs (deadhead vs optimal turn delay)
   7  Capacity-feasible trips (exact split), random-circuit shift distribution
   8  Turn-penalty sensitivity (30 exact re-solves)
   9  Payload sensitivity and split-aware local search
  10  Benchmark on 60 random street grids (exact vs heuristic vs baseline,
      branch-and-bound cross-check on small instances)

src/ folder
  ward_network.m        case-study ward
  apsp_floyd.m, sp_path.m
  min_matching_dp.m, all_matchings.m
  eulerize_closed.m
  ts_prepare.m          half-edges, turn angles, delays, vertex matchings, LB
  ts_apply.m, ts_trace.m, ts_eval.m, ts_mincost.m, ts_local_partner.m,
  ts_components.m       transition-system utilities
  ts_ilp.m, milp_min.m  exact programme with lazy connectivity cuts
  ts_heuristic.m        merge-and-switch heuristic
  ts_bnb.m              branch and bound (verification on small instances)
  hierholzer_partner.m  distance-only baseline circuit
  circuit_turns.m       per-transition class / delay / angle
  ts_sample.m, wh_rand.m  uniform sampling, Wichmann-Hill generator
  split_tour.m, split_ls.m  exact split, split-aware local search
  random_ward.m, solve_instance.m  benchmark instances and pipeline
  logmsg.m, make_figures.m

Reproducibility
  All randomness uses the Wichmann-Hill (AS 183) generator coded in
  src/wh_rand.m, so every number except wall-clock time is identical in
  MATLAB and Octave. Reported numbers: GNU Octave 8.4.0 with GLPK.
