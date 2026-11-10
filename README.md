# Minimal Compound Combinations of Traditional Chinese Medicine Formulas

## Description

Code and data for the paper "Multi-Objective Evolutionary Algorithms and Exact Epsilon-Constraint Optimization for Minimal Compound Combinations of Traditional Chinese Medicine Formulas".

The selection of a compound combination from a formula is formulated as a three-objective sparse 0-1 problem: the efficacy `f1` (negative, bioavailability-weighted coverage of COPD targets plus a reward for complementary compound pairs), the number of compounds `f2`, and the toxicity `f3` (worst single-compound toxicity plus CYP co-inhibition of compound pairs), with 1 to 10 compounds per combination. The weight `lam` of the pairwise terms takes the values 0, 0.5, and 1.

## Dataset information

All input data are third-party and public; the derived problem instances are included in `data/`.

- Compounds of the herbs and compound-target activities: NPASS, 2026 update (https://bidd.group/NPASS/), and CMAUP, 2024 update (https://bidd.group/CMAUP/). A compound-target pair is known when an activity value of at most 10 uM is reported.
- COPD targets and their weights: Open Targets Platform, release 26.09 (https://platform.opentargets.org/), disease MONDO_0005002, clinical and animal-model evidence.
- Oral bioavailability, toxicity, and CYP inhibition: predictions of ADMET-AI 2.0.1 (https://github.com/swansonk14/admet_ai).
- Network distances: the human interactome of Cheng et al. (2019), Supplementary Data 1 of https://doi.org/10.1038/s41467-019-09186-x.
- Herbs of the formulas and their botanical sources: Chinese Pharmacopoeia (2020 edition).

The terms of use of each source apply to the data derived from it.

## Requirements

- Python 3.9 or later with NumPy.
- PySCIPOpt with SCIP (only for the exact method, `scripts/run_exact.py`). The experiments used SCIP 10.0 and PySCIPOpt 6.3.0.

```
pip install -r requirements.txt
```

## Code information

```
tcmsel/
  problem.py        problem definition and instance loader
  exact.py          augmented epsilon-constraint method with SCIP
  core.py           budgeted evaluation and population container
  operators.py      sorting, selection, and variation operators
  algorithms/       SPEA2, NSGA-II, NSGA-III, SparseEA, MOEA/PSL, PM-MOEA, random search
  greedy.py         degree ranking and cumulative contribution
  indicators.py     non-dominated filter, exact 3-D hypervolume, IGD+
  runner.py         one MOEA run with greedy initial solutions
scripts/
  run_moea.py       run one algorithm on one instance
  run_exact.py      compute the Pareto front with the exact method
  run_greedy.py     evaluate the two greedy rules
  evaluate.py       relative HV, HV ratio, and IGD+ of saved runs
data/
  instances/        27 problem instances
  exact/            points returned by the exact method and a summary
  reference_sets/   reference sets used for the indicators
  statistics/       results of the statistical tests reported in the paper
```

### Data files

**Instances** (`data/instances/<name>.npz`). `F1`, `F2`, `F3`, `F4`, `F6`, and `F8` are the candidate pools of Bufei, Xiaoqinglong, Suzi Jiangqi, Dingchuan, Maxing Shigan, and Erchen decoctions; `ALL` is their merged pool; `ALL_n40_r0` to `ALL_n40_r9` and `ALL_n60_r0` to `ALL_n60_r9` are random sub-pools of 40 and 60 compounds. Each file contains

| Key | Shape | Meaning |
|---|---|---|
| `ik14` | N | first 14 characters of the standard InChIKey of each compound |
| `targets` | T | Ensembl gene identifiers of the coverable COPD targets |
| `B` | N x T | 1 if the compound and the target form a known pair (activity of at most 10 uM) |
| `w` | T | target weights |
| `e` | N | exposure weights (predicted probability of an oral bioavailability above 20%) |
| `tox` | N | toxicity of each compound |
| `CE` | N x N | 1 if two compounds are complementary in the interactome |
| `DDI` | N x N | CYP co-inhibition of two compounds |
| `Z_ce`, `Z_ddi` | scalar | normalizing constants of the pairwise terms (45) |
| `z` | N | z-score of the closest distance of each compound to the disease module |

**Exact method** (`data/exact/<name>_lam<lam>.csv`). Every point returned for the size limits k = 1 to 10, with its objective values, relative gap, solver status, and the indices of the selected compounds. `summary.csv` states for each instance and `lam` whether the front is complete, the largest size limit `k_c` up to which it is complete, and the hypervolume of the exact points relative to the reference set.

**Reference sets** (`data/reference_sets/<name>_lam<lam>.npz`). Non-dominated points `F`, their combinations `X`, the ideal and nadir points `lo` and `hi` used for normalization, and a flag for points found by the exact method.

**Statistical tests** (`data/statistics/`). `pairwise_tests.csv` compares the HV of every algorithm with that of the base algorithm, the MOEA with the highest median relative HV after 20,000 evaluations (SPEA2), on each instance and budget: medians, the Vargha-Delaney effect size `A12`, the p-value of the two-sided Wilcoxon rank-sum test, the Holm-adjusted p-value over the six comparisons of an instance, and the verdict (`+`/`-`: significantly higher/lower HV than the base algorithm, `=`: no significant difference). `friedman_tests.csv` gives the Friedman and Iman-Davenport statistics, the Nemenyi critical difference, and the mean ranks for each budget. The instance set `lambda>0` contains the six formulas with `lam` = 0.5 and 1, and `lambda=0` the six formulas with `lam` = 0.

## Usage instructions

Run SparseEA once on Maxing Shigan decoction with `lam = 1` and 20,000 evaluations, and score the run:

```
python scripts/run_moea.py --instance F6 --lam 1 --alg SparseEA --budget 20000 --run 0
python scripts/evaluate.py results/moea/F6/lam1/SparseEA/B20000
```

Algorithms: `SPEA2`, `NSGA-II`, `NSGA-III`, `SparseEA`, `MOEA/PSL`, `PM-MOEA`, `Random`. Budgets: 2000, 5000, 10,000, and 20,000 evaluations. The experiments used runs 0 to 29 for each algorithm, instance, `lam`, and budget; the run index and the other arguments determine the random seed, so the reported runs are reproduced exactly.

Compute the exact front of an instance (ten processes, one per size limit):

```
python scripts/run_exact.py --instance F6 --lam 0
python scripts/run_exact.py --instance F6 --lam 1 --time-per-model 3600 --time-per-k 14400
```

Evaluate the greedy rules of network pharmacology:

```
python scripts/run_greedy.py --instance F4 --lam 1
```

The option `--no-exposure` of `run_moea.py` sets all exposure weights to 1, and `--no-greedy` starts the algorithms without greedy initial solutions.

## Methodology

1. Problem. For a candidate pool of N compounds, a combination is a 0-1 vector with 1 to K = 10 selected compounds. Three objectives are minimized: `f1`, the negative of the sum of the bioavailability-weighted coverage of the weighted COPD targets and `lam` times the number of complementary compound pairs divided by 45; `f2`, the number of compounds; `f3`, the largest single-compound toxicity plus `lam` times the CYP co-inhibition of the compound pairs divided by 45.
2. Exact method. The augmented epsilon-constraint method (AUGMECON2) solves one mixed-integer model per size limit k and toxicity bound with SCIP; each model has a relative optimality gap of 1e-4 and a limit of one hour, and each size limit a limit of four hours; the front is complete when every model reaches this gap within the limits.
3. Evolutionary algorithms. SPEA2, NSGA-II, NSGA-III, SparseEA, MOEA/PSL, PM-MOEA, and random search with a population of 100 (91 for NSGA-III), binary coding, and constraint handling that ranks infeasible solutions behind feasible ones. The initial population contains greedy solutions whose evaluations, at most 10% of the budget, are charged to the budget. Budgets of 2000, 5000, 10,000, and 20,000 evaluations, 30 runs each.
4. Indicators. Objectives are normalized with the ideal and nadir points of the reference set; the hypervolume (HV) is computed exactly with the reference point (1 + 1/12)(1, 1, 1). Reported are the HV relative to the reference set, the HV ratio to a complete exact front, and IGD+.
5. Statistics. Two-sided Wilcoxon rank-sum tests against the MOEA with the highest median relative HV, Holm correction over the comparisons of an instance, the Vargha-Delaney effect size A12, and the Friedman test with the Iman-Davenport correction and the Nemenyi critical difference over instances.

## Citations

If you use this code or data, please cite the article above and the methods and data sources it builds on, in particular:

- Mavrotas G, Florios K. An improved version of the augmented epsilon-constraint method (AUGMECON2) for finding the exact Pareto set in multi-objective integer programming problems. Applied Mathematics and Computation 2013;219:9652-9669. https://doi.org/10.1016/j.amc.2013.03.002
- Achterberg T. SCIP: solving constraint integer programs. Mathematical Programming Computation 2009;1:1-41. https://doi.org/10.1007/s12532-008-0001-1
- Tian Y, Cheng R, Zhang X, Jin Y. PlatEMO: a MATLAB platform for evolutionary multi-objective optimization. IEEE Computational Intelligence Magazine 2017;12:73-87. https://doi.org/10.1109/MCI.2017.2742868
- Cheng F, Kovacs IA, Barabasi AL. Network-based prediction of drug combinations. Nature Communications 2019;10:1197. https://doi.org/10.1038/s41467-019-09186-x
- Swanson K, et al. ADMET-AI: a machine learning ADMET platform for evaluation of large-scale chemical libraries. Bioinformatics 2024;40:btae416. https://doi.org/10.1093/bioinformatics/btae416

## License and contribution guidelines

The code and data are released under the MIT License (see `LICENSE`). Data derived from third-party sources remain subject to the terms of those sources.

Bug reports and suggestions are welcome through the issue tracker of the repository. Contributions should be submitted as pull requests with a short description of the change.
