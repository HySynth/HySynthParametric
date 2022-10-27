# HySynthParametric

*HySynthParametric* is a tool to automatically synthesize a mathematical model
for a dynamical system.
*HySynthParametric* receives data in the form of (a collection of) *time series*
as input and outputs a *hybrid automaton* model and a value epsilon with the
guarantee that the model contains trajectories that are epsilon-close to the
time-series data.

The theoretical background of the synthesis algorithm implemented in
*HySynthParametric* is explained in [1].



## Instructions to install HySynthParametric

In addition to the code in this repository, you only need a
[Julia compiler](https://julialang.org/downloads/) (we used v1.8.1).

Make sure `julia` is available in your `PATH`.
On Linux you can add the following to `~/.profile` and log out and in again:

```shell
PATH=$PATH:PATH_TO_JULIA/bin
```


## Instructions to run HySynthParametric

This section explains a quick start to *HySynthParametric*.
As example we look at the benchmark script to reproduce the results from [1].

* Open a terminal in the *HySynthParametric* folder.

* Run the Julia REPL via `julia --project=.`

* Instantiate the package via `import Pkg; Pkg.instantiate()`.

* Run the benchmark script via `include("examples/run_benchmarks.jl")`.

Note that the package and all its dependencies get precompiled the first time
they are loaded.
This may take a long time.

The outputs are written to the directory where Julia was run from.


## References

[1] Miriam García Soto, Thomas A. Henzinger, and Christian Schilling:
*Synthesis of parametric hybrid automata from time series*.
Proceedings of the 20th International Symposium on Automated Technology for
Verification and Analysis (ATVA) 2022.
[DOI](https://doi.org/10.1007/978-3-031-19992-9_22).
[PDF](https://arxiv.org/abs/2208.06383).

```bibtex
@inproceedings{GarciaHS22,
  author    = {Miriam Garc{\'{\i}}a Soto and
               Thomas A. Henzinger and
               Christian Schilling},
  editor    = {Ahmed Bouajjani and
               Luk{\'{a}}s Hol{\'{\i}}k and
               Zhilin Wu},
  title     = {Synthesis of parametric hybrid automata from time series},
  booktitle = {{ATVA}},
  series    = {LNCS},
  volume    = {13505},
  pages     = {337--353},
  publisher = {Springer},
  year      = {2022},
  url       = {https://doi.org/10.1007/978-3-031-19992-9_22},
  doi       = {10.1007/978-3-031-19992-9_22}
}
```
