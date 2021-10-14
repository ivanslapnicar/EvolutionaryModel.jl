# Universal Evolutionary Model for Periodical Organisms

Accompanying Julia code for the paper "Universal evolutionary model for periodical organisms"

```
Eric Goles, Ivan Slapničar, Marco A. Lardies, "Universal Evolutionary Model for Periodical Species", Complexity, vol. 2021, Article ID 2976351, 15 pages, 2021. https://doi.org/10.1155/2021/2976351
```

The manuscript can be found [here](https://arxiv.org/abs/2010.00940).

The code in the file `Simulation.jl` (also in the Pluto notebook `Simulation_p.jl`) was used to run all simulations and produce all figures in the paper.

## Running the simulations
There are three ways to run the simulations.

1. Open the file [Simulation_p.jl.html](https://ivanslapnicar.github.io/EvolutionaryModel.jl/Simulation_p.jl.html),  press `Edit or run this notebook` button and then choose `binder`. This will downlaod necessary packages and star the notebook (in few minutes).

2. Download the repository either as zip-file or clone it to your local file system by running
```
git clone git@github.com:ivanslapnicar/EvolutionaryModel.jl.git
```
At Julia prompt, run
`include("Simulation.jl")`

3. Alternatevely, after dowloading you can use the Pluto notebook `Simulation_p.jl`.

## Citation

If you use any of this code and wish to cite it, please cite the following publication
```
@article{GSL20,
author = {Eric Goles, Ivan Slapničar and Marco A. Lardies},
title = {Universal Evolutionary Model for Periodical Organisms},
journal = {Complexity},
volume = {2021},
year = {2020},
number = {2976351},
pages = {1-15},
url = {https://www.hindawi.com/journals/complexity/2021/2976351/},
doi = {10.1155/2021/2976351}
}
```
