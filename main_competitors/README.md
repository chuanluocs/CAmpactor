# Main Competitors of *CAmpactor*

The implementations of the three PCAG solvers are in the same directory as this document. The directories `AutoCCAG/`, `FastCA/` and `TCA/` contains the implementation of *AutoCCAG*, *FastCA* and *TCA*, respectively. 

## Instructions for Building the Solvers

```
sh build.sh
```

Note that the three PCAG solvers should be built on a 64-bit GNU/Linux operating system. 

## Instructions for Running the Solvers

We provide `simple_run.py`, a Python3 script, as a unified interface for running these PCAG solvers. It allows the user to run *AutoCCAG*, *FastCA* or *TCA* with an initial PCA, which is either automatically generated by *SamplingCA* or explicitly provided by the user. For its detailed usage, users may execute it with the `-h` flag on. 

Note that the initial PCA should be in the same format as the PCA constructed by *SamplingCA*. 

**Tip**: The command for running *SamplingCA* is embedded in this script. Once there is unexpected error in running *SamplingCA*, users may refer to [SamplingCA/README.md](../SamplingCA/README.md) for solutions. 

## Example Command for Running the Solvers

An example of using `simple_run.py` to run *FastCA*: 

```
python3 simple_run.py -s F -i ../cnf_benchmarks/linux.cnf -o linux_PCA_by_FastCA.out --seed 1 --cutoff 20
```

The command above first calls *SamplingCA* to solve the instance `../cnf_benchmarks/linux.cnf` with default hyper-parameter settings, since the initial PCA is not provided by the user. Then *FastCA* is called to optimize the initial PCA and the cutoff time is set to 20 seconds. The result is stored in `linux_PCA_by_FastCA.out`. Here both *SamplingCA* and *FastCA* use the random seed of 1. 

The console output is expected to be similar with the following text:

```
c running SamplingCA ...
...
...  <-- output of SamplingCA
...
c running FastCA ...
1.54254 110     0
3.49662 109     41
3.63432 108     46
4.08852 107     54
5.54425 106     110
8.32681 105     169
8.44476 104     173
13.2817 103     278
13.3068 102     279
c done. 
```

The output between the lines `c running FastCA ...` and `c done.` is the output of *FastCA*. The second column represents the progress of PCA optimization, its size decreasing from 110 to 102. According to the console output, we expect the size of the PCA in `linux_PCA_by_FastCA.out` to be exactly 102. 

The command is almost the same for running *AutoCCAG* and *TCA*, except that the argument `-s` (or its full name `--solver`) should take `A` for *AutoCCAG* and `T` for *TCA*. And the console output of *AutoCCAG* and *TCA* is also in similar format with *FastCA*. 

## Related Experimental Results

We list the experimental results, which are related to the main competitors in this directory, as follows. 
- [Results_of_CAmpactor_and_its_SOTA_competitors.csv](../experimental_results/Results_of_CAmpactor_and_its_SOTA_competitors.csv): Results of *CAmpactor* (with *SamplingCA* as its initialization algorithm) and its state-of-the-art competitors on all testing instances. 
- [Results_of_CAmpactor_on_its_generality.csv](../experimental_results/Results_of_CAmpactor_on_its_generality.csv): Results of the generality of *CAmpactor* (with *AutoCCAG*, *FastCA* and *TCA* as its initialization algorithms) on all testing instances.
- [Results_of_two_versions_of_AutoCCAG_FastCA_TCA.csv](../experimental_results/Results_of_two_versions_of_AutoCCAG_FastCA_TCA.csv): Results of the versions of *AutoCCAG*, *FastCA* and *TCA*, which employ *SamplingCA* as their initialization algorithms, and the original versions of *AutoCCAG*, *FastCA* and *TCA*, which use their own, original initialization algorithms.