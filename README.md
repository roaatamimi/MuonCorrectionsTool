# MuonCorrections

## Description

This page contains instructions for using the Rochester corrections for muon momentum scale and resolution. The corrections are used to compensate for the misalignments in data and Monte Carlo that the CMS reconstruction software does not fully correct. The misalignments for data and MC are different and corrections have been extracted for both.

These instructions are for the Run1 corrections and a [2012 dataset](http://opendata.cern.ch/record/12341) is used as an example. `RochesterCorrections` contains the official code for the corrections. `Test` contains `Analysis.C` which reads the dataset, applies the corrections to a muon pair, computes invariant mass and produces an output-file with the corrected data. `Plot.C` creates the plot below which shows that the corrections have been applied correctly. `Plot.C` is based on [this](https://cms-opendata-workshop.github.io/workshop-lesson-tagandprobe/index.html) Tag and Probe Method tutorial.

**ADD PLOT**

## Usage instructions
1. Open ROOT in terminal
```
root
```

2. Compile `muresolution.cc`, `rochcor2012wasym.cc` and `Analysis.C`
```
.L RochesterCorrections/muresolution.cc++
.L RochesterCorrections/rochcor2012wasym.cc++
.L RochesterCorrections/Test/Analysis.C+
```

3. Run the main function
```
Analysis pf
pf.main()
```
4. Compile `Plot.C` and run the main function
```
.L RochesterCorrections/Test/Plot.C+
main()
```
