# MuonCorrections

## Description

This page contains instructions for using the [Rochester corrections](https://twiki.cern.ch/twiki/bin/viewauth/CMS/RochcorMuon) for muon momentum scale and resolution. These instructions are for the Run1 corrections and a [2012 dataset](http://opendata.cern.ch/record/12341) is used as an example.

`RochesterCorrections` contains the official code for the corrections. `Analysis.C` reads the data, applies the corrections and produces an output-file with the corrected data. `Plot.C` creates the plot below that shows that the corrections have been applied correctly.

**ADD PLOT**

## Usage instructions
1. Open ROOT in terminal
```
root -l
```

2. Compile *muresolution.cc*, *rochcor2012wasym.cc* and *Analysis.C*
```
.L muresolution.cc++
.L rochcor2012wasym.cc++
.L Analysis.C+
```

3. Run the main function
```
Analysis pf
pf.main()
```
