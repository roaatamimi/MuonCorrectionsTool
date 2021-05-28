# MuonCorrections

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
