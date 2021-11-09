# Crystallography_based_shifting
## TODOs
- [x] Data Structure (7/30)
- [ ] Test cases for optimization
- [x] Using Optimization API
- [x] Figure out how to construct residual function with an array of phases
- [x] Implement function for q calculation for different crystals

## Development log
### 8/1
* Capable of creating input file with python
* Using crystal and peak info to reconstruct spectrum

### 8/3
* Updated optimize.jl
* Tested peak position
* Start implementing how to calculate residual with spectrum and params

### 8/25
* Test on real data
* Each spectrum takes ~4 secs
* Speedup needed
* Want result object

### 11/9
* Bug fixes
* Still need to find proper values for regularization