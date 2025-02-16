# Ta-Sn-O dataset
Many dataset here are too large to be on github. Please refer to the data repository mentioned in the paper.

## TaSnO_data.npy
Full data of the Ta-Sn-O data. There are 625 stripes, each having 201 XRD spectra. Each spectra has 1024 pixels.

## TaSnO_DwellTpeak.csv
Includes the dwell and peak temperature of each anneal, and the cation ratio of each stripe.

## TaSnO_conds.json
Stores x position, y position, dwell and tpeak of each of the 625 streiips.

## TaSnO_EM=5.json
Phase label result of CrystalShift using expectation maximization loop of 5.

## TaSnO_HW.npy
NMF basis matrices (625 * 4 * 1024) of the data in `TaSnO_data.npy`.

## TaSnO_Hs.npy
NMF activation matrices (625 * 4 * 201) of the data in `TaSnO_data.npy`.

## TaSnO_q.npy
Q valule of the `TaSnO_data.npy`

## sticks.csv
Input file for creating `CrystalPhase` object.

## gt.yaml
Includes the hand labeled ground truth phase labels.
