# Scripts and Data accompanying the article 'Improving RNA branching predictions: advances and limitations' by Svetlana Poznanovic ́, Carson Wood, Michael Cloer, and Christine E. Heitsch

## Description

This repository contains the scripts created for the RNA Brancing parameters project as well as the resulting data. 
`bb_main.sage` is the main script. The input for various sequences are in the `~/Data` folder.
Output will be placed in the 
For more details, please see the associated paper (link will be made available).

## Dependencies

These scripts are written for SageMath in accompaniement with Python 3. They are known to 
work with Sagemath 9.2 and Python 3.7, although they likely are functional with earlier versions. 

## Usage

To run, use the bb_main.sage script. For example:

```sage ~/bb_main.sage```

This will find the optimal parameters for the sequences listed in the file specified as `seq_list` in the main code.
The optimal parameters are the parameters achieving the best average accuracy value across the specified data set.
The region containing those parameters will be specified in `~/MergeData/FinalOutput.txt`.
Additionally, data about algorithm computations and runtimes will be listed in `~/MergeData/SummativeData.txt`.

## Data 

We are providing input data used for testing for 50 tRNA and 50 5S sequences that can be used to run the code.

