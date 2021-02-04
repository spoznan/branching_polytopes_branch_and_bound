# Scripts and Data for branch-and-bound on branching polytopes

## Description

This repository contains the code and data used for the article 'Improving RNA branching predictions: advances and limitations' by Svetlana Poznanovic ́, Carson Wood, Michael Cloer, and Christine E. Heitsch (link will be made available). 


## Dependencies

These scripts are written for SageMath in accompaniement with Python 3. They are known to 
work with SageMath 9.2 and Python 3.7, although they likely are functional with earlier versions too. 

## Usage

To run, use the bb_main.sage script. For example:

```sage ~/bb_main.sage```

This will find the optimal parameters for the sequences listed in the file specified as `seq_list` in the main code.
The optimal parameters are the parameters achieving the best average accuracy value across the specified data set.
The region containing those parameters will be specified in `~/MergeData/FinalOutput.txt`.
Additionally, data about algorithm computations and runtimes will be listed in `~/MergeData/SummativeData.txt`.

## Data 

We are providing input data used for testing for 50 tRNA and 50 5S sequences that can be used to run the code in the `~/Data` folder. These include the branching polytopes in the `.rnapoly` files computed using the [pmfe](https://github.com/gtDMMB/pmfe) and the optimal structures computed for all vertices of the branching polytopes.

