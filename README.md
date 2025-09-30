# Welcome to coran, my COrrelation ANalysis repository!
This will basically *always* be a work in progress...


## Design Decisions
coran is effectively a library, as it provides a set of common utilities for performing angular correlation analyses (especially ones that require neutral particle/V0 reconstruction)

These utilities come in the form of classes, where each class has a "job" that is performed during the initialization of the class. The initialization input contains the bare minimum amount of data required to accomplish the job, and the job's "output" is accessible through @property getters for each class. 

As a concrete example, let's consider the DeltaPhiFit class under

I'm not a huge fan of OOP in general, and this design choice is *very* much object-oriented.
