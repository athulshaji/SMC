Source Separation with Non Negative Matrix Factorization

In this project, a simple source separation system will be built in MATLAB based on non-negative matrix factorization (NMF).

First, build a system that can perform NMF for a given file. 

Then, use a number of clean, individual source signals (using, for example, the single-source signals in the EBU SQAM database) to train a dictionary for each source. 

Next, we will use those dictionaries for separating mixtures of sources by combining the trained dictionaries into a big dictionary (i.e., the union). 

Mix some sources artificially and use the principles of the NMF for finding the activations (i.e., only the step finding the activations). 

Use these activations to estimate the individual sources and evaluate the performance of the source separation scheme. 

Things to investigate include how the performance depends on the signal-to-interference ratio, how performance depends on knowledge of exactly which instruments are to be separated, and which instruments are easy to separate from one-another and which are hard.
