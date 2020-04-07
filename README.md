# alpha-chenghui-monte-carlo

This is Chenghui's Monte Carlo code from the previous alpha measurement. The main questions it answers are related to the French effect, and also to contrast loss. It takes in real laser intensity images and uses the local intensity during each pulse for each atom by interpolating between pixel values. 

The code is a combination of c and python notebooks. The integration for Bragg is compiled in C to make it very fast, then there are python wrappers (Cython) for using these integration functions in python. The actual monte carlo simulation was run in C to make it easier to run, but can optionally be rewritten entirely in python.

Things the code does not include: transverse momentum kicks, diffraction phases. 
