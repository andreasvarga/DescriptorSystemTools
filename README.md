# **DSTOOLS - Descriptor System Tools for MATLAB**  

## About 

`DSTOOLS` is a collection of MATLAB functions for the operation on and manipulation of rational transfer-function matrices via their
descriptor system realizations. The `DSTOOLS` collection relies on the
Control System Toolbox and several mex-functions based on the Systems and Control Library
[SLICOT](http://slicot.org/). The underlying mex-functions have been implemented in collaboration with Vasile Sima.

Many of the implemented functions are based on the computational procedures described in Chapter 10 of the book:

* Andreas Varga, "[Solving Fault Diagnosis Problems - Linear Synthesis Techniques](http://www.springer.com/us/book/9783319515588)", vol. 84 of Studies in Systems, Decision and Control, Springer International Publishing, xxviii+394, 2017.

The User's Guide of the current version of the `DSTOOLS` collection is available in the DSTOOLS repository, in the file `dstoolsdoc.pdf` (the documentation of the version V0.71 is available on [arXiv](https://arxiv.org/abs/1707.07140)).  Additionally, the M-files of the functions are self-documenting and a detailed documentation of each function can be obtained online by typing help with the corresponding M-file name.
The above book provides additional information on the mathematical background on rational matrices and descriptor systems,
and gives detailed descriptions of most of the underlying procedures.

The current release of `DSTOOLS` is version 0.75, dated January 31, 2021.

## Requirements

The codes have been developed under MATLAB 2015b and have been tested with MATLAB 2016a through 2020b. To use the functions, the Control System Toolbox must be installed in MATLAB running under 64-bit Windows 7, 8, 8.1 or 10.

## License

* See `[LICENSE](https://github.com/andreasvarga/DescriptorSystemTools/blob/main/LICENSE)` for licensing information.

* Please cite `DSTOOLS` as "A. Varga. DSTOOLS - The Descriptor System Tools for MATLAB.
[https://sites.google.com/view/andreasvarga/home/software/dstools](https://sites.google.com/view/andreasvarga/home/software/dstools), 2019."

* Please cite the documentation of `DSTOOLS` as "A. Varga. Descriptor System Tools (DSTOOLS) User's Guide, ArXiv eprint [arXiv:1707.07140](https://arxiv.org/abs/1707.07140), September 2018."
