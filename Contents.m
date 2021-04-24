% DSTOOLS - Descriptor System Tools.
% Version 1.0.0.2            24-April-2021
% Copyright (c) 2016-2021 by A. Varga
%
% Demonstration.
%   DSToolsdemo - Demonstration of DSTOOLS.
%
% System analysis.
%   gpole     - Poles of a LTI descriptor system.
%   gzero     - Zeros of a LTI descriptor system. 
%   gnrank    - Normal rank of the transfer function matrix of a LTI system.
%   ghanorm   - Hankel norm of a proper and stable LTI descriptor system.
%   gnugap    - Nu-gap distance between two LTI systems.
%
% System order reduction.
%   gir       - Reduced order realizations of LTI descriptor systems.
%   gminreal  - Minimal realization of a LTI descriptor system.
%   gbalmr    - Balancing-based model reduction of a LTI descriptor system.
%   gss2ss    - Conversions to SVD-like forms without non-dynamic modes.
%
% Operations on transfer function matrices.
%   grnull    - Right nullspace basis of a transfer function matrix.
%   glnull    - Left nullspace basis of a transfer function matrix.
%   grange    - Range space basis of a transfer function matrix. 
%   gcrange   - Coimage space basis of a transfer function matrix. 
%   grsol     - Solution of the linear rational matrix equation G*X = F.
%   glsol     - Solution of the linear rational matrix equation X*G = F.
%   ginv      - Generalized inverses. 
%   gsdec     - Generalized additive spectral decompositions.
%   grmcover1 - Right minimum dynamic cover of Type 1 based order reduction.
%   glmcover1 - Left minimum dynamic cover of Type 1 based order reduction.
%   grmcover2 - Right minimum dynamic cover of Type 2 based order reduction.
%   glmcover2 - Left minimum dynamic cover of Type 2 based order reduction.
%   gbilin    - Generalized bilinear transformation.
%   gbilin1   - Transfer functions of commonly used bilinear transformations. 
%
% Factorizations of transfer function matrices.
%   grcf      - Right coprime factorization with proper and stable factors.
%   glcf      - Left coprime factorization with proper and stable factors.
%   grcfid    - Right coprime factorization with inner denominator.
%   glcfid    - Left coprime factorization with inner denominator.
%   gnrcf     - Normalized right coprime factorization.
%   gnlcf     - Normalized left coprime factorization.
%   giofac    - Inner-outer/QR-like factorization.
%   goifac    - Co-outer-co-inner/RQ-like factorization.
%   grsfg     - Right spectral factorization of gamma^2*I-G'*G.
%   glsfg     - Left spectral factorization of gamma^2*I-G*G'.
%
% Model-matching problem.
%   grasol    - Approximate solution of the linear rational matrix equation G*X = F.
%   glasol    - Approximate solution of the linear rational matrix equation X*G = F.
%   glinfldp  - Solution of the least distance problem min||G1-X G2||_inf.
%   gnehari   - Generalized Nehari approximation.
%
% Feedback stabilization. 
%   gsfstab    - Generalized state-feedback stabilization.
%   eigselect1 - Selection of a real eigenvalue to be assigned.
%   eigselect2 - Selection of a pair of eigenvalues to be assigned.
%   galoc2     - Generalized pole allocation for second order systems.
%
% Pencil similarity transformations
%   gklf      - Kronecker-like staircase forms of a linear matrix pencil. 
%   gsklf     - Special Kronecker-like form of a system matrix pencil. 
%   gsorsf    - Specially ordered generalized real Schur form.
%
% SLICOT-based mex-functions.
%   sl_gstra  - Descriptor system coordinate transformations.
%   sl_gminr  - Minimal realization of descriptor systems.
%   sl_gsep   - Descriptor system additive spectral decompositions.
%   sl_gzero  - Computation of system zeros and Kronecker structure.
%   sl_klf    - Pencil reduction to Kronecker-like forms.
%   sl_glme   - Solution of generalized linear matrix equations. 
