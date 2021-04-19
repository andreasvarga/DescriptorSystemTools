function makemex
%
% Function for generating gateway functions from MEX files for descriptor systems (Linux version, 64 bit).
% This file should be located in the directory containing the source MEX files, which is a subdirectory
% of the directory where the libraries slicot.a and descript.a (if any) are located.
%
flags = 'FFLAGS="$FFLAGS -fPIC -fno-omit-frame-pointer -fdefault-integer-8" -largeArrayDims';
descript_mex_src = '';
libslicot = '../slicot.a';
libdescript = '../descript.a';
descript_mex = {
    'sl_glme', ...
    'sl_gminr', ... 
    'sl_gsep', ...
    'sl_gstra', ...
    'sl_gzero', ...
    'sl_klf', ... 
    'sl_place', ...
    };
%
for k = 1:length(descript_mex)
    file = descript_mex{k};
    fprintf( 'mex %s %s%s.F %s %s -lmwlapack -lmwblas\n', flags, descript_mex_src, file, libdescript, libslicot );
    %eval( sprintf( 'mex %s %s%s.F %s %s -lmwlapack -lmwblas\n', flags, descript_mex_src, file, libdescript, libslicot ) );
end

%!-----

