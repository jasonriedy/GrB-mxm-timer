GrB-mxm-timer
=============

A very simple timer for GrB_mxm that emulates a "k-hop" sequence of
operations as in the [TigerGraph Graph Database
Benchmark](https://github.com/tigergraph/graph-database-benchmark).

Usage
=====

    Usage: GrB-mxm-timer [OPTION]...

    Times the iteration of B = A*B

      -h, --help               Print help and exit
      -V, --version            Print version and exit
      -s, --scale=INT          Scale (log2 # vertices in A)  (default=`20')
      -e, --edgefactor=INT     Edge factor, so # edges = ef * 2^scale
                                 (default=`8')
      -A, --A=FLOAT            R-MAT upper left quadrant probability
                                 (default=`0.55')
      -B, --B=FLOAT            R-MAT upper right & lower left quadrant probability
                                 (default=`0.1')
      -N, --noisefact=FLOAT    Noise factor on each recursion  (default=`0.1')

      -c, --b-ncols=INT        Number of columns in B  (default=`16')
      -C, --b-used-ncols=INT   Number of columns actually used in the initial B
                                 (default=`1')
      -E, --b-nents-col=INT    Number of entries per column in the initial B
                                 (default=`1')

      -k, --khops=STRING       Number of iterations / hops (can be a space-delim
                                 list)  (default=`4')

          --NE-chunk-size=INT  Number of edges to generate in a chunk.
                                 (default=`3000')
          --verbose[=INT]      Provide status updates via stdout.  (default=`1')


Building
========

This should build out of the box on a generic host system.

To link with LucataGraphBLAS, override LDLIBS and possibly LDFLAGS to
point to your build.

To build completely as a MWX, add similar overrides (as well as the
compiler), and call `make TARGET_MWX=1`.

