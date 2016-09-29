# midas2ROOT
Modified version of a code to convert UK MIDAS data into ROOT TTrees

This code can be run from the command line by doing:

root "makeRootTree.C(\"RXXX_0\",\"\",0,0)"

where 'XXX' should be replaced by the run number for the current run.

Alternatively, from ROOT, do:

.L makeRootTree.C
makeRootTree("RXXX_0","",0,0)

If you set the second argument above to be "cal", you can set gains and offsets for electronics channels without too much of a problem.
