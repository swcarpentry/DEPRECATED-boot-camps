Compile each program without optimization first.

For simpleTest.cc, run this line to see errors in this code. 

::

  valgrind --track-origins=yes --leak-check=full ./simpleTest 300 300


We also have a cache test line. Run this line to see the cache errors.

::

  valgrind --tool=cachegrind ./a.out 0 1000 100000

There are two paths in this code. If the first input is 1, it runs a cache-sensitive version of the loop. 
If it is 0, it runs a cache-insensitive version.

FYI: on the Trieste lab machines, this is what cache looks like:

::

  guy ~>dmesg | grep cache
  CPU: L1 I cache: 32K, L1 D cache: 32K
  CPU: L2 cache: 6144K
  CPU: L1 I cache: 32K, L1 D cache: 32K
  CPU: L2 cache: 6144K

You can run the same command to see cache on your linux machine. Another way to see the exact cache setup that 
valgrind found is the following:

::

  cg_annotate --auto=yes cachegrind.out.21960

Note that your cachegrind.out will have a different number. This command is also handy because it shows which functions caused cache
misses.



