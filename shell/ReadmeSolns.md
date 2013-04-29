# Solutions to execises in Introduction to Shell

* * * *
**Short Exercise**

1. Use the manual page for `ls` to guess what you would expect from
using the arguments `-l`, '-t', '-r' at the same time.

``` 
This will produce a long listing `-l`, sorted by the time that the
file last changed `-t`, and in the reverse order `-r`.
```

2. Try the following and see if you can figure out what they do, either by examining the results or consulting the manual page.
   * `ls -lS` (equivalent to `ls -l -S`)
   * `ls -lt` (equivalent to `ls -l -t`)
   * `ls -1`  (that's the number one, not a letter 'ell')

```
In this order, these will produce;
* `ls -lS` lists files in a long listing sorted by their size
* `ls -lt` lists files in a long listing sorted by its time
* `ls -1`  lists files in a brief listing, but 1 per line
```

* * * *
