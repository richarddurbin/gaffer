# syncmer_index

Code for building an enhanced suffix array of a set of syncmer-converted sequences.

## Building

From this directory: make all

## Using

Paths relative to current directory.

Build syncmers from fasta:

> ../syng -o mysync test.fa

Build index from syncmers: 

> ./syncmer_index build mysync.1readsyn myindex

Generates files

myindex.sa - suffix array

myindex.lcp - LCP array

myindex.da - document array (which sequence the symbol came from)


ASCII dump of index contents:

> ./syncmer_index build mysync.1readsyn myindex


## Acknowledgements

The use of gsa-is https://github.com/felipelouza/gsa-is is gratefully acknowledged.

See gsa-is/LICENSE for license terms (MIT) and gsa-is/README.md for citation details.
