# gaffer - a GFA toolkit

The main product here is **gaffer**, with currently a single source file gaffer.c.
Its initial task is to read in a HiFiAsm unitig file, carry out various QC and cleanup operations, read and write assemblies into ONEcode files, and convert from an overlap to a blunt-ended graph.

The project also contains relevant versions of various Durbin package utility programs:

- **ONEview** to view ONEcode files, and convert them between binary and ascii.  By convention .1seq files are binary, .seq files are ascii.
- **ONEstat** to give stats for 1code files.
- **composition** to give information about the data in DNA sequence files: fasta[.gz], fastq[.gz], 1seq (if compiled with -DONEIO, as by default in this project), SAM/BAM/CRAM (if compiled with -DBAMIO, see below), and a custom binary sequence file format (now deprecated with 1seq preferred).
- **seqconvert** to convert between DNA file types.  This also will homopolymer compress (-H "hoco"), and if writing to a 1seq file while doing this will store the offsets in the original file of each new sequence base position, allowing to unhoco (-U) back to the original sequence.  Because ONEcode files compress all lists, this 

For each program, running it without any arguments gives usage information.

## Building
```
  git clone https://github.com/richarddurbin/gaffer.git
  make
```

If you want to be able to read SAM/BAM/CRAM files then you need to install htslib in a parallel directory and use Makefile.bam:
```
  cd ..
  git clone https://github.com/samtools/htslib.git
  cd htslib
  autoreconf -i  # Build the configure script and install files it uses
  ./configure    # Optional but recommended, for choosing extra functionality
  make
  make install
  cd ../gaffer
  make clean
  make -f Makefile.bam
```

## Libraries

There are various useful libraries, with header files:

- **ONElib.[hc]** supports ONEcode file reading and writing, with implementation in the single file ONElib.c with no dependendencies.  See https://github.com/VGP/vgp-tools/tree/master/Core for further information.
- **seqio.[hc]** supports reading, writing DNA files with a few other basic operations.  Implementation in seqio.c, with dependencies on utils.[hc], libz, ONElib and htslib depending on compile operations.
- **utils.[hc]** some very low level type definitions (e.g. I8 to I64 and U8 to U64), die(), warn(), new(), new0(), and a timing package.  No dependencies beyond normal C run time library.  NB there is a handy fzopen() which will silently open .gz files as standard files, but this depends on funopen() which is not available on all systems.  If this does not compile/link then you will need to undefine WITH_ZLIB to link.
- **array.[ch], dict.[ch], hash.[ch]** respectively implement advanced language style extendable arrays, dictionaries (hashes of strings) and general hashes of basic types (up to 64-bit).

## Synopsis

An example usage pattern is given below.  Note that as of 17/8/2022 the bluntification step has bugs, so the result is not correct.

```
> gaffer
Usage: gaffer <commands>
   -read  <stem>   : read ONE files
   -write <stem>   : write ONE files
   -schema         : print ONE schema
   -readGfa <gfa file> : only reads S and L lines for now
   -dna <sequence file matching S line names>
   -removeBadLinks : (for now) remove imperfect overlaps
   -blunt          : makes new non-overlapping graph
   -depth <gfa file with A lines>
   -extend         : adds match blocks for shared incoming edges - needs seq
       
> gaffer -readGfa TEST/daStaPalu1.r_utg.gfa.gz -write daStaPalu1-r
WARNING: ** self link from/to utg010600c 21198 overlap 0 - ignoring
WARNING: ** self link from/to utg010600c 21199 overlap 0 - ignoring
WARNING: ** self link from/to utg017205c 34408 overlap 0 - ignoring
WARNING: ** self link from/to utg017205c 34409 overlap 0 - ignoring
WARNING: ** self link from/to utg018266c 36530 overlap 0 - ignoring
WARNING: ** self link from/to utg018266c 36531 overlap 0 - ignoring
read 48356 S lines and 131656 L lines from GFA file TEST/daStaPalu1.r_utg.gfa.gz
removed 62412 duplicate links, leaving 69244
readGfa: user   122.342262      system  1.100602        max_RSS 1079476224      memory  5966385532
wrote 24178 objects to daStaPalu1-r.1seg    // 1code file of the segments
wrote 24178 objects to daStaPalu1-r.1sgs    // 1code seq file of the DNA sequences of the segments
wrote 69244 objects to daStaPalu1-r.1lnk    // 1code file of the links
write: user     2.277581        system  1.384311        max_RSS 2242887680      memory  5966385549
total: user     124.629004      system  3.074676        max_RSS 3326738432      memory  5966385549

> gaffer -read daStaPalu1-r -removeBadLinks -blunt -write daStaPalu1-blunt
read file daStaPalu1-r.1lnk with 24178 segs
read file daStaPalu1-r.1lnk with 69244 links
read file daStaPalu1-r.1sgs with 24178 sequences
read: user      9.287462        system  0.607038        max_RSS 1794195456      memory  8941142529
10414 imperfect overlaps removed, 58830 remain
removeBadLinks: user    1.296138        system  0.227495        max_RSS 0       memory  8941142529
243080 initial new links before compression
made blunt gfa with 231179 seqs, 160533 links and 48356 walks
all walks check out
blunt: user     34.994954       system  1.648977        max_RSS 2273067008      memory  9007107609
wrote 115590 objects to daStaPalu1-blunt.1seg
wrote 160533 objects to daStaPalu1-blunt.1lnk
wrote 48356 objects to daStaPalu1-blunt.1wlk
write: user     0.113076        system  0.462926        max_RSS 0       memory  9007108270
total: user     45.695350       system  2.947524        max_RSS 4067262464      memory  9007108270

> seqconvert -o daStaPalu1-r.fa.gz daStaPalu1-r.1sgs    // -fa -z are implied by the output filename
reading from file type onecode  with 24178 sequences totLen 2979016643
written 24178 sequences to file type fasta, total length 2979016643, max length 6381111
user    484.026538      system  2.124296        max_RSS 54951936        memory  46323317    // writing gzip files is very slow

> seqconvert -1 -H -o daStaPalu1-r.hoco.1seq daStaPalu1-r.1sgs 
reading from file type onecode  with 24178 sequences totLen 2979016643
written 24178 sequences to file type onecode, total length 2067571198, max length 4464555
user    47.692720       system  2.497245        max_RSS 2053259264      memory  24086839758

> seqconvert -fa -z -U -o daStaPalu1-r.unhoco.fa.gz daStaPalu1-r.hoco.1seq
reading from file type onecode  with 24178 sequences totLen 2067571198
written 24178 sequences to file type fasta, total length 2979016643, max length 6381111
user    559.758163      system  2.878647        max_RSS 107282432       memory  108422711

> ls -l TEST/daStaPalu1.r_utg.gfa.gz daStaPalu* 
-rw-r--r--  1 rd  staff  847123532 Mar 23 13:34 TEST/daStaPalu1.r_utg.gfa.gz
-rw-r--r--  1 rd  staff    1707152 Aug 17 12:35 daStaPalu1-blunt.1lnk
-rw-r--r--  1 rd  staff     460153 Aug 17 12:35 daStaPalu1-blunt.1seg
-rw-r--r--  1 rd  staff    1188974 Aug 17 12:35 daStaPalu1-blunt.1wlk
-rw-r--r--  1 rd  staff     903589 Aug 17 12:28 daStaPalu1-r.1lnk
-rw-r--r--  1 rd  staff     365580 Aug 17 12:28 daStaPalu1-r.1seg
-rw-r--r--  1 rd  staff  745015260 Aug 17 12:28 daStaPalu1-r.1sgs
-rw-r--r--  1 rd  staff  815204880 Aug 17 13:17 daStaPalu1-r.fa.gz
-rw-r--r--  1 rd  staff  890157073 Aug 17 14:22 daStaPalu1-r.hoco.1seq
-rw-r--r--  1 rd  staff  815204880 Aug 17 14:57 daStaPalu1-r.unhoco.fa.gz

> ONEview -h daStaPalu1-r.1lnk | head  // -h does not show the header, just the body
L 0 > 1767 <
O 12648
L 0 < 18185 <
O 10209
L 0 < 20414 <
O 4869
L 1 > 841 >
O 9344
L 1 > 2020 >
O 7852

> ONEview -H daStaPalu1-r.1sgs   // -H shows only the header, with the schema and object statistics
1 3 seq 1 1
2 3 sgs
.
~ O S 1 3 DNA             sequence: the DNA string
~ D I 1 6 STRING          id - sequence identifier; unnecessary for segments
~ D P 1 8 INT_LIST        kmer ploidy estimates - 0 for error, 1,2... for ploidy, R for repeat
~ D Q 1 6 STRING          base quality - Q values (ascii string = q+33)
.
# S 24178         // count of S lines
@ S 6381111       // maximum length of list (in this case DNA) in an S line
+ S 2979016643    // total lengths of all lists in S lines
.

> composition -t daStaPalu1-r.1sgs
onecode file, 24178 sequences >= 0, 2979016643 total, 123211.87 average, 1884 min, 6381111 max
user    2.566312        system  0.088715        max_RSS 38862848        memory  29542744

> composition -t daStaPalu1-r.fa.gz
fasta file, 24178 sequences >= 0, 2979016643 total, 123211.87 average, 1884 min, 6381111 max
user    19.503632       system  0.221412        max_RSS 16875520        memory  16780520      // note that this is much slower
```
Note that not only are the 1seq files smaller than compressed fasta, they also read much faster because parsing is trivial (composition reads in each sequence, even though the basic counts are available in the 1seq header.  composition can all give other information such as base composition (-b), approximate length distribution (-l), quality score distributions (-q) if they are present.
