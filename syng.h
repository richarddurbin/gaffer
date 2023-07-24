/*  File: syng.h
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2023
 *-------------------------------------------------------------------
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Jul 24 18:08 2023 (rd109)
 * Created: Mon May 29 08:19:18 2023 (rd109)
 *-------------------------------------------------------------------
 */

#include "utils.h"
#include "array.h"
#include "dict.h"
#include "hash.h"

typedef struct {
  int n ;		// number of copies of the syncmer
} Sync ;

typedef struct {
  int sa, sb ;          // syncs, encode as -x-1 if -ve orientation, made canonical as below
  int n ;		// number of copies of the link
  int overlap ;
} Link ;

Array syncs ;
Array links ;
DICT  *syncDict ;
Hash  linkHash ;

static char *syngSchemaText =
  "1 3 def 1 0               schema for syng\n"
  ".\n"
  "P 3 seg                   SEGMENT file\n"
  "S 7 syncmer               file of syncmers\n"
  "D h 3 3 INT 3 INT 3 INT   k, w, seed for the seqhash - here k = smer length, w = syncmer length+1\n"
  "O S 1 3 INT               length\n"
  "D K 1 3 INT               kmer count\n"
  //  "D L 1 8 INT_LIST          list of reads containing the syncmer\n"
  //  "D O 1 6 STRING            orientation of syncmer in read '+'|'-'\n"
  ".\n"
  "P 3 seq                   SEQUENCE\n"
  "S 6 segseq                segment sequence - objects are 1:1 with those in seg file\n"
  "S 7 syncseq               syncmer sequence\n"
  "S 7 readseq               read sequence\n"
  "S 9 contigseq             contig sequence\n"
  "O S 1 3 DNA               sequence of the syncmer\n"
  ".\n"
  "P 6 seqsyn                sequences of syncmers\n"
  "S 7 readsyn               read sequence in syncmers\n"
  "S 9 contigsyn             contig sequence in syncmers\n"
  "D h 3 3 INT 3 INT 3 INT   k, w, seed for the seqhash - here k = smer length, w = syncmer length+1\n"
  "O S 1 8 INT_LIST          sequence: list of syncmer seg ids\n"
  "D P 1 8 INT_LIST          positions of the syncmers\n"
  "D O 1 6 STRING            orientations of the syncmers\n"
  "D R 1 3 INT               index of read in original read file\n"
  ".\n"
  "P 4 link                            LINK (default, or a JUMP if J is present)\n"
  "O L 4 3 INT 4 CHAR 3 INT 4 CHAR     s1 dir1 s2 dir2 - s1,2 are indices in seg file, dir='+'|'-'\n"
  "D O 1 3 INT                         overlap - else presume abut\n"
  "D J 1 3 INT                         a JUMP not a link - int is the gap, 0 if unknown (GFA '*'), -ve for possible overlap\n"
  "D G 1 6 STRING                      cigar string - else presume exact (only meaningful if overlapping)\n"
  "D Q 1 3 INT                         MQ mapping quality\n"
  "D M 1 3 INT                         NM number of mismatches\n"
  "D R 1 3 INT                         RC read count\n"
  "D F 1 3 INT                         FC fragment count\n"
  "D K 1 3 INT                         KC k-mer count\n"
  "D I 1 6 STRING                      ID edge identifier (deprecated)\n"
  "D C 0                               SC if present then a shortcut - see spec - used for scaffolds\n" ;

/****************** end of file ********************/
