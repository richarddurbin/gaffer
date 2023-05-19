/*  File: syng.c
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2023
 *-------------------------------------------------------------------
 * Description: syncmer-based graph assembler
 * Exported functions:
 * HISTORY:
 * Last edited: May 19 08:01 2023 (rd109)
 * Created: Thu May 18 11:57:13 2023 (rd109)
 *-------------------------------------------------------------------
 */

#include "seqio.h"
#include "seqhash.h"
#include "ONElib.h"

#include "utils.h"
#include "array.h"
#include "dict.h"

typedef struct {
  int n ;		// number of copies of the syncmer
  
} Syncmer ;

Array syncmer ;
DICT  *syncDict ;

static char *syngSchemaText ; // at end of file

static char usage[] =
  "Usage: syng <options> <read sequence file>\n"
  "  -w <syncmer length>    : [1023]\n"
  "  -k <smer length>       : [16] must be under 32\n"
  "  -seed <seed>           : [7] for the hashing function\n"
  "  -o <outfile prefix>    : [syng]\n" ;

int main (int argc, char *argv[])
{
  int w = 1023 ;
  int k = 16 ;
  int seed = 7 ;
  char *outPrefix = "syng-out" ;
  SeqIO *sio ;
  OneSchema *schema = oneSchemaCreateFromText (syngSchemaText) ;
  OneFile *seg, *segseq, *seqsyn, *link ;
  char *acgt = "acgt", *tgca = "tgca" ; // normal and complemented

  timeUpdate (0) ;

  syncDict = dictCreate (1<<20) ;
  syncmer = arrayCreate (1<<20, Syncmer) ;

  argc-- ; ++argv ;
  if (!argc) { printf ("%s",usage) ; exit (0) ; }

  while (argc > 0 && **argv == '-')
    if (!strcmp (*argv, "-w") && argc > 1) { w = atoi(argv[1]) ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (*argv, "-k") && argc > 1) { k = atoi(argv[1]) ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (*argv, "-seed") && argc > 1) { seed = atoi(argv[1]) ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (*argv, "-o") && argc > 1) { outPrefix = argv[1] ; argc -= 2 ; argv += 2 ; }
    else die ("unknown parameter %s\n%s", *argv, usage) ;

  Seqhash *sh = seqhashCreate (k, w, seed) ;

  if (argc != 1) die ("need to give a read sequence file\n%s", usage) ;
  if (!(sio = seqIOopenRead (*argv, dna2indexConv, false)))
    die ("failed to open read sequence file %s", *argv) ;
  char *syncString = new0 (w+k+1, char) ;
  if (!(seqsyn = oneFileOpenWriteNew (fnameTag (outPrefix,"1seqsyn"), schema, "seqsyn", true, 1)))
      die ("failed to open %s.1seqsyn to write", outPrefix) ;
  Array readSync = arrayCreate(64,I64) ;
  Array readPos = arrayCreate(64,I64) ;
  Array readDir = arrayCreate(64,char) ;
  U64 nSeq = 0, totSeq = 0, totSync = 0 ;
  while (seqIOread (sio))
    { char* seq = sqioSeq(sio) ;
      SeqhashIterator *sit = syncmerIterator (sh, seq, sio->seqLen) ;
      int pos, index ;
      // printf ("%s %" PRId64 " :", sqioId(sio), sio->seqLen) ;
      arrayMax(readSync) = 0 ; arrayMax(readPos) = 0 ; arrayMax(readDir) = 0 ;
      while (syncmerNext (sit, 0, &pos, 0))
	{ array(readPos,arrayMax(readPos),I64) = pos ;
	  int i = pos-1, j = pos+w+k ; while (seq[++i] == 3 - seq[--j]) ;
	  if (seq[i] < seq[j])
	    { array(readDir,arrayMax(readDir),char) = '+' ;
	      for (i = 0 ; i < w+k ; ++i) syncString[i] = acgt[seq[pos+i]] ;
	      dictAdd (syncDict, syncString, &index) ;
	      arrayp(syncmer,index,Syncmer)->n++ ;
	      array(readSync,arrayMax(readSync),I64) = index ;
	      // printf (" %d+%d", index, pos) ;
	    }
	  else
	    { array(readDir,arrayMax(readDir),char) = '-' ;
	      for (i = 0 ; i < w+k ; ++i) syncString[w+k-i-1] = tgca[seq[pos+i]] ;
	      dictAdd (syncDict, syncString, &index) ;
	      arrayp(syncmer,index,Syncmer)->n++ ;
	      array(readSync,arrayMax(readSync),I64) = index ;
	      //      printf (" %d-%d", index, pos) ;
	    }
	}
      //      printf ("\n") ;
      oneWriteLine (seqsyn, 'S', arrayMax(readSync), arrp(readSync,0,I64)) ;
      oneWriteLine (seqsyn, 'P', arrayMax(readPos), arrp(readPos,0,I64)) ;
      oneWriteLine (seqsyn, 'O', arrayMax(readDir), arrp(readDir,0,char)) ;
      nSeq++ ; totSeq += sio->seqLen ; totSync += arrayMax(readSync) ;
    }

  printf ("read %" PRIu64 " sequences, total length %" PRIu64 ""
          ", yielding %" PRIu64 " instances of %d syncmers, average %.2f coverage\n", 
          nSeq, totSeq, totSync, arrayMax(syncmer), totSync / (double)arrayMax(syncmer)) ; 
  
  oneFileClose (seqsyn) ;
  arrayDestroy (readSync) ;
  arrayDestroy (readPos) ;
  arrayDestroy (readDir) ;

  if (!(seg = oneFileOpenWriteNew (fnameTag (outPrefix,"1seg"), schema, "seg", true, 1)))
    die ("failed to open %s.1seg to write", outPrefix) ;
  int i ; 
  for (i = 0 ; i < arrayMax(syncmer) ; ++i)
    { Syncmer *s = arrp(syncmer, i, Syncmer) ;
      oneInt(seg,0) = w+k ; oneWriteLine (seg, 'S', 0, 0) ;
      oneInt(seg,0) = s->n ; oneWriteLine (seg, 'R', 0, 0) ;
    }
  oneFileClose (seg) ;

  if (!(segseq = oneFileOpenWriteNew (fnameTag (outPrefix,"1segseq"), schema, "segseq", true, 1)))
    die ("failed to open %s.1segseq to write", outPrefix) ;
  for (i = 0 ; i < dictMax(syncDict) ; ++i)
    oneWriteLine (segseq, 'S', w+k, dictName (syncDict, i)) ;
  oneFileClose (segseq) ;

  fprintf (stderr, "total: ") ; timeTotal (stderr) ;
  return 0 ;
}

static char *syngSchemaText =
  "1 3 def 1 0               schema for syng\n"
  ".\n"
  "P 3 seg                   SEGMENT file\n"
  "S 7 syncmer               file of syncmers\n"
  "O S 1 3 INT               length\n"
  "D R 1 3 INT               read count\n"
  "D L 1 8 INT_LIST          list of reads containing the syncmer\n"
  "D P 1 8 INT_LIST          position of syncmer in each read\n"
  "D O 1 6 STRING            orientation of syncmer in read '+'|'-'\n"
  ".\n"
  "P 3 seq                   SEQUENCE\n"
  "S 6 segseq                segment sequences - objects are 1:1 with those in seg file\n"
  "O S 1 3 DNA               sequence of the syncmer\n"
  ".\n"
  "P 6 seqsyn                sequences of syncmers\n"
  "O S 1 8 INT_LIST          sequence: list of syncmer seg ids\n"
  "D P 1 8 INT_LIST          positions of the syncmers\n"
  "D O 1 6 STRING            orientations of the syncmers\n"
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

/*********************** end of file **********************/
