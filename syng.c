/*  File: syng.c
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2023
 *-------------------------------------------------------------------
 * Description: syncmer-based graph assembler
 * Exported functions:
 * HISTORY:
 * Last edited: May 29 13:45 2023 (rd109)
 * Created: Thu May 18 11:57:13 2023 (rd109)
 *-------------------------------------------------------------------
 */

// AGENDA
// need to add in edges, with counts
// need to consolidate edges to form unitigs
// when do I convert back to non-hoco?
// need some ability to remove nodes or reads, or error-correct

#include "seqio.h"
#include "seqhash.h"
#include "ONElib.h"

#include "syng.h"

static void canonicaliseLink (int ia, int ib, int *sa, int *sb)
{ int ka = (ia >= 0) ? ia : -ia-1 ;
  int kb = (ib >= 0) ? ib : -ib-1 ;
  if (ka < kb) { *sa = ia ; *sb = ib ; }
  else if (kb > ka) { *sa = -ib-1 ; *sb = -ia-1 ; }
  else // they are equal
    if (ia >= 0) { *sa = ia ; *sb = ib ; }
    else { *sa = -ib-1 ; *sb = -ib-1 ; }
}

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

  syncs = arrayCreate (1<<20, Sync) ;
  syncDict = dictCreate (1<<20) ;
  links = arrayCreate (1<<20, Link) ;
  linkHash = hashCreate (1<<24) ;

  storeCommandLine (argc, argv) ;
  argc-- ; ++argv ;
  if (!argc) { printf ("%s",usage) ; exit (0) ; }

  while (argc > 0 && **argv == '-')
    if (!strcmp (*argv, "-w") && argc > 1) { w = atoi(argv[1]) ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (*argv, "-k") && argc > 1) { k = atoi(argv[1]) ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (*argv, "-seed") && argc > 1) { seed = atoi(argv[1]) ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (*argv, "-o") && argc > 1) { outPrefix = argv[1] ; argc -= 2 ; argv += 2 ; }
    else die ("unknown parameter %s\n%s", *argv, usage) ;

  Seqhash *sh = seqhashCreate (k, w+1, seed) ; // need the +1 here - not nice, could fix

  if (argc != 1) die ("need to give a read sequence file\n%s", usage) ;
  if (!(sio = seqIOopenRead (*argv, dna2indexConv, false)))
    die ("failed to open read sequence file %s", *argv) ;
  char *syncString = new0 (w+k+1, char) ;
  if (!(seqsyn = oneFileOpenWriteNew (fnameTag (outPrefix,"1seqsyn"), schema, "seqsyn", true, 1)))
      die ("failed to open %s.1seqsyn to write", outPrefix) ;
  oneAddProvenance (seqsyn, "syng", "0.0", getCommandLine()) ;
  oneAddReference (seqsyn, *argv, 0) ; // I don't know what to put for count - maybe 0 is correct?
  Array readSync = arrayCreate(64,I64) ;
  Array readPos = arrayCreate(64,I64) ;
  Array readDir = arrayCreate(64,char) ;
  U64 nSeq = 0, totSeq = 0, totSync = 0 ;
  while (seqIOread (sio))
    { char* seq = sqioSeq(sio) ;
      SeqhashIterator *sit = syncmerIterator (sh, seq, sio->seqLen) ;
      int pos, iSync, iLastSync, iLink ;
      arrayMax(readSync) = 0 ; arrayMax(readPos) = 0 ; arrayMax(readDir) = 0 ;
      while (syncmerNext (sit, 0, &pos, 0))
	{ array(readPos,arrayMax(readPos),I64) = pos ;
	  int i = pos-1, j = pos+w+k ; while (seq[++i] == 3 - seq[--j]) ;
	  if (seq[i] < 3 - seq[j])
	    { array(readDir,arrayMax(readDir),char) = '+' ;
	      //	      printf ("+ pos %d i %d %c %c\n", pos, i, acgt[seq[i]], acgt[seq[j]]) ;
	      for (i = 0 ; i < w+k ; ++i) syncString[i] = acgt[seq[pos+i]] ;
	      dictAdd (syncDict, syncString, &iSync) ;
	      arrayp(syncs,iSync,Sync)->n++ ;
	      array(readSync,arrayMax(readSync),I64) = iSync ;
	    }
	  else
	    { array(readDir,arrayMax(readDir),char) = '-' ;
	      // printf ("- pos %d i %d %c %c\n", pos, i, acgt[seq[i]], acgt[seq[j]]) ;
	      for (i = 0 ; i < w+k ; ++i) syncString[w+k-i-1] = tgca[seq[pos+i]] ;
	      dictAdd (syncDict, syncString, &iSync) ;
	      arrayp(syncs,iSync,Sync)->n++ ;
	      array(readSync,arrayMax(readSync),I64) = iSync ;
	      iSync = -iSync-1 ;
	    }
	  if (arrayMax(readPos) > 1) // add the link
	    { int sa, sb ;
	      canonicaliseLink (iLastSync, iSync, &sa, &sb) ;
	      Link *link ;
	      if (hashAdd (linkHash, HASH_INT2(sa,sb), &iLink))
		{ link = arrayp(links,iLink,Link) ;
		  link->sa = sa ; link->sb = sb ;
		}
	      else link = arrp(links,iLink,Link) ;
	      ++link->n ; 
	    }
	  iLastSync = iSync ;
	}
      oneWriteLine (seqsyn, 'S', arrayMax(readSync), arrp(readSync,0,I64)) ;
      oneWriteLine (seqsyn, 'P', arrayMax(readPos), arrp(readPos,0,I64)) ;
      oneWriteLine (seqsyn, 'O', arrayMax(readDir), arrp(readDir,0,char)) ;
      oneInt(seqsyn, 0) = nSeq ; oneWriteLine (seqsyn, 'R', 0, 0) ;
      nSeq++ ; totSeq += sio->seqLen ; totSync += arrayMax(readSync) ;
    }

  printf ("read %" PRIu64 " sequences, total length %" PRIu64 ""
          ", yielding %" PRIu64 " instances of %d syncmers, average %.2f coverage\n", 
          nSeq, totSeq, totSync, arrayMax(syncs), totSync / (double)arrayMax(syncs)) ; 
  
  oneFileClose (seqsyn) ;
  arrayDestroy (readSync) ;
  arrayDestroy (readPos) ;
  arrayDestroy (readDir) ;

  if (!(seg = oneFileOpenWriteNew (fnameTag (outPrefix,"1seg"), schema, "seg", true, 1)))
    die ("failed to open %s.1seg to write", outPrefix) ;
  oneAddProvenance (seg, "syng", "0.0", getCommandLine()) ;
  int i ; 
  for (i = 0 ; i < arrayMax(syncs) ; ++i)
    { Sync *s = arrp(syncs, i, Sync) ;
      oneInt(seg,0) = w+k ; oneWriteLine (seg, 'S', 0, 0) ;
      oneInt(seg,0) = s->n ; oneWriteLine (seg, 'K', 0, 0) ;
    }
  oneFileClose (seg) ;

  if (!(segseq = oneFileOpenWriteNew (fnameTag (outPrefix,"1segseq"), schema, "segseq", true, 1)))
    die ("failed to open %s.1segseq to write", outPrefix) ;
  oneAddProvenance (segseq, "syng", "0.0", getCommandLine()) ;
  for (i = 0 ; i < dictMax(syncDict) ; ++i)
    oneWriteLine (segseq, 'S', w+k, dictName (syncDict, i)) ;
  oneFileClose (segseq) ;

  arrayDestroy (syncs) ;
  dictDestroy (syncDict) ;
  arrayDestroy (links) ;
  hashDestroy (linkHash) ;
	  
  fprintf (stderr, "total: ") ; timeTotal (stderr) ;
  return 0 ;
}

/*********************** end of file **********************/
