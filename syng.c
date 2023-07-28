/*  File: syng.c
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2023
 *-------------------------------------------------------------------
 * Description: syncmer-based graph assembler
 * Exported functions:
 * HISTORY:
 * Last edited: Jul 27 16:02 2023 (rd109)
 * Created: Thu May 18 11:57:13 2023 (rd109)
 *-------------------------------------------------------------------
 */

// AGENDA
// need to add in edges, with counts
// need to consolidate edges to form unitigs
// when do I convert back to non-hoco?
// need some ability to remove nodes or reads, or error-correct

#include <pthread.h>

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
    else { *sa = -ib-1 ; *sb = -ia-1 ; }
}

typedef struct {
  Seqhash *sh ;
  Array seq ;		// of char, input: concatenated sequences in index 0..3
  Array seqLen ;	// of I64, input: per sequence
  Array pos ;		// of I64, output: concatenated start positions of syncmers
  Array posLen ;	// of I64, output: number of syncmers per sequence
} ThreadInfo ;

static void *threadProcess (void* arg) // find the start positions of all the syncmers
{
  ThreadInfo *ti = (ThreadInfo*) arg ;
  int i ;
  I64 seqStart = 0 ;

  arrayMax(ti->pos) = 0 ;
  arrayMax(ti->posLen) = 0 ;
  
  for (i = 0 ; i < arrayMax(ti->seqLen) ; ++i)
    { I64 seqLen = arr(ti->seqLen, i, I64) ;
      char *seq = arrp(ti->seq, seqStart, char) ;
      seqStart += seqLen ;
      int pos, posStart = arrayMax(ti->pos) ;
      SeqhashIterator *sit = syncmerIterator (ti->sh, seq, seqLen) ;
      while (syncmerNext (sit, 0, &pos, 0)) array(ti->pos, arrayMax(ti->pos), I64) = pos ;
      array(ti->posLen, i, I64) = arrayMax(ti->pos) - posStart ;
    }
  return 0 ;
}

static char usage[] =
  "Usage: syng <options> <read sequence file>\n"
  "  -w <syncmer length>    : [1023]\n"
  "  -k <smer length>       : [16] must be under 32\n"
  "  -seed <seed>           : [7] for the hashing function\n"
  "  -o <outfile prefix>    : [syng]\n"
  "  -T <threads>           : [4] number of threads\n" ;

int main (int argc, char *argv[])
{
  int w = 1023 ;
  int k = 16 ;
  int seed = 7 ;
  char *outPrefix = "syng-out" ;
  int nThread = 4 ;
  pthread_t *threads ;
  ThreadInfo *threadInfo ;
  SeqIO *sio ;
  OneSchema *schema = oneSchemaCreateFromText (syngSchemaText) ;
  OneFile *segseq, *seqsyn, *segFile, *linkFile ;
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
    else if (!strcmp (*argv, "-T") && argc > 1) { nThread = atoi(argv[1]) ; argc -=2 ; argv +=2 ; }
    else die ("unknown parameter %s\n%s", *argv, usage) ;

  Seqhash *sh = seqhashCreate (k, w+1, seed) ; // need the +1 here - not nice, could fix

  if (nThread < 1) die ("number of threads %d must be at least 1", nThread) ;
  threads = new (nThread, pthread_t) ;
  threadInfo = new0 (nThread, ThreadInfo) ;
  int i ;
  for (i = 0 ; i < nThread ; ++i)
    { threadInfo[i].sh = sh ;
      threadInfo[i].seq = arrayCreate (101<<20, char) ; // 101 Mb
      threadInfo[i].seqLen = arrayCreate (20000, I64) ;
      threadInfo[i].pos = arrayCreate (1<<20, I64) ;
      threadInfo[i].posLen = arrayCreate (20000, I64) ;
    }

  if (argc != 1) die ("need to give a read sequence file\n%s", usage) ;
  if (!(sio = seqIOopenRead (*argv, dna2indexConv, false)))
    die ("failed to open read sequence file %s", *argv) ;
  char *syncString = new0 (w+k+1, char) ;
  if (!(seqsyn = oneFileOpenWriteNew (fnameTag (outPrefix,"1readsyn"), schema, "readsyn", true, 1)))
      die ("failed to open %s.1readsyn to write", outPrefix) ;
  oneAddProvenance (seqsyn, "syng", "0.0", getCommandLine()) ;
  oneAddReference (seqsyn, *argv, 0) ; // I don't know what to put for count - maybe 0 is correct?
  Array readSync = arrayCreate(64,I64) ;
  Array readDir = arrayCreate(64,char) ;
  U64 nSeq = 0, totSeq = 0, totSync = 0 ;
  bool isDone = false ;
  while (!isDone) 
    { int i ;
      for (i = 0 ; i < nThread ; ++i)  // read 100Mb DNA per thread and find syncs in parallel
	{ ThreadInfo *ti = &threadInfo[i] ;
	  arrayMax(ti->seq) = 0 ;
	  arrayMax(ti->seqLen) = 0 ;
	  int seqStart = 0 ;
	  while (arrayMax(ti->seq) < 100<<20 && seqIOread (sio))
	    { array(ti->seqLen, arrayMax(ti->seqLen), I64) = sio->seqLen ;
	      array(ti->seq, seqStart+sio->seqLen, char) = 0 ;
	      memcpy (arrp(ti->seq, seqStart, char), sqioSeq(sio), sio->seqLen) ;
	      seqStart += sio->seqLen ;
	    }
	  if (!arrayMax(ti->seq)) // we are done
	    { nThread = i ;
	      isDone = true ;
	      break ;
	    }
	}
      for (i = 0 ; i < nThread ; ++i) // create threads
	pthread_create (&threads[i], 0, threadProcess, &threadInfo[i]) ;
      for (i = 0 ; i < nThread ; ++i)
	pthread_join (threads[i], 0) ; // wait for threads to complete
      for (i = 0 ; i < nThread ; ++i) //
	{ ThreadInfo *ti = threadInfo + i ;
	  I64 *posList = arrp(ti->pos, 0, I64), lastPos ;
	  char *seq = arrp(ti->seq, 0, char) ;
	  int iRead, iPos, iSync, iLastSync ;
	  for (iRead = 0 ; iRead < arrayMax(ti->posLen) ; ++iRead)
	    { I64 posLen = arr(ti->posLen, iRead, I64) ;
	      arrayMax(readSync) = 0 ; arrayMax(readDir) = 0 ;
	      for (iPos = 0 ; iPos < posLen ; ++iPos)
		{ I64 pos = posList[iPos] ;
		  int x = pos-1, y = pos+w+k ; while (seq[++x] == 3 - seq[--y]) ;
		  if (seq[x] < 3 - seq[y])
		    { array(readDir,iPos,char) = '+' ;
		      for (x = 0 ; x < w+k ; ++x) syncString[x] = acgt[seq[pos+x]] ;
		      dictAdd (syncDict, syncString, &iSync) ;
		      arrayp(syncs,iSync,Sync)->n++ ;
		      array(readSync,iPos,I64) = iSync ;
		    }
		  else
		    { array(readDir,iPos,char) = '-' ;
		      for (x = 0 ; x < w+k ; ++x) syncString[w+k-x-1] = tgca[seq[pos+x]] ;
		      dictAdd (syncDict, syncString, &iSync) ;
		      arrayp(syncs,iSync,Sync)->n++ ;
		      array(readSync,iPos,I64) = iSync ;
		      iSync = -iSync-1 ;
		    }
		  if (iPos > 1) // add the link
		    { int sa, sb, iLink ;
		      canonicaliseLink (iLastSync, iSync, &sa, &sb) ;
		      Link *link ;
		      if (hashAdd (linkHash, HASH_INT2(sa,sb), &iLink))
			{ link = arrayp(links,iLink,Link) ;
			  link->sa = sa ; link->sb = sb ;
			  link->overlap = w + k - (pos - lastPos) ;
			}
		      else
			{ link = arrp(links,iLink,Link) ;
			  if (w + k - (pos - lastPos) != link->overlap)
			    fprintf (stderr, "inconsistent overlap segment %d to %d: %d != %d\n",
				     sa, sb, w + k - (int)(pos - lastPos), link->overlap) ;
			}
		      ++link->n ; 
		    }
		  iLastSync = iSync ;
		  lastPos = pos ;
		} // iPos
	      if (posLen)
		{ oneWriteLine (seqsyn, 'S', posLen, arrp(readSync,0,I64)) ;
		  oneWriteLine (seqsyn, 'P', posLen, posList) ;
		  oneWriteLine (seqsyn, 'O', posLen, arrp(readDir,0,char)) ;
		  oneInt(seqsyn, 0) = nSeq ; oneWriteLine (seqsyn, 'R', 0, 0) ;
		}
	      nSeq++ ; totSeq += sio->seqLen ; totSync += arrayMax(readSync) ;
	      seq += arr(ti->seqLen, iRead, I64) ;
	      posList += posLen ;
	    } // read
	} // thread
    }

  printf ("read %" PRIu64 " sequences, total length %" PRIu64 ""
          ", yielding %" PRIu64 " instances of %d syncmers, average %.2f coverage\n", 
          nSeq, totSeq, totSync, arrayMax(syncs), totSync / (double)arrayMax(syncs)) ; 
  
  oneFileClose (seqsyn) ;
  arrayDestroy (readSync) ;
  arrayDestroy (readDir) ;

  if (!(segFile = oneFileOpenWriteNew (fnameTag (outPrefix,"1seg"), schema, "syncmer", true, 1)))
    die ("failed to open %s.1seg to write syncmer segs", outPrefix) ;
  oneAddProvenance (segFile, "syng", "0.0", getCommandLine()) ;
  for (i = 0 ; i < arrayMax(syncs) ; ++i)
    { Sync *s = arrp(syncs, i, Sync) ;
      oneInt(segFile,0) = w+k ; oneWriteLine (segFile, 'S', 0, 0) ;
      oneInt(segFile,0) = s->n ; oneWriteLine (segFile, 'K', 0, 0) ;
    }
  oneFileClose (segFile) ;

  if (!(linkFile = oneFileOpenWriteNew (fnameTag (outPrefix, "1link"), schema, "link", true, 1)))
    die ("failed to open %s.1link to write links", outPrefix) ;
  oneAddProvenance (segseq, "syng", "0.0", getCommandLine()) ;
  for (i = 0 ; i < arrayMax(links) ; ++i)
    { Link *l = arrp(links, i, Link) ;
      if (l->sa >= 0) { oneInt(linkFile,0) = l->sa ; oneChar(linkFile,1) = '+' ; }
      else { oneInt(linkFile,0) = -l->sa - 1 ; oneChar(linkFile,3) = '-' ; }
      if (l->sb >= 0) { oneInt(linkFile,2) = l->sb ; oneChar(linkFile,1) = '+' ; }
      else { oneInt(linkFile,2) = -l->sb - 1 ; oneChar(linkFile,3) = '-' ; }
      oneWriteLine (linkFile, 'L', 0, 0) ;
      oneInt(linkFile,0) = l->overlap ; oneWriteLine (linkFile, 'O', 0, 0) ;
      oneInt(linkFile,0) = l->n ; oneWriteLine (linkFile, 'R', 0, 0) ;
    }
  oneFileClose (linkFile) ;

  if (!(segseq = oneFileOpenWriteNew (fnameTag (outPrefix,"1syncseq"), schema, "syncseq", true, 1)))
    die ("failed to open %s.1syncseq to write syncmer sequences", outPrefix) ;
  oneAddProvenance (segseq, "syng", "0.0", getCommandLine()) ;
  fprintf (stderr, "w+k %d dictMax %d\n", w+k, dictMax(syncDict)) ;
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
