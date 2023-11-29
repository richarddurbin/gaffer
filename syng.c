/*  File: syng.c
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2023
 *-------------------------------------------------------------------
 * Description: syncmer-based graph assembler
 * Exported functions:
 * HISTORY:
 * Last edited: Nov 27 18:02 2023 (rd109)
 * Created: Thu May 18 11:57:13 2023 (rd109)
 *-------------------------------------------------------------------
 */

#define SYNC_VERSION "0.0"

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

static void writeParams (OneFile *of)
{
  oneInt(of,0) = params.k ;
  oneInt(of,1) = params.w ;
  oneInt(of,2) = params.seed ;
  oneWriteLine (of, 'h', 0, 0) ;
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
  char *outPrefix = "syng-out" ;
  int nThread = 4 ;
  pthread_t *threads ;
  ThreadInfo *threadInfo ;
  SeqIO *sio ;
  OneSchema *schema = oneSchemaCreateFromText (syngSchemaText) ;
  OneFile *segseq, *seqsyn, *segFile ;
  char *acgt = "acgt", *tgca = "tgca" ; // normal and complemented

  timeUpdate (0) ;

  syncs = arrayCreate (1<<20, Sync) ;
  syncDict = dictCreate (1<<20) ;

  params.k = PARAMS_K_DEFAULT ;
  params.w = PARAMS_W_DEFAULT ;
  params.seed = PARAMS_SEED_DEFAULT ;
  
  storeCommandLine (argc, argv) ;
  argc-- ; ++argv ;
  if (!argc) { printf ("%s",usage) ; exit (0) ; }
  while (argc > 0 && **argv == '-')
    if (!strcmp (*argv, "-w") && argc > 1) { params.w = atoi(argv[1]) ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (*argv, "-k") && argc > 1) { params.k = atoi(argv[1]) ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (*argv, "-seed") && argc > 1) { params.seed = atoi(argv[1]) ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (*argv, "-o") && argc > 1) { outPrefix = argv[1] ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (*argv, "-T") && argc > 1) { nThread = atoi(argv[1]) ; argc -=2 ; argv +=2 ; }
    else die ("unknown parameter %s\n%s", *argv, usage) ;

  fprintf (stderr, "k, w, seed are %d %d %d\n", params.k, params.w, params.seed) ;
  Seqhash *sh = seqhashCreate (params.k, params.w+1, params.seed) ; // need the +1 here - not nice, could fix

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
  char *syncString = new0 (params.w+params.k+1, char) ;
  if (!(seqsyn = oneFileOpenWriteNew (fnameTag (outPrefix,"1readsyn"), schema, "readsyn", true, 1)))
      die ("failed to open %s.1readsyn to write", outPrefix) ;
  oneAddProvenance (seqsyn, "syng", SYNC_VERSION, getCommandLine()) ;
  oneAddReference (seqsyn, *argv, 0) ; // I don't know what to put for count - maybe 0 is correct?
  writeParams (seqsyn) ;
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
	      totSeq += sio->seqLen ; 
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
	  I64 *posList = arrp(ti->pos, 0, I64) ;
	  char *seq = arrp(ti->seq, 0, char) ;
	  int iRead, iPos ;
    U64 iSync ;
	  for (iRead = 0 ; iRead < arrayMax(ti->posLen) ; ++iRead)
	    { I64 posLen = arr(ti->posLen, iRead, I64) ;
	      arrayMax(readSync) = 0 ; arrayMax(readDir) = 0 ;
	      for (iPos = 0 ; iPos < posLen ; ++iPos)
		{ I64 pos = posList[iPos] ;
		  int x = pos-1, y = pos+params.w+params.k ; while (seq[++x] == 3 - seq[--y]) ;
		  if (seq[x] < 3 - seq[y])
		    { array(readDir,iPos,char) = '+' ;
		      for (x = 0 ; x < params.w+params.k ; ++x) syncString[x] = acgt[seq[pos+x]] ;
		      dictAdd (syncDict, syncString, &iSync) ;
		      arrayp(syncs,iSync,Sync)->n++ ;
		      array(readSync,iPos,I64) = iSync ;
		    }
		  else
		    { array(readDir,iPos,char) = '-' ;
		      for (x = 0 ; x < params.w+params.k ; ++x)
			syncString[params.w+params.k-x-1] = tgca[seq[pos+x]] ;
		      dictAdd (syncDict, syncString, &iSync) ;
		      arrayp(syncs,iSync,Sync)->n++ ;
		      array(readSync,iPos,I64) = iSync ;
		    }
		} // iPos
	      if (posLen)
		{ oneWriteLine (seqsyn, 'S', posLen, arrp(readSync,0,I64)) ;
		  oneWriteLine (seqsyn, 'P', posLen, posList) ;
		  oneWriteLine (seqsyn, 'O', posLen, arrp(readDir,0,char)) ;
		  oneInt(seqsyn, 0) = nSeq ; oneWriteLine (seqsyn, 'R', 0, 0) ;
		}
	      nSeq++ ; totSync += arrayMax(readSync) ;
	      seq += arr(ti->seqLen, iRead, I64) ;
	      posList += posLen ;
	    } // read
	} // thread
    }

  printf ("read %" PRIu64 " sequences, total length %" PRIu64 ""
          ", yielding %" PRIu64 " instances of %" PRIu64 " syncmers, average %.2f coverage\n", 
          nSeq, totSeq, totSync, arrayMax(syncs), totSync / (double)arrayMax(syncs)) ; 
  
  oneFileClose (seqsyn) ;
  arrayDestroy (readSync) ;
  arrayDestroy (readDir) ;

  if (!(segFile = oneFileOpenWriteNew (fnameTag (outPrefix,"1seg"), schema, "syncmer", true, 1)))
    die ("failed to open %s.1seg to write syncmer segs", outPrefix) ;
  oneAddProvenance (segFile, "syng", SYNC_VERSION, getCommandLine()) ;
  for (i = 0 ; i < arrayMax(syncs) ; ++i)
    { Sync *s = arrp(syncs, i, Sync) ;
      oneInt(segFile,0) = params.w+params.k ; oneWriteLine (segFile, 'S', 0, 0) ;
      oneInt(segFile,0) = s->n ; oneWriteLine (segFile, 'K', 0, 0) ;
    }
  oneFileClose (segFile) ;

  if (!(segseq = oneFileOpenWriteNew (fnameTag (outPrefix,"1syncseq"), schema, "syncseq", true, 1)))
    die ("failed to open %s.1syncseq to write syncmer sequences", outPrefix) ;
  oneAddProvenance (segseq, "syng", SYNC_VERSION, getCommandLine()) ;
  writeParams (segseq) ;
  for (i = 0 ; i < dictMax(syncDict) ; ++i)
    oneWriteLine (segseq, 'S', params.w+params.k, dictName (syncDict, i)) ;
  oneFileClose (segseq) ;

  arrayDestroy (syncs) ;
  dictDestroy (syncDict) ;
	  
  fprintf (stderr, "total: ") ; timeTotal (stderr) ;
  return 0 ;
}

/*********************** end of file **********************/
