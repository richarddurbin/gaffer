/*  File: syng.c
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2023
 *-------------------------------------------------------------------
 * Description: syncmer-based graph assembler
 * Exported functions:
 * HISTORY:
 * Last edited: Sep 29 22:59 2024 (rd109)
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
  Seqhash  *sh ;
  KmerHash *kh ;        // read-only here
  Array     seq ;	// of char, input: concatenated sequences in index 0..3
  Array     seqLen ;	// of I64, input: per sequence
  Array     pos ;	// of I64, output: concatenated start positions of syncmers
  Array     posLen ;	// of I64, output: number of syncmers per sequence
  Array     sync ;      // of I64, if non-zero, found in kh; -ve if reverse direction
} ThreadInfo ;

static void *threadProcess (void* arg) // find the start positions of all the syncmers
{
  ThreadInfo *ti = (ThreadInfo*) arg ;
  int i ;
  I64 seqStart = 0 ;
  I64 iSync ;
  U64 *uBuf = new(ti->kh->plen,U64) ;

  arrayMax(ti->pos) = 0 ;
  arrayMax(ti->posLen) = 0 ;
  
  for (i = 0 ; i < arrayMax(ti->seqLen) ; ++i)
    { I64 seqLen = arr(ti->seqLen, i, I64) ;
      char *seq = arrp(ti->seq, seqStart, char) ;
      seqStart += seqLen ;
      int pos, posStart = arrayMax(ti->pos) ;
      SeqhashIterator *sit = syncmerIterator (ti->sh, seq, seqLen) ;
      while (syncmerNext (sit, 0, &pos, 0))
	{ I64 iPos = arrayMax(ti->pos) ;
	  array(ti->pos, iPos, I64) = pos ;
	  iSync = 0 ;
	  kmerHashFindThreadSafe (ti->kh, seq+pos, &iSync, uBuf) ;
	  array(ti->sync,iPos,I64) = iSync ; // will be 0 if not found
	}
      seqhashIteratorDestroy (sit) ;
      array(ti->posLen, i, I64) = arrayMax(ti->pos) - posStart ;
    }

  newFree (uBuf, ti->kh->plen, U64) ;
  return 0 ;
}

static void writeParams (OneFile *of)
{
  oneInt(of,0) = params.k ;
  oneInt(of,1) = params.w ;
  oneInt(of,2) = params.seed ;
  oneWriteLine (of, 'h', 0, 0) ;
}

static void checkParams (OneFile *of)
{
  if (oneInt(of,0) != params.k ||
      oneInt(of,1) != params.w ||
      oneInt(of,2) != params.seed)
    die ("hash parameters mismatch: (k,w,s) file (%d,%d,%d) != code (%d,%d,%d)",
	 oneInt(of,0), oneInt(of,1), oneInt(of,2), params.k, params.w, params.seed) ;
}

static char usage[] =
  "Usage: syng <options> <read sequence file>+\n"
  "  -w <syncmer length>    : [1023]\n"
  "  -k <smer length>       : [16] must be under 32\n"
  "  -seed <seed>           : [7] for the hashing function\n"
  "  -o <outfile prefix>    : [syngOut]\n"
  "  -T <threads>           : [8] number of threads\n"
  "  -A <sync file>         : add to this sync file\n" ;

int main (int argc, char *argv[])
{
  char       *outPrefix = "syngOut" ;
  char       *syncFileName = 0 ;
  int         nThread = 8 ;
  pthread_t  *threads ;
  ThreadInfo *threadInfo ;
  SeqIO      *sio ;
  OneSchema  *schema = oneSchemaCreateFromText (syngSchemaText) ;
  KmerHash   *kh ;

  timeUpdate (0) ;

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
    else if (!strcmp (*argv, "-A") && argc > 1) { syncFileName = argv[1] ; argc -= 2 ; argv +=2 ; }
    else die ("unknown parameter %s\n%s", *argv, usage) ;

  fprintf (stderr, "k, w, seed are %d %d %d\n", params.k, params.w, params.seed) ;
  Seqhash *sh = seqhashCreate (params.k, params.w+1, params.seed) ; // need the +1 here - not nice, could fix

  U64 nSeq = 0, totSeq = 0, totSync = 0 ;

  OneFile *ofSync ;
  if (syncFileName) // open it, check hash match and read in the syncmer hash and the node counts
    { OneFile *of = oneFileOpenRead (syncFileName, schema, "khash", 1) ;
      if (!of) die ("failed to open syncmer file %s to read", syncFileName) ;

      // read the syncmer hash parameters
      while (oneReadLine (of) && of->lineType != 'h') ;
      if (of->lineType == 'h') checkParams (of) ;
      else die ("sync file %s has no 'h' parameters record", syncFileName) ;

      // read in the kmerhash
      kh = kmerHashReadOneFile (of) ;
      if (kh->len != params.w + params.k)
	die ("syncmer len mismatch %d != %d + %d", kh->len, params.w, params.k) ;

      // read in the node counts 
      if (of->lineType != 'C' || oneLen(of) != kmerHashMax(kh))
	die ("failed to find expected C line in khash file") ;
      nodes = arrayCreate (kh->max+1, Node) ;
      I64 i, *counts = oneIntList(of) ;
      for (i = 1 ; i <= kmerHashMax(kh) ; ++i) totSync += (array(nodes,i,Node).n = counts[i-1]) ;

      printf ("read %llu nodes from %s\n", kmerHashMax(kh), syncFileName) ;

      ofSync = oneFileOpenWriteFrom (syncFileName, of, true, 1) ;
      if (!ofSync) die ("failed to reopen syncmer file %s to write", syncFileName) ;
      oneFileClose (of) ; // close the old file
      timeUpdate (stdout) ;
    }
  else
    { kh = kmerHashCreate (0, params.w + params.k) ; 
      nodes = arrayCreate (1<<20, Node) ;
      arrayp(nodes, 0, Node)->n = 0 ; // because first sync is at 1 not 0
      ofSync = oneFileOpenWriteNew (fnameTag (outPrefix,"1khash"), schema, "khash", true, 1) ;
      if (!ofSync) die ("failed to open %s.1khash to write syncmers", outPrefix) ;
    }
  oneAddProvenance (ofSync, "syng", SYNC_VERSION, getCommandLine()) ;
  writeParams (ofSync) ;
    
  OneFile *ofSeq = oneFileOpenWriteNew (fnameTag (outPrefix,"1syncseq"), schema, "syncseq", true, 1) ;
  if (!ofSeq) die ("failed to open %s.1nodeseq to write", outPrefix) ;
  oneAddProvenance (ofSeq, "syng", SYNC_VERSION, getCommandLine()) ;
  int i ;
  for (i = 0 ; i < argc ; ++i) oneAddReference (ofSeq, argv[i], i+1) ; // write the sources
  writeParams (ofSeq) ;
  
  Array readSync = arrayCreate(64,I64) ;
  Array readDir = arrayCreate(64,char) ;
  U64 nSeq0 = 0, totSeq0 = 0, totSync0 = totSync, syncMax0 = kmerHashMax(kh) ;
  
  if (nThread < 1) die ("number of threads %d must be at least 1", nThread) ;
  threads = new (nThread, pthread_t) ;
  threadInfo = new0 (nThread, ThreadInfo) ;
  for (i = 0 ; i < nThread ; ++i)
    { threadInfo[i].sh = sh ;
      threadInfo[i].kh = kh ;
      threadInfo[i].seq = arrayCreate (101<<20, char) ; // 101 Mb
      threadInfo[i].seqLen = arrayCreate (20000, I64) ;
      threadInfo[i].pos = arrayCreate (1<<20, I64) ;
      threadInfo[i].sync = arrayCreate (1<<20, I64) ;
      threadInfo[i].posLen = arrayCreate (20000, I64) ;
    }

  int sourceFile = 0 ;
  while (argc--)
    { ++sourceFile ;
      if (!(sio = seqIOopenRead (*argv, dna2indexConv, false)))
	die ("failed to open read sequence file %s", *argv) ;
      printf ("sequence file %d %s ", sourceFile, *argv) ; fflush (stdout) ;
      bool isDone = false ;
      while (!isDone) 
	{ for (i = 0 ; i < nThread ; ++i)  // read 100Mb DNA per thread and find nodes in parallel
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
		{ // nThread = i ;
		  isDone = true ;
		}
	    }
	  for (i = 0 ; i < nThread ; ++i) // create threads
	    pthread_create (&threads[i], 0, threadProcess, &threadInfo[i]) ;
	  for (i = 0 ; i < nThread ; ++i)
	    pthread_join (threads[i], 0) ; // wait for threads to complete
	  for (i = 0 ; i < nThread ; ++i) //
	    { ThreadInfo *ti = threadInfo + i ;
	      I64 *posList = arrp(ti->pos, 0, I64) ;
	      I64 *syncList = arrp(ti->sync, 0, I64) ;
	      char *seq = arrp(ti->seq, 0, char) ;
	      int iRead, iPos ;
	      for (iRead = 0 ; iRead < arrayMax(ti->posLen) ; ++iRead)
		{ I64 posLen = arr(ti->posLen, iRead, I64) ;
		  arrayMax(readSync) = 0 ; arrayMax(readDir) = 0 ;
		  for (iPos = 0 ; iPos < posLen ; ++iPos)
		    { I64 pos = posList[iPos] ;
		      I64 iSync = syncList[iPos] ;
		      if (!iSync) kmerHashAdd (kh, seq+pos, &iSync) ;
		      if (iSync > 0) array(readDir,iPos,char) = '+' ;
		      else { array(readDir,iPos,char) = '+' ; iSync = -iSync ; }
		      array(readSync,iPos,I64) = iSync ;
		      array(nodes,iSync,Node).n++ ;
		    } // iPos
		  if (posLen)
		    { oneWriteLine (ofSeq, 'S', posLen, arrp(readSync,0,I64)) ;
		      oneWriteLine (ofSeq, 'P', posLen, posList) ;
		      oneWriteLine (ofSeq, 'D', posLen, arrp(readDir,0,char)) ;
		      oneInt(ofSeq, 0) = sourceFile ;
		      oneInt(ofSeq, 1) = nSeq+1 ;
		      oneWriteLine (ofSeq, 'R', 0, 0) ;
		    }
		  nSeq++ ; totSync += arrayMax(readSync) ;
		  seq += arr(ti->seqLen, iRead, I64) ;
		  posList += posLen ;
		  syncList += posLen ;
		} // read
	    } // thread
	} // isDone: end of file
      seqIOclose (sio) ;
      printf ("had %llu sequences %llu bp, yielding %llu nodes including %llu extra syncmers\n",
	      nSeq-nSeq0, totSeq-totSeq0, totSync-totSync0, kmerHashMax(kh)-syncMax0) ;
      timeUpdate (stdout) ;
      nSeq0 = nSeq ; totSeq0 = totSeq ; totSync0 = totSync ; syncMax0 = kmerHashMax(kh) ;
      ++argv ;
    } // source file

  printf ("Total for this run %llu sequences, total length %llu\n", nSeq, totSeq) ;
  printf ("Overall total %llu instances of %llu syncmers, average %.2f coverage\n", 
	  totSync, kmerHashMax(kh), totSync / (double)kmerHashMax(kh)) ; 
  
  oneFileClose (ofSeq) ;
  arrayDestroy (readSync) ;
  arrayDestroy (readDir) ;

  kmerHashWriteOneFile (kh, ofSync) ;
  I64 *counts = new (kmerHashMax(kh), I64) ;
  for (i = 1 ; i <= kmerHashMax(kh)  ; ++i) counts[i-1] = arrp(nodes,i,Node)->n ;
  oneWriteLine (ofSync, 'C', kmerHashMax(kh), counts) ;
  newFree (counts, kmerHashMax(kh), I64) ;
  oneFileClose (ofSync) ;

  arrayDestroy (nodes) ;
  kmerHashDestroy (kh) ;
  seqhashDestroy (sh) ;
	  
  fprintf (stderr, "total: ") ; timeTotal (stderr) ;
  return 0 ;
}

/*********************** end of file **********************/
