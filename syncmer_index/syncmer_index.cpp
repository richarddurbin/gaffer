#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "gsacak.h"
#include "experiments/external/malloc_count/malloc_count.h" //memory counter

#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <cassert>
#include <algorithm>



extern "C"
{
//#include "../utils.h"
#include "../ONElib.h"
}


//extern "C" oneFileOpenRead(char const*, OneSchema*, char const*, int);

#ifndef DEBUG
	#define DEBUG 0
#endif

// Each syncmer is encoded as (i<<1)|sense + 2, where i is the syncmer index from the ONE file
// and sense is the corresponding strand: '+' equals 0 and '-' equals 1
// The "+2" is necessary because the suffix array code reserves letters 0 and 1
// for "end of entire dataset" and "end of read" respectively


const int_text END_READ=1;
const int_text END_DATA=0;

typedef std::vector<int_text> SyncmerRead;

void addRead( SyncmerRead& v, const SyncmerRead& thisRead )
{
	v.insert( v.end(), thisRead.begin(), thisRead.end() );
	v.push_back(END_READ);
}

std::string decode( int_text s )
{
	if (s==END_DATA) return "#";
	else if (s==END_READ) return "$";

	s-=2;
	char strand=(s&1)?'-':'+';
	s>>=1;
	std::string d=std::to_string(s)+strand;
	return d;
}

struct SyncmerIndex
{
	SyncmerIndex(size_t n) : SA(n), LCP(n), DA(n) {}
	void save( const char* prefix );
	void load( const char* prefix );


	std::vector<uint_t> SA;
	std::vector<int_t> LCP;
	std::vector<int_da> DA;
};

//gsacak_int((int_text*)&v[0], (uint_t*)&s.SA[0], (int_t*)&s.LCP[0], (int_da*)&s.DA[0], n, alphabetSize);

void SyncmerIndex::save( const char* prefix )
{
	FILE* pOut;
	std::string prefixString(prefix), name;
	
	name=prefixString+".sa";
	std::cerr << "Writing index component to " << name << std::endl; 
	pOut=fopen(name.c_str(),"wb");
	assert(pOut!=NULL);
	assert(fwrite(&SA[0],sizeof(uint_t),SA.size(),pOut)==SA.size());
	fclose(pOut);

	name=prefixString+".lcp";
	std::cerr << "Writing index component to " << name << std::endl; 
	pOut=fopen(name.c_str(),"wb");
	assert(pOut!=NULL);
	assert(fwrite(&LCP[0],sizeof(int_t),LCP.size(),pOut)==LCP.size());
	fclose(pOut);

	name=prefixString+".da";
	std::cerr << "Writing index component to " << name << std::endl; 
	pOut=fopen(name.c_str(),"wb");
	assert(pOut!=NULL);
	assert(fwrite(&DA[0],sizeof(int_da),DA.size(),pOut)==DA.size());
	fclose(pOut);


}

void SyncmerIndex::load( const char* prefix )
{
	FILE* pIn;
	std::string prefixString(prefix), name;
	size_t s;

	name=prefixString+".sa";
	std::cerr << "Reading index component from " << name << std::endl; 
	pIn=fopen(name.c_str(),"rb");
	assert(pIn!=NULL);
	fseek(pIn, 0L, SEEK_END); s=ftell(pIn); rewind(pIn);
	SA.resize(s/sizeof(uint_t));
	assert(fread(&SA[0],sizeof(uint_t),SA.size(),pIn)==SA.size());
	fclose(pIn);

	name=prefixString+".lcp";
	std::cerr << "Reading index component from " << name << std::endl; 
	pIn=fopen(name.c_str(),"rb");
	assert(pIn!=NULL);
	fseek(pIn, 0L, SEEK_END); s=ftell(pIn); rewind(pIn);
	LCP.resize(s/sizeof(int_t));
	assert(fread(&LCP[0],sizeof(int_t),LCP.size(),pIn)==LCP.size());
	fclose(pIn);

	name=prefixString+".da";
	std::cerr << "Reading index component from " << name << std::endl; 
	pIn=fopen(name.c_str(),"rb");
	assert(pIn!=NULL);
	fseek(pIn, 0L, SEEK_END); s=ftell(pIn); rewind(pIn);
	DA.resize(s/sizeof(int_da));
	assert(fread(&DA[0],sizeof(int_da),DA.size(),pIn)==DA.size());
	fclose(pIn);


}




// reverse complement a read. For syncmers, the operation consists of:
// 1. reversing the order of the syncmers
// 2. flipping the sense bit of each syncmer (XOR of least significant bit)
// NB this doesn't work if thisRead is a collection as it will mess up the end markers
void reverseComplement( SyncmerRead& thisRead )
{
	reverse( thisRead.begin(), thisRead.end());
	for (auto& i : thisRead ) i ^= 1;
}


uint_t parseFileONE( const char* fileName, std::vector<int_text>& v )
{
	OneSchema* vs=NULL;
	char* fileType=NULL;

	OneFile *in = oneFileOpenRead(fileName, vs, fileType, 1);
   	if (in == NULL) 
 	{ fprintf(stderr,"Can't open sequence file %s to read\n",fileName);
	   exit (1);
	}
	int numEntries=0;
	int numSyncmers=0;
	int numReads=0;
	SyncmerRead thisRead;
	char* strand;
	int_text maxSyncmer=0;

	while (oneReadLine(in))
	{     
	//	std::cout << in->lineType << std::endl;
		if (in->lineType == 'S')
       	{ 
			numReads++;
			numEntries=oneLen(in);
			thisRead.clear();
			I64* r=oneIntList(in);
			for (int i(0);i<numEntries;i++,r++ )
			{
				thisRead.push_back(*r);
				numSyncmers++;
			}
       	//  reverseComplement(oneString(in), oneLen(in)); // user-provided, assume acts in place
        	// oneWriteLine(out, 'S', oneLen(in), oneString(in));
       	}
     	else if (in->lineType == 'O')
		{
			strand=oneString(in);
			for (std::size_t i=0;i<thisRead.size();i++)
			{
				if (thisRead[i]>maxSyncmer) maxSyncmer=thisRead[i];
				assert( (strand[i]=='+')||(strand[i]=='-'));
				thisRead[i]<<=1;
				thisRead[i]|=(1*(strand[i]=='-'));
				thisRead[i]+=2;
			}
			addRead( v, thisRead);
			reverseComplement( thisRead );
			addRead( v, thisRead);
		}

	}

	v.push_back(END_DATA); // end of data marker
	std::cout << "read " << numSyncmers << " syncmers from " << numReads << " reads" << std::endl;

	oneFileClose(in);

	maxSyncmer<<=1;
	maxSyncmer|=1;
	maxSyncmer+=3; // I think this is END_READ and END_DATA plus 1 because syncmers start at 0
	std::cout << "alphabet size is " << maxSyncmer << std::endl;
	return maxSyncmer;
}


void parseFile( const char* fileName, std::vector<int_text>& v )
{
	std::ifstream inputFile(fileName);
    if (!inputFile) {
        std::cerr << "Failed to open the file." << std::endl;
        exit (-1);
    }

	// very basic parse of text version of one format
	std::string field;
	int_text syncmer;
    std::string line, strand, numEntries;
	int numSyncmers=0;
	int numReads=0;
	SyncmerRead thisRead;
    while (std::getline(inputFile, line)) 
	{
		field.clear(); // otherwise a blank line can trigger an erroneous extra read
        std::istringstream iss(line);
		iss >> field;	
		if (field=="S")
		{ 
			numReads++;
			iss >> numEntries; // this contains num of syncmers
			thisRead.clear();	
        	while (iss >> field) 
			{
				try { thisRead.push_back(std::stoi(field)); }
				catch (const std::invalid_argument& e) 
				{ std::cerr << "Invalid integer: " << field << std::endl; } 
				catch (const std::out_of_range& e) 
				{ std::cerr << "Integer out of range: " << field << std::endl; }
				numSyncmers++;
			}
//			std::cout << field << " " << thisRead.size() << std::endl;
			assert(std::stoi(numEntries)==thisRead.size());
//			v.back()+=2; // 0 and 1 are reserved letters so shift everything up by 2
        }
		else if (field=="O")
		{	
			iss >> numEntries;
			iss >> strand;
			assert(std::stoi(numEntries)==thisRead.size());
			assert(std::stoi(numEntries)==strand.size());
			for (std::size_t i=0;i<thisRead.size();i++)
			{
				assert( (strand[i]=='+')||(strand[i]=='-'));
				thisRead[i]<<=1;
				thisRead[i]|=(1*(strand[i]=='-'));
				thisRead[i]+=2;
			}
			addRead( v, thisRead);
			reverseComplement( thisRead );
			addRead( v, thisRead);
        }
//		v.push_back(1); // end of string marker
    }
	v.push_back(END_DATA); // end of data marker
	std::cout << "read " << numSyncmers << " syncmers from " << numReads << " reads" << std::endl;


}


void scanIndex( const SyncmerIndex& s, const SyncmerRead& v )
{
//if (n>100) n=100; // avoid printing too much;
	std::cout << "# i\tSA\tDA\tLCP\tBWT\tsorted_char\tnext_char\tlcp_char\textend_char\tstring_char" << std::endl;
	std::cout << "# i\t- index of entry" << std::endl;
	std::cout << "# SA[i]\t - suffix array (position of i-th smallest suffix)" << std::endl;
	std::cout << "# DA[i]\t - document array (which read i-th smallest suffix is from)" << std::endl;
	std::cout << "# LCP[i]\t - LCP array (common starting chars between i-th and (i-1)th smallest suffix)" << std::endl;
	std::cout << "# first[i]\t - starting char of i-th smallest suffix" << std::endl;
	std::cout << "# second[i]\t - second char of i-th smallest suffix" << std::endl;
	std::cout << "# lcp_next[i]\t - first char in i-th smallest that differs from (i+1)st - will be $ for maximal match" << std::endl;
	std::cout << "# lcp_prev[i]\t - first char in i-th smallest that differs from (i-1)th" << std::endl;
	// fields in order

	size_t n(s.SA.size());

	for(int i = 0; i < n; ++i)
	{
	//    char j = (SA[i])? Text[SA[i]-1]:'#';
	//    if(j==1) j = '$';
	  int_text j = (s.SA[i])? v[s.SA[i]-1]:END_DATA;
	  int_text k = (s.SA[i]+1==n)?END_DATA:v[s.SA[i]+1];
	  int_text l;
	  std::string dl;
	  if (i==n-1)
	    l=END_DATA;
	  else if (s.SA[i]+s.LCP[i+1]==n)
	    l=END_DATA;
	  else
	    l=v[s.SA[i]+s.LCP[i+1]];
	  dl=decode(l); if (dl!="$" && (s.LCP[i+1]!=0)) dl="ZZ"+dl;

	  int_text ll;
	  if (s.SA[i]+s.LCP[i]==n)
	    ll=END_DATA;
	  else
	    ll=v[s.SA[i]+s.LCP[i]];
	  //  dl=decode(l); if (dl!="$" && (LCP[i+1]!=0)) dl="ZZ"+dl;
	  
#ifdef XXX	  
	  std::cout << i << "\t" 
		    << SA[i] << "\t" 
		    << DA[i] << "\t" 
		    << LCP[i] << "\t"
		    << j << " " << decode(j) << "\t"
		    << v[SA[i]] << " " << decode(v[SA[i]]) << "\t"
		    << k << " " << decode(k) << "\t"
		    << l << " " << dl << "\t"
		    << ll << " " << decode(ll) << "\t"
		    << v[i] << " " << decode(v[i]) << std::endl;
#endif
	  
		// more brief output - just show decoded syncmers and don't show v[i]

	  std::cout << i << "\t" 
		    << s.SA[i] << "\t" 
		    << s.DA[i] << "\t" 
		    << s.LCP[i] << "\t"
		    << decode(j) << "\t"
		    << decode(v[s.SA[i]]) << "\t"
		    << decode(k) << "\t"
		    << dl << "\t"
		    << decode(ll) << "\t"
		    << std::endl;


	}


}



int main(int argc, char *argv[])
{

	if ((argc!=4))
	{
		std::cerr << "Usage: " << argv[0] << " build seqFile.1 outPrefix: build index from oneFile, store index in outPrefix.*" << std::endl;
		std::cerr << "Usage: " << argv[0] << " scan seqFile.1 outPrefix: load index from outPrefix.*, output ASCII analysis" << std::endl;
		return -1;
	}


//else if (argc==3)
//{
//
//}
	#if DEBUG
		printf("sizeof(int_t) = %zu bytes\n", sizeof(int_t));
	#endif

	unsigned char *Text;
	uint_t n=0;

// gsacak(s, SA, NULL, NULL, n) //computes only SA
// gsacak(s, SA, LCP,  NULL, n) //computes SA and LCP
// gsacak(s, SA, NULL, DA, n)   //computes SA and DA
// gsacak(s, SA, LCP,  DA, n)   //computes SA, LCP and DA

	std::vector<int_text> v;
	uint_t alphabetSize=parseFileONE( argv[2], v );
//	alphabetSize+=2;
/** @brief Computes the suffix array SA (LCP, DA) of T^cat in s[0..n-1]
 *
 *  @param s		input concatenated string, using separators s[i]=1 and with s[n-1]=0
 *  @param SA		Suffix array
 *  @param LCP	LCP array
 *  @param DA		Document array
 *  @param n		String length
 *  @param k    alphabet size+2 (0 and 1 are reserved)
 *
 *  @return depth of the recursive calls.
int gsacak_int(int_text *s, uint_t *SA, int_t *LCP, int_da *DA, uint_t n, uint_t k);

 */

	if (strcmp(argv[1],"build")==0)
	{	
	n=v.size();
	std::cout << "# Read " << n << " characters in all" <<std::endl;

	// allocate
	SyncmerIndex s(n);

	std::cout << "# About to build index" << std::endl;
	// sort
	gsacak_int((int_text*)&v[0], (uint_t*)&s.SA[0], (int_t*)&s.LCP[0], (int_da*)&s.DA[0], n, alphabetSize);

	s.save(argv[3]);


//	scanIndex(s, v);
	}
	else if (strcmp(argv[1],"scan")==0)
	{
		SyncmerIndex s(1);
		s.load(argv[3]);
		scanIndex(s, v);
	}
	else
	{
		std::cerr << "could not understand command " << argv[1] << std::endl;
		return -1;
	}

	return 0;
}

