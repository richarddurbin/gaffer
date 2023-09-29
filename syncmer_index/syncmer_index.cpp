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
	if (s==0) return "#";
	else if (s==1) return "$";

	s-=2;
	char strand=(s&1)?'-':'+';
	s>>=1;
	std::string d=std::to_string(s)+strand;
	return d;
}





// reverse complement a read. For syncmers, the operation consists of:
// 1. reversing the order of the syncmers
// 2. flipping the sense bit of each syncmer (XOR of least significant bit)
// NB this doesn't work if this Read is a collection as it will mess up the end markers
void reverseComplement( SyncmerRead& thisRead )
{
	reverse( thisRead.begin(), thisRead.end());
	for (auto& i : thisRead ) i ^= 1;
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
        std::istringstream iss(line);
		iss >> field;	
		if (field=="S")
		{ 
			numReads++;
			iss >> numEntries; // this contains num of syncmers, TBD check
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
	std::cerr << "read " << numSyncmers << " syncmers from " << numReads << " reads" << std::endl;


}


int main(int argc, char *argv[]){

if (argc!=2)
{
	std::cerr << "Usage: " << argv[0] << " oneFile" << std::endl;
	return -1;
}


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
uint_t alphabetSize=999999;

parseFile( argv[1], v );


//v.push_back(1);
//v.push_back(2);

#ifdef XXX
v.push_back(3);
v.push_back(4);
v.push_back(5);
v.push_back(1);
v.push_back(4);
v.push_back(5);
v.push_back(6);
v.push_back(7);
v.push_back(1);
v.push_back(6);
v.push_back(7);
v.push_back(8);
v.push_back(1);
v.push_back(0); // end of data marker
#endif

//v.push_back(45);
//v.push_back(46);
//v.push_back(46);
//v.push_back(1);
//v.push_back(0);

//printf("%lu\n",strlen(argv[1]));

#ifdef XXX
	// intput data
	if(argc>=2){
		//concatenate all strings s_1, s_2, .., s_d in s_1$s_2$..%s_d$#
		int i = 2, sum=0;
		for(; i<= argc; i++){
			sum += strlen((argv[i-1]))+1;
		}
		n = sum+1;
		Text = malloc(n*sizeof(unsigned char));
		sum=0;
		for(i=2; i<= argc; i++){
			sscanf(argv[i-1], "%s", &Text[sum]);
			sum += strlen((argv[i-1]))+1;
			Text[sum-1]=1;//separator
		}
		Text[n-1]=0;
		printf("N = %d\n", n);
	}
	else{
		fprintf(stderr, "Please, insert at least one string.\n");
		exit(-1);
	}

	int i, j;
	printf("T^{cat} = ");
	for(i=0;i<n-1;i++){
		if(Text[i]==1) printf("$");
		else printf("%c", Text[i]);
	}
	printf("#\n");
#endif

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
	n=v.size();
	std::cout << "Read " << n << " characters in all" <<std::endl;

	// allocate
	uint_t *SA = (uint_t *)malloc(n * sizeof(uint_t));
	int_t *LCP = (int_t *)malloc(n * sizeof(int_t));
	int_da *DA = (int_da *)malloc(n * sizeof(int_da));
	printf("about to run\n");
	// sort
	gsacak_int((int_text*)&v[0], (uint_t*)SA, LCP, DA, n, alphabetSize);
//	gsacak_int((int_text*)&v[0], (uint_t*)SA, NULL, NULL, n, 99);

if (n>100) n=100; // avoid printing too much;
std::cout << "i\tSA\tDA\tLCP\tBWT\tstring\n" << std::endl;
	for(int i = 0; i < n; ++i) {
	//    char j = (SA[i])? Text[SA[i]-1]:'#';
	//    if(j==1) j = '$';
	    int_text j = (SA[i])? v[SA[i]-1]:999;
		std::cout << i << "\t" 
					<< SA[i] << "\t" 
					<< DA[i] << "\t" 
					<< LCP[i] << "\t"
					<< j << " " << decode(j) << " "
					<< v[i] << " " << decode(v[i]) << std::endl;
//	    printf("%d\t%d\t%d\t%d\t%d\t%d\n",i, SA[i], DA[i], LCP[i], j, v[i]);
//	    for(j = SA[i]; j < n; ++j) {
//		if(Text[j]==1) printf("$");
//		else printf("%c", Text[j]);
//	    }
//	    printf("#\n");
	}
#ifdef XXX
	// output
	printf("i\tSA\tDA\tLCP\tBWT\tsuffixes\n");
	for(i = 0; i < n; ++i) {
	    char j = (SA[i])? Text[SA[i]-1]:'#';
	    if(j==1) j = '$';
	    printf("%d\t%d\t%d\t%d\t%c\t",i, SA[i], DA[i], LCP[i], j);
	    for(j = SA[i]; j < n; ++j) {
		if(Text[j]==1) printf("$");
		else printf("%c", Text[j]);
	    }
	    printf("#\n");
	}
#endif
	printf("about to free\n");
	// deallocate
	free(SA);
	free(DA);
	free(LCP);
//	free(Text);

return 0;
}

