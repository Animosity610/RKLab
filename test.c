#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <memory.h>
#include <ctype.h>

/* a large prime for RK hash (BIG_PRIME*256 does not overflow)*/
long long BIG_PRIME = 5003943032159437; 

/* constants used for printing debug information */
const int PRINT_RK_HASH = 5;
const int PRINT_BLOOM_BITS = 160;

/* modulo addition */
long long
madd(long long a, long long b)
{
	return ((a+b)>BIG_PRIME?(a+b-BIG_PRIME):(a+b));
}

/* modulo substraction */
long long
mdel(long long a, long long b)
{
	return ((a>b)?(a-b):(a+BIG_PRIME-b));
}

/* modulo multiplication*/
long long
mmul(long long a, long long b)
{
	return ((a*b) % BIG_PRIME);
}
int doublecheck(const char *first, int slen, const char *second)
{
	int i = 0;
	for (i = 0; i < slen; i++) if (first[i] != second[i]) return 0; // checks if 2 strings of the same length are equal
	return 1;
}

int
rabin_karp_match(const char *ps,	/* the query string */
								 int k, 					/* the length of the query string */
								 const char *ts,	/* the document string (Y) */ 
								 int n						/* the length of the document Y */ )
{
	long long QVALUE = 0; // variable for hash value of ps
	long long CVALUE[(n-k)+1];	//this is the array of hash values of n-k+1 substrings of length k of ts
	long long PVALUE[k];	//this is the value of the powers of 256
	PVALUE[0]=1;	//256^0 = 1
	CVALUE[0]=ts[k-1];
	QVALUE=ps[k-1];

	int i;

	for (i=1; i < k;i++){
		PVALUE[i]=mmul(PVALUE[i-1],256);
		QVALUE = madd(QVALUE,mmul(PVALUE[i],ps[k-i-1]));
		CVALUE[0] = madd(CVALUE[0],mmul(PVALUE[i],ts[k-i-1]));
	}

	for (i=1;i<n-k+1;i++){
		CVALUE[i] = madd(ts[(i+k)-1],mmul(256,mdel(CVALUE[i-1],mmul(PVALUE[k-1],ts[i-1]))));
	}
    printf("%lli\n",QVALUE);
    printf("%lli\n",CVALUE[0]);
	for (i=0;i<n-k+1;i++){
		if (QVALUE == CVALUE[i] && doublecheck(ps,k,&ts[i]))
			return 1;
	}
	return 0;
}


int main(){
    char s1[]="Long";
    char s2[]="Long is Long";
    
    int a = rabin_karp_match(s1,4,s2,12);
    printf("%d\n",a);
    return 0;
}