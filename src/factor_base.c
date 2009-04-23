/*
 *THE Code in this file deals about building the factor base for the quadratic sieve.
 *Functions:
 *	build_factor_base	
 */

#include<stdio.h>
#include<stdlib.h>
#include<gmp.h>
#include<math.h>
#include<string.h>

/*For some common definations*/
#include "qfs.h"

//Maxprimes is the number of primes in the file which should be present in the PWD
#define MAXPRIMES 18384


static int * compress(int *p, int count);

/*
 *Builds the factor base for the sieve
 *returns the no of elements in the factor base
 */

unsigned int build_factor_base(int **p,  int prime_bound,  int prime_multiple)
{
	int i=0,count=0,temp;
	int limit;
	mpz_t tmp;
	FILE *fp;

	/* just to make full use of all the bits in the vector field
	 * of struct node 
	 */
	limit = (prime_multiple * prime_bound) + WORD_SIZE - (prime_multiple * prime_bound) % WORD_SIZE;
	
	(*p) = malloc (sizeof(int) * limit);
	bzero (*p, sizeof(int) * prime_multiple * prime_bound);
	
	mpz_init (tmp);
	
	/*file "primes.txt" should contain lots of primes and should be in pwd */
	if ( (fp = fopen ("primes.txt", "r")) == NULL) {
		fprintf(stderr, "\n\nERROR:(build_factor_base) primes.txt not found\n\n");
		exit(1);
	}

	
	for (i=0; i < MAXPRIMES; i++) {
		if (fscanf (fp, "%d", &temp) == EOF) {
			printf ("\nERROR:(build_factor_base) Could not find MAXPRIMES"
					"primes in primes.txt\n\n");
			exit (1);
		}
		
		//set tmp=temp
		mpz_set_si (tmp,temp);
		
		/*add the legendre symbols to the factor base
		 *i.e if (n/p)=1  ==> add
		 *else skip
		 */
		if ( mpz_legendre(qfs_num,tmp)==1  || temp < 3 ){
			*(*p+count++) = temp;
		}
		if(count == limit)
			break;
	}
	
	if(count == MAXPRIMES)
		printf("\n\nWarning:Maxprimes Reached add more primes to prime.txt\n\n");
	
	fclose(fp);
	
	(*p)=compress(*p,count);
	mpz_clear(tmp);
	
	return count;
}




/*
 *Only to get rid of additional memory allocated  to *p in the build factor base
 */
static int * compress (int *p,int count)
{
	int *q,i;
	
	//allocate to accomdate all primes
	q=malloc(sizeof(int)*count);
	
	for (i=0;i<count;i++) {
		*(q+i)=*(p+i);
	}
	
	//free the array initially allocated
	free(p);
	return q;
}

/*
 *TODO:Think of some other way to get the same funtionality or better
 *i feel the design of the functions is quiet messed up
*/

