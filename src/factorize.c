
/*
 *All the functions in this file deals with factorizing using erthostenes sieve 
 *Functions:
  build_arr,factorize,group,print_table
 */

#include<stdio.h>
#include<stdlib.h>
#include<gmp.h>
#include<math.h>
#include<string.h>
#include <pthread.h>

/*For some common definations*/
#include "qfs.h"

static inline void repete (register struct node **ptr, const int offset, const int step, const int pr,
			  const int idx, const int mask);

static inline void trail (struct node *ptr, const int pr, const int idx, const int mask);

static void eratosthenes (struct node **start, int *prime, int B);

static void *threaded_factorize (void *param);


extern unsigned long int sieve_offset;
extern int nthread;		//defined in main.c contains the number of threads to be created


//Structure to pass data to each thread

/*thread_idx is a unique number assigned to each thread manually 
and its not pthread_t or value assigned by pthread_create*/

struct thread_data{
	unsigned int thread_idx;
	struct node **start;
	int *prime;
	int B;
};

/*
 *NOTE:Different threads are handling differnt elements of struct node array
 *     so there can be no race condition => no mutex required
 */


/*
 *initializes and allocates memory for the required members of the structure node......
 *creats the array and initializes q_x to (seed+x)^2-qfs_num
 *vector node->v to 00000
 *returns a pointer to the node array
 */

struct node **build_arr (mpz_t offset, int v_size)
{
	struct node **p;
	register unsigned int i;
	mpz_t seed;
	
	//allocation for the array of pointers of main seive array
	p = malloc (sizeof(struct node *) * (sieve_offset));
	
	//set seed
	mpz_init (seed);
	mpz_sqrt (seed, qfs_num);
	
	//Allocate and initialize the various elements of struct
	for (i=0 ;i<sieve_offset ;i++){
		register struct node *t;
		//allocate memory for structure
		t = malloc (sizeof(struct node));
		mpz_init (t->q_x);
		mpz_init (t->m);
		
		//set m
		mpz_set (t->m, offset);
		mpz_add_ui ((t)->m, (t)->m, i);
		
		//calculate q_x
		mpz_add (t->q_x, seed, t->m);
		mpz_mul (t->q_x, t->q_x, t->q_x);
		mpz_sub (t->q_x, t->q_x, qfs_num);
		
		//allocate for the vector v and set it to zero
		(t)->v = calloc(v_size, sizeof(int));
		*(p+i)=t;
	}
	
	mpz_clear (seed);
	return p;
}


/*
 *Factorizes the numbers in the list using sieve of Eratosthenes
 *start=>pointer to the node corresponding to x=-M
 *prime=>pointer to the factor base
 *B=>no of elements in the factor base
 */

void factorize (struct node **start, int *prime, int B)
{
	int i=0;
	pthread_t *thread_id;			//array of thread id's
	struct thread_data **thread_args;	//thread_arguments
		
	
	//run eratosthenes algorithm 
	eratosthenes (start, prime, B);
	
	//allocate space to store the thread id's
	thread_id = malloc (sizeof (pthread_t) * nthread);
	//allocate space for the pointers of thread data
	thread_args = malloc (sizeof (struct thread_data *) * nthread);
	
	
	/*Create as many threads as requested by the user*/
	for (i=0 ; i < nthread ;i++) {
		(*(thread_args + i)) = malloc (sizeof(struct thread_data));
		
		(*(thread_args +i))->thread_idx = i;
		(*(thread_args +i))->start = start;
		(*(thread_args +i))->prime = prime;
		(*(thread_args +i))->B = B;
		
		//create thread
		pthread_create (&thread_id[i], NULL, &threaded_factorize, *(thread_args + i));
	}
	
	
	/*join all the threads*/
	for (i=0 ; i<nthread ;i++) {
		pthread_join (thread_id[i], NULL);
		//free the thread data structures
		free (*(thread_args+i));
	}
	
	
	free (thread_id);
	free (thread_args);
	//done
}


/*
 *This function implements eratosthenes sieve
 */

static void eratosthenes (struct node **start, int *prime, int B)
{
	int i=0,j=0;
	int level,step;
	int t1;
	
	//current word accessing in v
	int idx=0;
	//masks the current bit accessed in a word
	int mask=1;

	
	struct node **temp;
	
	for (i=0 ; i < B ; i++){
		temp = start;
		
		/*if sieve_offset is smaller than the largest prime in prime set
		 *then we may have a segmentation fault
		 *
		 *i!=0 is reqd bcoz say sieve_offset=10000 and *(prime+i)=-1 (treated as unsigned int in comparision)
		 *in the above case the we get no factors
		 */
		if (sieve_offset < *(prime+i) && i!=0 )
			break;
		//On wordsize boundaries execute the loop
		if ((i % WORD_SIZE == 0) && (i != 0) ) {
			//icrement to point to next word
			idx++;
			//reinitialize mask to point to bit 0
			mask=1;
		}
		
		if (!i) {
			/*if a number is negative make it positive and update vector*/
			t1 = sieve_offset;
			while (t1--) {
				if (mpz_sgn ((*temp)->q_x) < 0) {
					mpz_neg ((*temp)->q_x, (*temp)->q_x);		//update q_x
					(*temp)->v[idx] = (*temp)->v[idx] ^ mask;	//update v
				}
				temp++;
			}
		}
		else {
			level=0;
			for (j=0; j < *(prime+i) ;j++){
				//if number is not divisable continue
				if (!mpz_divisible_ui_p((*(temp+j))->q_x, *(prime+i)))
					continue;
				//the number is divisble
				else {
					/*level is required to determine the number of steps or
					 *elements to be skipped for the erathosthenes sieve
					 */
					level=0;
					
					while ((!level) || mpz_divisible_ui_p ((*(temp+j))->q_x, *(prime+i)) ) {
						
						level++;
						step = pow (*(prime+i),level);
						
						//divide the rest of the series
						repete (temp ,j ,step ,*(prime+i),idx, mask);
					}
				}
			}
		}
		mask = mask << 1;	//update mask for next prime
	}
}




static void *threaded_factorize (void *param)
{
	unsigned long int i;
	unsigned long int j;
	int mask,idx;
	
	struct thread_data *td;
	
	td = (struct thread_data *) param;
		
	for (i=td->thread_idx ; i < sieve_offset ; i += nthread){
		//value of ptr is const and not the value pointed by it
		register struct node * const ptr = *(td->start + i);
		
		/*j=0 corresponds to -1 so skip it*/
		idx = 0;
		//mask =2 as we are skipping -1 corresponding to mask=1
		mask = 2;
		
		for (j=1 ; j < td->B ; j++) {
			
			//we need not check for j!=0 as in eratosthenes for obvoius reasons
			if (j % WORD_SIZE == 0 ) {
				idx++;
				mask=1;
			}
			trail (ptr, *(td->prime + j), idx,mask);
			mask = mask << 1;
		}
	}
	return NULL;
}


/*
 *Code to do the trail division
 *ptr ---pointing to element corresponding to -M
 *offset --from above element
 *pr --present prime we r dealing
 */

static inline void trail (struct node * const ptr,  const int pr,
			  const int idx,  const int mask)
{
	while( mpz_divisible_ui_p (ptr->q_x, pr)) {
		//update q_x
		mpz_divexact_ui ( (ptr)->q_x,  (ptr)->q_x  ,pr);
		//update v
		(ptr)->v[idx] = (ptr)->v[idx] ^ mask;
	}
}




/*
 *Divides all the q_x of struct node in array p and updates the structure members
 *ptr =>  address of 0th element of main array
 *offset => position of starting point
 *step => no of elements to be skipped
 *M => sieving limit
 *pr => present base number we r dealing with...........
 */

static inline void repete (register struct node **ptr,const int offset,const int step,const int pr,
			   const int idx,const int mask)
{
	struct node **limit;
	
	limit = ptr + sieve_offset;
	ptr = ptr + offset;
	
	while (  ptr < limit ) {
		//update q_x
		mpz_divexact_ui ((*ptr)->q_x, (*ptr)->q_x, pr);
		//update q
		(*ptr)->v[idx] = (*ptr)->v[idx] ^ mask;
		ptr = ptr + step;
	}
}


#ifdef POIU_DEBUG
/*
 *This function is very useful while debugging and unnecessary while factorizing
 *Stores the table to specified file
 *num is the number of elements to print
 */

void print_table (struct node **ptr, int *prime, char *s, int num, int B)
{
	int i,j,sm_cnt=0;
	FILE *fp;
	if (! (fp = fopen (s, "w"))){
		fprintf (stderr, "\nERROR:Cannot write to the file %s\n", s);
		exit (0);
	}
	
	gmp_fprintf (fp, "Primes Used:\n");
	
	for (i=0 ; i<B ; i++){
		gmp_fprintf (fp, "%10d ", *(prime+i));
		if(i % 8 ==0) fprintf(fp, "\n");
	}
	
	gmp_fprintf (fp, "\n\n\n     m       Q(x)             vectors\n");
	
	for(i=0 ; (i<num) && ((*(ptr+i))!=NULL) ;i++){
		
		gmp_fprintf (fp, "%0Zd %0Zd\t\t", (*(ptr+i))->m, (*(ptr+i))->q_x);
		if(!mpz_cmp_ui ((*(ptr+i))->q_x ,1))
			sm_cnt++;
		
		for(j=0;j < B/WORD_SIZE+1 ;j++){
			gmp_fprintf (fp, "%08X ", (*(ptr+i))->v[j]);
		}
		gmp_fprintf (fp, "\n");
	}
	
	gmp_fprintf (fp, "Total no of smooth numbers %d", sm_cnt);
	fclose (fp);
}

#endif

