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

static inline void repete (struct node **ptr, const int offset, const int step, const int pr,
			  const int idx, const int mask);
			  
static inline void trail (struct node **ptr, const int offset, const int pr,
			const int idx, const int mask);
			
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

/*A mutex to prevent race condition*/
pthread_mutex_t mutex_a = PTHREAD_MUTEX_INITIALIZER;

/*
 *NOTE:Different threads are handling differnt primes so it is enough to lock the mutex 
 *     only while updating the values associated with struct node array
 *
 *	i.e if a number X is divisible by say thread1 and it is also divisible by thread2 even then
 *	there cannot be a race condition as the factors by which it is divisible wiil be different
 *	set_of_primes(thread1) INTERSECTION set_of_primes(thread2) = NULL
 */


/*
 *initializes and allocates memory for the required members of the structure node......
 *creats the array and initializes q_x to (seed+x)^2-qfs_num
 *vector node->v to 00000
 *returns a pointer to the node array
 */

struct node **build_arr (mpz_t offset, int v_size)
{
	struct node **p,**t;
	register unsigned int i;
	mpz_t seed;
	
	//allocation for the array of pointers of main seive array
	p = malloc (sizeof(struct node *) * (sieve_offset));
	
	//set seed
	mpz_init (seed);
	mpz_sqrt (seed, qfs_num);
	
	//Allocate and initialize the various elements of struct
	for (i=0 ;i<sieve_offset ;i++){
		t=p+i;
		//allocate memory for structure
		(*t) = malloc (sizeof(struct node));
		mpz_init ((*t)->q_x);
		mpz_init ((*t)->m);
		
		//set m
		mpz_set ((*t)->m, offset);
		mpz_add_ui ((*t)->m, (*t)->m, i);
		
		//calculate q_x
		mpz_add ((*t)->q_x, seed, (*t)->m);
		mpz_mul ((*t)->q_x, (*t)->q_x, (*t)->q_x);
		mpz_sub ((*t)->q_x, (*t)->q_x, qfs_num);
		
		//allocate for the vector v and set it to zero
		(*t)->v = calloc(v_size, sizeof(int));
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
	for (i=0 ; i<nthread ;i++) 
		pthread_join (thread_id[i], NULL);
	
	//free the thread data structures
	for (i=0 ;i<nthread ;i++)
		free (*(thread_args+i));
	
	free (thread_id);
	free (thread_args);
	//done
}


/*
 *threaded function control will be passsed to this function on thread creation
 *NOTE:Each threaded function handles only a unique set of primes
 */

static void *threaded_factorize (void *param)
{
	int i=0,j=0;
	int level,step;
	int t1;
	const struct thread_data *td = (struct thread_data *)param;
	
	//current word accessing in v
	int idx;
	//masks the current bit accessed in a word
	int mask;

	
	struct node **temp;
	
	idx = td->thread_idx;
	mask = 1;
	i = td->thread_idx * WORD_SIZE;
	
	
	for ( ; i < td->B ; i++){
		temp = td->start;

		//On wordsize boundaries execute the loop
		if ((i % WORD_SIZE == 0) && (i != (td->thread_idx * WORD_SIZE) ) ) {
			//go to next word or element in v[] which should be handled by this thread
			idx = idx + nthread;
			
			mask=1;		//reinitialize mask to point to bit 0
			
			//skip the primes which should not be handled by this thread
			i = i + (nthread-1) * WORD_SIZE;
			
			//this thread  has completed its job
			if (i >= td->B)
				return NULL;
		}
		
		
		if (!i) {
			/*if a number is negative make it positive and update vector*/
			t1 = sieve_offset;
			while (t1--) {
				if (mpz_sgn ((*temp)->q_x) < 0) {
					///mutex LOCK
					pthread_mutex_lock (& mutex_a);
					
					mpz_neg ((*temp)->q_x, (*temp)->q_x);		//update q_x
					(*temp)->v[idx] = (*temp)->v[idx] ^ mask;	//update v
					
					///UNLOCK mutex
					pthread_mutex_unlock (& mutex_a);
				}
				temp++;
			}
		}
		else {
			level=0;
			for (j=0; j<sieve_offset ;j++){
				//if number is not divisable continue
				if (!mpz_divisible_ui_p((*(temp+j))->q_x, *(td->prime+i)))
					continue;
				//the number is divisble
				else {
					/*level is required to determine the number of steps or
					 *elements to be skipped for the erathosthenes sieve
					 */
					level=0;
					
					while ((!level) || mpz_divisible_ui_p ((*(temp+j))->q_x, *(td->prime+i)) ) {
						level++;
						step = pow (*(td->prime+i),level);
						
						if (j < *(td->prime+i)){
							//Erathosthenes sieve
							repete (temp ,j ,step ,*(td->prime+i)
									,idx, mask);
						}
						else {
							//trail division
							trail (temp ,j ,*(td->prime+i), idx, mask);
							break;
						}
						
					}
				}
			}
		}
		mask = mask << 1;		//update mask for next prime
	}
	return NULL;
}



/*
 *Code to do the trail division
 *ptr ---pointing to element corresponding to -M
 *offset --from above element
 *pr --present prime we r dealing
 */

static inline void trail(struct node **ptr,const int offset,const int pr,
			  const int idx,const int mask)
{
	while( mpz_divisible_ui_p ((*(ptr+offset))->q_x, pr)) {
		
		///LOCK the mutex as we changing the value realted to struct node array
		pthread_mutex_lock (& mutex_a);
		
		mpz_tdiv_q_ui ((*(ptr+offset))->q_x,  (*(ptr+offset))->q_x  ,pr);		//update q_x
		(*(ptr+offset))->v[idx] = (*(ptr+offset))->v[idx] ^ mask;	//update v
		
		///UNLOCK the mutex
		pthread_mutex_unlock (& mutex_a);
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

static inline void repete (struct node **ptr,const int offset,const int step,const int pr,
			   const int idx,const int mask)
{
	struct node **limit;
	limit = ptr+sieve_offset;
	
	///LOCK the mutex as struct node array will be modified
	pthread_mutex_lock (&mutex_a);
	
	for (ptr=ptr+offset ; ptr < limit  ; ptr=ptr+step) {
		mpz_tdiv_q_ui ((*ptr)->q_x, (*ptr)->q_x, pr);			//update q_x
		(*ptr)->v[idx] = (*ptr)->v[idx] ^ mask;				//update v
	}
	///UNLOCK the mutex
	pthread_mutex_unlock (&mutex_a);
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

