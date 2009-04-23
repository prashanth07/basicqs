/*
 *Functions in this file deals with gaussian elimination
 *finding the factors of qfs_num
 *functions
 	write_resule,  gaussian
 */

#include<stdio.h>
#include<stdlib.h>
#include<gmp.h>
#include<math.h>
#include<string.h>

/*For some common definations*/
#include "qfs.h"

extern int inp_base;

/*
 *v is the copy of the vector
 *pr is the no the array with no of bits greater than or equal to smooth_cnt
 */

struct gauss{
	int *v;
	int *pr;
};

static struct gauss ** copy_gauss (struct node **source,int smooth_cnt,int prime_cnt);

static void find_factors (struct gauss **ptr,  struct node **start, int prime_cnt  ,int smooth_cnt);

static void free_gauss (struct gauss **target, int smooth_cnt);

static void clean_up (struct node **smooth_ptr, unsigned int smooth_cnt);

void print_gauss(char *,struct gauss **ptr , int prime_cnt ,int smooth_cnt);
/*
 *Does the gaussian Elimination on the gauss structure
 */
 
void gaussian (struct node **start  ,int prime_cnt  ,int smooth_cnt)
{
	int mask;		//mask for the bit accessed
	int word;		//current word being accessed
	int i,j,k;
	struct gauss **ptr,*temp;
	int found;
	
	ptr=copy_gauss(start,smooth_cnt,prime_cnt);
	
	mask=1,word=0;
	
	for (i=0 ;i<prime_cnt ;i++) {
		
		if (i % WORD_SIZE == 0 && i!=0){	//update word and mask
			mask=1;
			word++;
		}
		
		/*get a vectror with corresponding bit = 1
		 *bit is identified by a combination of word and bit
		 */
		
		found=0;
		for (j=i; j<smooth_cnt ;j++){		
			if ( ((*(ptr+j))->v[word] & mask) ){
				if (i == j){
 					found=1;	//vector found
					break;
				}
				else{
					temp = *(ptr+i);
					*(ptr+i) = *(ptr+j);
					*(ptr+j) = temp;
					found=1;
				}
			}
		}
		
		if(!found){
			mask = mask<<1;
			continue;
		}
		
		/*
		 *Do Gaussian elimination using the above found vector
		 *set that corresponding bit to zero in the remaining vectors
		 *using the vectors found above
		 */
		
		for (j=i+1 ; j<smooth_cnt  ;j++){
			if(  (*(ptr+j))->v[word] & mask  ){
				//the corresponding bit in the vector = 1
				for (k=0 ; k < (prime_cnt/WORD_SIZE)+1 ;k++) {
					(*(ptr+j))->v[k] = (*(ptr+j))->v[k] ^ (*(ptr+i))->v[k]; 
				}
				for (k=0 ; k < (smooth_cnt/WORD_SIZE)+1 ;k++) {
					(*(ptr+j))->pr[k] = (*(ptr+j))->pr[k] ^ (*(ptr+i))->pr[k];
				}
			}
		}
		mask = mask << 1;
	}
	
	//find the factors
	find_factors(ptr,  start , prime_cnt  , smooth_cnt);
	
	//free the memory allocated
	free_gauss(ptr,smooth_cnt);
}



/*
 *Finds the factors after gaussian elimination is done
 *if the code below is confusing check the output of
 *print gauss before and after the gaussian elimination
 */

static void find_factors(struct gauss **ptr,  struct node **start,
			 int prime_cnt, int smooth_cnt)
{
	int i,j,mask,word;
	mpz_t a,b,temp,seed;
	
	mpz_init (a);
	mpz_init (b);
	mpz_init (temp);
	mpz_init (seed);
	
	//set seed
	mpz_sqrt (seed, qfs_num);
	
	for (i = prime_cnt  ; i < smooth_cnt ;i++) {
		mask=1;
		word=0;
		
		mpz_set_si (a,1);
		mpz_set_si (b,1);
		mpz_set_si (temp,1);
		
		for (j=0 ; j< smooth_cnt ;j++) {
			
			if(j % WORD_SIZE == 0 && j!=0) {
				word++;
				mask=1;
				if( ! (((*(ptr+i))->pr[word]) & (~0))  ){
					j = j+WORD_SIZE-1;
					continue;
				}
			}
			
			if ( (*(ptr+i))->pr[word] & mask){
				//j th smooth number in **start is involved 
				//a=a * ((seed+(*(start+i))->m)^2-qfs_num)
				mpz_set (temp, (*(start+j))->m);
				mpz_add (temp, seed, temp);
				mpz_mul (b, b, temp);
				
				mpz_mul (temp, temp, temp);
				mpz_sub (temp, temp, qfs_num);
				mpz_mul (a, a, temp);
			}
			
			mask = mask << 1;
		}
		
		mpz_sqrt (a, a);
		mpz_sub (temp, a, b);
		
		//factor = gcd_EA(a-b,qfs_num);
		mpz_gcd (temp, temp, qfs_num);

		if ( (mpz_cmpabs_ui (temp, 1)) && (mpz_cmpabs (temp, qfs_num)) ){
			//we have found a useful factor i.e != 1 or !=qfs_num
			mpz_tdiv_q (a, qfs_num, temp);
			
			write_result (temp);				//write result to a file
			
			mpz_clear (seed);
			mpz_clear (temp);
			mpz_clear (a);
			mpz_clear (b);
			
			//release all memory
			clean_up (start, smooth_cnt);
			return;
		}
	}
}



/*
 *Releases memory allocated to smooth_ptr
 */

static void clean_up (struct node **smooth_ptr, unsigned int smooth_cnt)
{
	int i;
	for(i=0;i<smooth_cnt;i++){
		mpz_clear ((*(smooth_ptr+i))->q_x);
		mpz_clear ((*(smooth_ptr+i))->m);
		free ((*(smooth_ptr+i))->v);
		free (*(smooth_ptr+i));
	}
	free (smooth_ptr);
}



/*
 *Copies the vectors from node to a new local structure(gauss structure)
 */

static struct gauss ** copy_gauss(struct node **source,int smooth_cnt,int prime_cnt)
{
	int i,j;
	int mask=1,word=0;
	struct gauss **target;
	
	target = malloc(sizeof(struct gauss *) * smooth_cnt);
	
	for(i=0;  i<smooth_cnt ;i++){
		(*(target+i)) = malloc(sizeof(struct gauss)*smooth_cnt);
		(*(target+i))->v = malloc (sizeof(int)*(prime_cnt/WORD_SIZE+1));
		(*(target+i))->pr = malloc (sizeof(int)*(smooth_cnt/WORD_SIZE+1));
		
		memset ((*(target+i))->pr, 0 ,sizeof(int)*(smooth_cnt/WORD_SIZE+1));
		
		for (j=0 ;j<(prime_cnt/WORD_SIZE + 1) ;j++)
			(*(target+i))->v[j] = (*(source+i))->v[j];
	}
	
	mask=1;
	word=0;
	
	for (i=0 ; i<smooth_cnt ; i++){
		if( i % WORD_SIZE == 0 && i!=0){
			mask=1;
			word++;
		}
		(*(target+i))->pr[word] = (*(target+i))->pr[word] ^ mask;
		mask = mask << 1;
	}
	
	return target;
}


#ifdef POIU_DEBUG
/*
 *Writes the Gauss structure to a file specified
 */

void print_gauss (char *s,struct gauss **ptr , int prime_cnt ,int smooth_cnt)
{
	int i,j;
	FILE *fp;
	
	if ((fp = fopen(s,"w")) == NULL){
		fprintf(stderr,"\n\nERROR:(print_gauss)Cannot create file %s\n\n", s);
		exit(1);
	}
	
	for (i=0 ; i<smooth_cnt ; i++){
		fprintf (fp, "\n");
		
		for (j=0 ; j<(prime_cnt/WORD_SIZE)+1 ;j++)
			fprintf(fp," %8X ",(*(ptr+i))->v[j]);
		fprintf(fp,"\t");
		
		for(j=0 ; j<(smooth_cnt/WORD_SIZE)+1 ;j++)
			fprintf(fp," %8X ",(*(ptr+i))->pr[j]);
	}
	
	fclose(fp);
}

#endif


/*
 *Writes the result to a file
 */

void write_result(mpz_t a)
{
	mpz_t b;
	FILE *fp;
	
	if(!(fp=fopen("output.txt","w"))){
		fprintf( stderr, "ERROR:Cannot write to the file output.txt");
		exit(1);
	}
	mpz_init(b);
	
	mpz_tdiv_q(b,qfs_num,a);
	gmp_fprintf (fp,"Factors of ",qfs_num);
	mpz_out_str(fp,inp_base,qfs_num);
	gmp_fprintf (fp," are\n");
	mpz_out_str(fp,inp_base,a);
	gmp_fprintf (fp," and ");
	mpz_out_str(fp,inp_base,b);
	fclose(fp);
	
	gmp_fprintf (stdout,"\nFactors of ",qfs_num);
	mpz_out_str(stdout,inp_base,qfs_num);
	gmp_fprintf (stdout," are\n");
	mpz_out_str(stdout,inp_base,a);
	gmp_fprintf (stdout," and ");
	mpz_out_str(stdout,inp_base,b);
	printf("\n");
	mpz_clear(b);
}



/*
 *Frees the gauss structure allocated during the call to gauss function
 */
 
static void free_gauss(struct gauss **target,int smooth_cnt)
{
	int i;
	register struct gauss **tmp=target;
	
	for(i=0;  i<smooth_cnt ;i++,tmp++){
		free( (*tmp)->v );
		free( (*tmp)->pr);
		free( *tmp );
	}
	free(target);
}

