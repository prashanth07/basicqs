/*******************************************************************************
		 Header for quadratic sieve implmentation - poiu
*******************************************************************************/

#ifndef QFS_GMP
#define QFS_GMP

/**word size is size of unsigned int in the target system**/
#define WORD_SIZE (sizeof(unsigned int)*8)

/*
 *define POIU_DEBUG to compile functions specifically for debugging
 *it has nothing else to do except provideing some data for debugging 
 */

//#define POIU_DEBUG 

extern int verbose;
extern mpz_t qfs_num;

struct node{
	int *v;
	mpz_t q_x;
	mpz_t m;
};

///function prototypes
struct node **build_arr(mpz_t offset,int v_size);

void factorize (struct node **start ,int *prime ,int B);

void gaussian(struct node **start  ,int prime_cnt  ,int smooth_cnt);

int group(struct  node **start_ptr,struct node **smooth_ptr,int smooth_cnt ,int prime_cnt);

unsigned int build_factor_base(int **p,int prime_bound ,int prime_multiple);

void write_result(mpz_t a);

//functions useful in debugging
void print_table(struct node **ptr,int *prime,char *s,int num,int B);

//another function print_gauss defined in gaussian.c 
#endif
