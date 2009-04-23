/*
 *The functions in this file takes input from user and sets up some important global
 *varibles used by other programs
 *functions : 
	main, sieve, group
 */

#include<stdio.h>
#include<stdlib.h>
#include<gmp.h>
#include<math.h>
#include<string.h>
#include<getopt.h>

/*For some common definations*/
#include "qfs.h"

static void print_usage (int exit_status);
static void sieve();

char *program_name;


int mem_limit = 8;
int prime_multiple = 2;
int inp_base = 10;
unsigned long int sieve_offset;
int smooth_extra = 32;
int nthread = 1;

int main (int argc, char *argv[])
{
	int next_option;
	char *inp_filename=NULL,  *inp_num=NULL;
	
	const struct option long_options[] = {
		{ "help",		0,NULL,'h'},
		{ "mem-limit",		1,NULL,'m'},
		{ "prime-multiple",	1,NULL,'p'},
		{ "infile",		1,NULL,'f'},
		{ "number",		1,NULL,'n'},
		{ "verbose",		0,NULL,'v'},
		{ "base",		1,NULL,'b'},
		{ "smooth-extra",	1,NULL,'s'},
  		{ "nthread",		1,NULL,'t'},
		{ NULL,			0,NULL,0}	//required for the end of array
	};
	
	const char short_options[]="hvm:p:f:n:b:s:t:";
	
	program_name=argv[0];	//store the program name
	
	do {
		next_option = getopt_long (argc, argv, short_options, long_options, NULL);
	
		switch (next_option)
		{
			case 'h':			//print help
				print_usage(0);
				break;
				
			case 'v':			//verbose mode
				verbose = 1;
				break;
				
			case 'm':			//set memory limit
				mem_limit = atoi(optarg);
				break;
				
			case 'p':			//set prime multiple
				prime_multiple = atoi(optarg);
				break;
				
			case 'f':			//specify input file name
				inp_filename = optarg;
				break;
				
			case 'n':			//specify input number
				inp_num = optarg;
				break;
				
			case 's':			//specify number of additional smooth nums
				smooth_extra = atoi(optarg);
				break;
				
			case 'b':			//base of input number
				inp_base = atoi(optarg);
				break;
				
			case 't':			//number of threads
				nthread = atoi (optarg);
				break;
				
			case -1:
				break;
				
			case '?':			//unknown option
				print_usage(1);
				break;
				
			default:			//something has gone wrong!!!!!!!
				abort();
				
		}
	} while (next_option != -1);
	
	
	if (inp_filename==NULL && inp_num==NULL){	//no input specified
		print_usage(1);
	} 
	else if (inp_filename!=NULL && inp_num!=NULL){	//more than one inputs specified
		print_usage(1);
	}
	
	/*initialize qfs_num and give control to sieve*/
	mpz_init(qfs_num);
	
	if (inp_filename){	//read from file
		FILE *fp;
		if ((fp = fopen(inp_filename,"r")) == NULL) {
			fprintf (stderr, "ERROR:Could not access the file %s\n", inp_filename);
			exit (1);
		}
		
		mpz_inp_str( qfs_num, fp, inp_base);
		fclose (fp);
	}
	else {
		mpz_set_str( qfs_num, inp_num, inp_base);
	}
	
	/*
	 *NOTE:any number of form x^y cannot be factorised by this sieve
	 *TODO:Think of some way to check a number of form x^y efficiently
	 *here we Check only whether the number is a perfect square 
	 *if yes output the result and exit
	 *Necessary because Q.S cant factorize the x^y kind of numbers
	 */
	if (mpz_perfect_square_p (qfs_num) ){
		mpz_t a;
		mpz_init(a);
		
		mpz_sqrt(a,qfs_num);
		write_result(a);
		
		mpz_clear(a);
		exit(0);
	}
	
	//so the number is not of form x^y continue
	sieve();
	
	return 0;
}


///Declarions that r global to the project and used almost by all functions

mpz_t qfs_num;

int verbose=0;

///done



/*
 *Function to do the sieving based on the input given
 */
static void sieve()
{
	float log_qfs_num;
	int smooth_cnt;
	int *prime,prime_cnt,B;
	
	struct node **start_ptr,**smooth_ptr;
	mpz_t offset1,offset2;	
	
	/*Building Factor Base*/
	//here we have implemented ln as (no of bits+1)/ln(2)
	
	log_qfs_num = mpz_sizeinbase(qfs_num,2)/(1.442695);
	
	B = pow( exp( sqrt( log_qfs_num*log(log_qfs_num) ) ),0.353553)+1;
		
	prime_cnt = build_factor_base(&prime,B,prime_multiple);
	//done
	printf("logn=%f\n",log_qfs_num);
	//setup sieve offset
	int mem_per_node;
	
	mem_per_node = sizeof(struct node);	//size of node
	mem_per_node += (prime_cnt/WORD_SIZE+1) * sizeof(unsigned int);		//allocated to node->v
	mem_per_node += ((mpz_sizeinbase(qfs_num,2)+1)*2)/sizeof(mp_limb_t) + 1;//allocated to node->q_x
	mem_per_node += (2*sizeof(mp_limb_t));					//allocated to node->m
	
	sieve_offset = (mem_limit * 1024 * 1024)/(mem_per_node + sizeof(struct node *));
	
	if(sieve_offset > (unsigned int)pow(B,3)){
		sieve_offset = pow(B,3);
	}
	//done
	
	if (verbose){	
		gmp_printf ("\nPrime multiple %d \nPrimes Bound: %d  \nNo of Primes:%d\n"
				"Smooth numbers required:%d\n"
				"Sieve offset calculated %u for memory_limit of %d MB \n"
				"Number of threads Selected %d"
				,prime_multiple,B,prime_cnt, prime_cnt+32, sieve_offset, mem_limit, nthread);
		
		gmp_printf ("\nData Collection Phase......\n");
	}
	
	//Collecting smooth numbers
	smooth_cnt=0;
	smooth_ptr=malloc(sizeof(struct node *) * (prime_cnt + smooth_extra + smooth_extra));
	
	mpz_init(offset1);
	mpz_init(offset2);
	
	
	/*
	*The loop below calculates the required number of smooth numbers and stores it to smooth_ptr**
	*section 1 goes along the positive direction - offset1 keeps track of it
	*section 2 goes along the negative direction - offset2 keeps track of it
	*/
	while(1){
		///section 1
		if(verbose){
			gmp_printf ("smooth_numbers_collected = %4d(%.2f complete)\t",
				    smooth_cnt,(float)smooth_cnt*100/(prime_cnt+32));
			gmp_printf("offset1:%Zd\n", offset1);
		}
		
		start_ptr = build_arr (offset1, prime_cnt/WORD_SIZE+1);
		
		factorize (start_ptr, prime, prime_cnt);
		
		smooth_cnt = group (start_ptr, smooth_ptr, smooth_cnt, prime_cnt);

		//we have got enough smooth numbers
		if (smooth_cnt > prime_cnt + smooth_extra)
			break;
		
		mpz_add_ui (offset1, offset1, sieve_offset);
		

		///section 2
		mpz_sub_ui(offset2,offset2,sieve_offset);
		
		if(verbose){
			gmp_printf ("smooth_numbers_collected = %4d(%.2f complete)\t",
				    smooth_cnt, (float)smooth_cnt*100/(prime_cnt+32));
			gmp_printf ("offset2:%Zd\n", offset2);
		}
		
		start_ptr = build_arr(offset2,prime_cnt/WORD_SIZE+1);
		
		factorize (start_ptr,prime,prime_cnt);
		
		smooth_cnt  = group (start_ptr, smooth_ptr, smooth_cnt, prime_cnt);
		
		//we have got enough smooth numbers
		if (smooth_cnt > prime_cnt + smooth_extra )
			break;
		
	}
	
	mpz_clear(offset1);
	mpz_clear(offset2);
	
	if(verbose){
		gmp_printf ("\nData collection completed ...%d smooth numbers collected",smooth_cnt); 
		gmp_printf ("\nData Analysis phase.....\n");
	}
	
	//Smooth_ptr will be de-allocated in gaussian.c
	gaussian(smooth_ptr  ,prime_cnt  ,smooth_cnt);
	
	free(prime);
	
	mpz_clear(qfs_num);
}



/*
 *Function that prints the details about the options
 */
static void print_usage(int exit_status)
{
	fprintf (stderr,
		"\nUsage: %s [OPTIONS] -n NUMBER.\n"
		"   or: %s [OPTIONS] -f INPUT_FILE_NAME.\n"
		"Factorizes the number specified in the file, or otherwise. using basic quadratic sieve.\n"
		"\nMandatory arguments to long options are mandatory for short options too.\n"
		" -b, --base=BASE		specify the input base[BASE] default 10\n"
		" -f, --infile=FILENAME		specify the input file containing number\n"
		" -h, --help			display this help and exit\n"
		" -m, --mem-limit=MEMSIZE	specify the memory limit(in MB) for program[refer README file]\n"
		" -n, --number=NUMBER_STRING	specify the number to be factorised\n"
		" -p, --prime-multiple=PRIME	specify the number of primes to use[refer README file]\n"
		" -s, --smooth-extra=SMT        specify the additional smooth numbers to calculate.\n" 
		" -t, --nthread=NTHREAD         specify the number of threads NTHREAD.\n"
		" -v, --verbose			prints all the info as program executes\n\n" 
		,program_name, program_name);
	
	exit (exit_status);
}




/*
 *collects all the numbers with q_x < limit to the top of the table
 *and arranges them wrt magnitude of m
 *returns the no of smooth numbers
 */

int group (struct  node **start_ptr,struct node **smooth_ptr,int smooth_cnt ,int prime_cnt)
{
	int i,v_size;
	struct node *tmp;
	//vector  s size
	v_size = (prime_cnt/WORD_SIZE + 1) * sizeof(unsigned int);
	
	for (i=0 ; i<sieve_offset ; i++) {
		
		if ( mpz_cmp_ui((*(start_ptr+i))->q_x,1) == 0 ) {

			smooth_ptr[smooth_cnt++] = start_ptr[i];
			
			//we have got enough smooth nunbers
			if (smooth_cnt > prime_cnt + smooth_extra + smooth_extra)
				break;
		}
		else {
			/*clear all the memory allocated to the structure as its not a smooth number*/
			tmp = *(start_ptr+i);
			
			mpz_clear(tmp->q_x);
			mpz_clear(tmp->m);
			free(tmp->v);
			free(tmp);
		}
	}
	
 	/*we have got enough smooth numbers but clean up the mess 
	 *i.e release the rest of memory allocated
	 */
	if(i < sieve_offset){
		for(i=i+1 ; i<sieve_offset ; i++){
			tmp = *(start_ptr+i);
			
			mpz_clear(tmp->q_x);
			mpz_clear(tmp->m);
			
			free(tmp->v);
			free(tmp);
		}
	}
	
	free (start_ptr);
	
	//return the number of smooth numbers found
	return smooth_cnt;
}
