/*
 *Generates composite numbers with factors of almost equal sizes
 *n = p * q;
 *n=>composite number
 *p and q  => primes
 */
 
#include<stdio.h>
#include<gmp.h>
#include<time.h>
#include<stdlib.h>

//Primality test Iteration 5-10 is pretty good
#define MAX_ITER 10

#define PRIME 1

/*to pick a public key among the foll.*/
unsigned long primes[]={2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,51,53,59,61,71,73,79,83,89,93,97};

struct RSA{
	mpz_t publickey;
	mpz_t privatekey;
	mpz_t n;
	mpz_t phi_n;
};


/*public:{publickey,n}
 *private:{privatekey,n}
 */

gmp_randstate_t STATE;


void generate_rsa_key(struct RSA *);
unsigned long seed_time(void);
inline void encrypt(mpz_t *,mpz_t *,struct RSA *);
inline void decrypt(mpz_t *,mpz_t *,struct RSA *);
inline void print_rsa_keys(struct RSA );
void store_rsa(char *,struct RSA *);
void load_rsa(char *,struct RSA *);

int keysize;

int main(int argc,char *argv[])
{
	struct RSA rsa;
	//initialize Random intgen generation
	unsigned long a;
	
	if (argc != 2){
		fprintf (stderr, "USAGE: %s <no of bits(key)>", argv[0]);
		exit (1);
	}
	keysize = atoi (argv[1]);
	
	a=seed_time();
	gmp_randinit_default(STATE);
	gmp_randseed_ui(STATE,a);
	
	//load_rsa("tmp",&rsa);
	generate_rsa_key(&rsa);
	
	printf("\n\n TESTING ENCRYPTION");
	
	store_rsa("tmp",&rsa);
	print_rsa_keys(rsa);
	//encryption
	mpz_t msg,encryp,decryp;
	mpz_init(msg);
	mpz_init(encryp);
	mpz_init(decryp);
	mpz_set_ui(msg,0xDDAABEFF);
	encrypt(&encryp,&msg,&rsa);
	gmp_printf("\n\nENCRYPTION:\n\nmsg  : %Zx\nencryp: %Zx",msg,encryp);
	decrypt(&decryp,&encryp,&rsa);
	gmp_printf("\n\nDECRYPTION:\n\ndecryp: %Zx\n",decryp);
	mpz_clear(msg);
	mpz_clear(encryp);
	mpz_clear(decryp);
	return 0;
}


/*
 *This function generates the rsa keys 
 *Both public and private keys are generated ......
 */

void generate_rsa_key(struct RSA *rsa)
{
	mpz_t p,q;                    //Primes
	
	mpz_t p_tmp,q_tmp,temp;       //temporary values for calculations
	int i,publickey;

	mpz_init(rsa->n);
	mpz_init(rsa->phi_n);
	mpz_init(rsa->publickey);
	mpz_init(rsa->privatekey);
	
	mpz_init2(p,keysize/2);
	mpz_init2(q,keysize/2);

	mpz_init(p_tmp);
	mpz_init(q_tmp);
	mpz_init(temp);
	
	//Generate primes
	do{
		mpz_urandomb(p,STATE,keysize/2);
	}while(mpz_probab_prime_p(p,MAX_ITER)<PRIME);
	do{
		mpz_urandomb(q,STATE,keysize/2);
	}while(mpz_probab_prime_p(q,MAX_ITER)<PRIME);
	
	mpz_mul(rsa->n,p,q);

	//claculate phi_n
	mpz_sub_ui(q_tmp,q,1);
	mpz_sub_ui(p_tmp,p,1);
	mpz_mul(rsa->phi_n,p_tmp,q_tmp);
	//phi_n done

	//Public key generation
	for(i=0 ; i < sizeof(primes)/sizeof(int) ; i++){
		if(mpz_gcd_ui(temp,rsa->phi_n,primes[i])==1){
			mpz_set_si(rsa->publickey,primes[i]);
			publickey=primes[i];
			break;
		}
	}
	//Public Key generation is Done

	//private Key
	for(i=1;i<1000;i++){
		mpz_mul_ui(p_tmp,rsa->phi_n,i);
		mpz_add_ui(temp,p_tmp,1);
		if(mpz_tdiv_r_ui(p_tmp,temp,publickey)==0){
			mpz_tdiv_q(rsa->privatekey,temp,rsa->publickey);
			break;
		}
	}

	gmp_printf("\n\n%Zd %Zd\n\n",p,q);
	//CLEAR ALL
	mpz_clear(temp);
	mpz_clear(p);
	mpz_clear(q);
	mpz_clear(p_tmp);
	mpz_clear(q_tmp);
}


/*
 *This function is used to set the seed for the generation of random numbers
 *These random numbers are used to generate the primes p and q
 *The seed is set using the time functions in time.h.
 */


unsigned long seed_time()
{
	static int is_this_first_time=1;
	int num;
	if(is_this_first_time){
		int a,b,c,d,e;
		char str1[4],str2[4];
		time_t now;
		now = time((time_t *)NULL);
		sscanf(ctime(&now),"%s %s %d %d:%d:%d %d",str1,str2,&a,&b,&c,&d,&e);
		num = num + (num *a) + (b+c)*num+(d+e)*num+a+b+c+d+e;
		is_this_first_time=0;
	}
	return num;
}


//Function to encrypt the message
//(plaintext^publickey )mod n
inline void encrypt(mpz_t *ciphertext,mpz_t *plaintext,struct RSA *rsa)
{
	mpz_powm(*ciphertext,*plaintext,rsa->publickey,rsa->n);
}


//Function to decrypt the message
//(ciphertext^privatekey)mod n
inline void decrypt(mpz_t *plaintext,mpz_t *ciphertext,struct RSA *rsa)
{
	mpz_powm(*plaintext,*ciphertext,rsa->privatekey,rsa->n);
}


//Function to print the keys
inline void print_rsa_keys(struct RSA rsa)
{
	gmp_printf("Publickey : %Zx\nPrivatekey: %Zx\nn\t  : %Zx\nphi(n)    : %Zx"
		   ,rsa.publickey,rsa.privatekey,rsa.n,rsa.phi_n);
}

/*
 *Function that stores the rsa keys to a specified file
 */
void store_rsa(char *filename,struct RSA *rsa)
{
	FILE *fp;
	if(!(fp=fopen(filename,"w"))){
		printf("ERROR:Could not creat file %s",filename);
		exit(1);
	}
	
	mpz_out_str(fp,16,rsa->publickey);
	fputc('\n',fp);
	mpz_out_str(fp,16,rsa->privatekey);
	fputc('\n',fp);
	mpz_out_str(fp,16,rsa->n);
	fputc('\n',fp);
	mpz_out_str(fp,16,rsa->phi_n);
	fputc('\n',fp);
	fclose(fp);
}


/*
 *Function that loads the rsa keys from the specified file 
 *the loaded keys can be used to do the encryption and Decryption
 */

void load_rsa(char *filename,struct RSA *rsa)
{
	FILE *fp;
	if(!(fp=fopen(filename,"r"))){
		printf("ERROR:Could not open the File");
		exit(1);
	}
	
	mpz_init(rsa->n);
	mpz_init(rsa->phi_n);
	mpz_init(rsa->publickey);
	mpz_init(rsa->privatekey);
	
	mpz_inp_str(rsa->publickey,fp,16);
	fgetc(fp);
	mpz_inp_str(rsa->privatekey,fp,16);
	fgetc(fp);
	mpz_inp_str(rsa->n,fp,16);
	fgetc(fp);
	mpz_inp_str(rsa->phi_n,fp,16);
	fgetc(fp);
	fclose(fp);
}
