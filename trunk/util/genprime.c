#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>

#define OUTPUT_BASE 10

int main(int argc, char *argv[])
{
	mpz_t seed,z;
	FILE *fp;
	int i,count=0;


	if (argc < 3) {
		fprintf (stderr,"Error\nUsage:%s <No.of.primes> <filename>\n",argv[0]);
		exit(-1);
	}

	count = atoi (argv[1]);
	fprintf(stderr,"Generating %d Primes.....",count);
	mpz_init2 (seed,-1);
	mpz_init (z);

	if(!(fp = fopen (argv[2],"w"))){
		printf ("ERROR:could not open the file");
		return -1;
	}

	fprintf (fp,"-1\n");
	for(i=0;i<count-1;i++) {
		mpz_nextprime (z,seed);
		mpz_out_str (fp,OUTPUT_BASE,z);
		fputc ('\n',fp);
		mpz_set (seed,z);
	}
	fprintf (stderr,"done\n");
	fclose(fp);
	return 0;
}
