##Make file for the qfs
#FIXME:  Very Badly written Makefile ....

#Source Directories
SRC_DIR=./src

#Utilites Directories (Other than main source)
UTIL_DIR=./util

basicqs : main.o factorize.o factor_base.o  gaussian.o
	gcc $(CFLAGS) -o basicqs $(SRC_DIR)/main.o $(SRC_DIR)/factorize.o $(SRC_DIR)/factor_base.o $(SRC_DIR)/gaussian.o -lm -lgmp -lpthread

main.o : $(SRC_DIR)/main.c $(SRC_DIR)/qfs.h
	gcc -c $(CFLAGS) -Wall $(SRC_DIR)/main.c

factorize.o : $(SRC_DIR)/factorize.c $(SRC_DIR)/qfs.h
	gcc -c $(CFLAGS) -Wall $(SRC_DIR)/factorize.c

factor_base.o : $(SRC_DIR)/factor_base.c $(SRC_DIR)/qfs.h
	gcc -c $(CFLAGS) -Wall $(SRC_DIR)/factor_base.c

gaussian.o : $(SRC_DIR)/gaussian.c $(SRC_DIR)/qfs.h
	gcc -c $(CFLAGS) -Wall $(SRC_DIR)/gaussian.c


###for generatin composite numbers for testing
rsa : $(UTIL_DIR)/rsa_main.c
	gcc $(CFLAGS) -Wall $(UTIL_DIR)/rsa_main.c -o rsa -lgmp

###for generating primes.txt file
genprime: $(UTIL_DIR)/genprime.c
	gcc $(CFLAGS) -Wall $(UTIL_DIR)/genprime.c -o genprime -lgmp

##clean usual
.PHONY : clean

clean :
	rm -f basicqs rsa genprime $(SRC_DIR)/*.o $(UTIL_DIR)/*.o *.o
