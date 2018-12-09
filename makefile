GCC_FLAGS = -Wall

SOURCE = interpol_gsl.c
DESTINATION = interpol_gsl
LIBRARY_FLAGS = -lgsl -lm -lgslcblas

all:	$(SOURCE)
	gcc $(GCC_FLAGS) $(SOURCE)  -o $(DESTINATION) $(LIBRARY_FLAGS) 

clean:
	rm main
