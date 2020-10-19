GCC = gcc
GCC_FLAGS = -Wall -O2 -lm
GCC_OBJ_FLAGS = -Wall -O2 -c

main: main.o bin-packing.o population.o chromosome.o
	$(GCC) $(GCC_FLAGS) main.o bin-packing.o population.o chromosome.o \
		-o main.out

bin-pack-test: bin-pack-test.o bin-packing.o population.o chromosome.o
	$(GCC) $(GCC_FLAGS) bin-pack-test.o bin-packing.o population.o \
		chromosome.o -o bin-pack-test.out

pop-test: pop-test.o population.o chromosome.o
	$(GCC) $(GCC_FLAGS) pop-test.o population.o chromosome.o \
		-o pop-test.out

chrom-test: chrom-test.o chromosome.o
	$(GCC) $(GCC_FLAGS) chrom-test.o chromosome.o \
		-o chrom-test.out

clean:
	rm main.o bin-packing.o population.o chromosome.o bin-pack-test.o \
		pop-test.o chrom-test.o main.out bin-pack-test.out \
		pop-test.out chrom-test.out

main.o: main.c
	$(GCC) $(GCC_OBJ_FLAGS) main.c

bin-packing.o: bin-packing.c
	$(GCC) $(GCC_OBJ_FLAGS) bin-packing.c

population.o: population.c
	$(GCC) $(GCC_OBJ_FLAGS) population.c

chromosome.o: chromosome.c
	$(GCC) $(GCC_OBJ_FLAGS) chromosome.c

bin-pack-test.o: bin-pack-test.c
	$(GCC) $(GCC_OBJ_FLAGS) bin-pack-test.c

chrom-test.o: chrom-test.c
	$(GCC) $(GCC_OBJ_FLAGS) chrom-test.c

pop-test.o: pop-test.c
	$(GCC) $(GCC_OBJ_FLAGS) pop-test.c

