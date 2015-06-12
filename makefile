run: main.c
	mpicc -O2 -o primegap  main.c -lm -lgmp -std=c99

