/*********************************************************************
primegap

DESCRIPTION:

Takes a range of numbers as two arguments
Calculates in parallel the largest gap between two prime numbers in that range
Prints Run info: num of processes, timings and number of partitions used to distribute range

Edward Guloien
0749435
guloiej

*********************
RUNNING INSTRUCTIONS:

If running small ranges, then becareful not to run on too many processes.

The program uses the number of proccessor create the number of partitions.
If 32 processes are used then 32x32 = 1024 partitions are made. Thus running
a range of 0 -> 1000 might cause an issue.

This is not a problem for larger ranges.

***********************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#include <gmp.h>
#include <math.h>

int main(int argc, char* argv[])
{
  	int         rank; 		      /* rank of process      */
  	int         p;         	      /* number of processes  */

  	MPI_Init(&argc, &argv);
  	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &p);

	double t1 = MPI_Wtime();

	long int input1 = atol(argv[1]);
	long int input2 = atol(argv[2]);
	long int length_a = input2 - input1;	

	//init mpz range	
	mpz_t in1, in2, len_a;
	mpz_init_set_ui(in1, input1);
	mpz_init_set_ui(in2, input2);
	mpz_init_set_ui(len_a, length_a);

	long int part_width = length_a / (p * p);
	long unsigned int max_gap; 
	long unsigned int final_max_gap; 
	double fastestTime;
	double slowestTime;	
	
	//init prime variables
	mpz_t range_start, range_end;
	mpz_t prime, nextprime, primegap, max_primegap;
	mpz_init(prime);
	mpz_init(primegap);
	mpz_init(max_primegap);
	mpz_init(nextprime);
	mpz_init(range_start);
	mpz_init(range_end);


	for (int i = 0; i < p*p; i++){
		if ( rank  == i % p){
			//sets range
			if (i == 0 && p > 1 ){
				mpz_set(range_start, in1 );
				mpz_set_ui(range_end, part_width + input1 );
			}
			else if (i < p*p-1){
				mpz_set_ui(range_start, part_width*i + input1);
				mpz_set_ui(range_end, part_width*(i+1) + input1);
			}	
			else {
				mpz_set_ui(range_start, part_width*i + input1);
				mpz_set(range_end, in2);
			}	

			//gmp_printf("P:%d\t[%Zd, %Zd]\t(rank+1+i)modp:%d\timodp:%d\n",rank, range_start, range_end, (rank+1+i)%p, i%p);

			//sets starting num to prime
			if( mpz_probab_prime_p(range_start,25) == 2 )
				mpz_init_set(prime, range_start);
			else if ( mpz_probab_prime_p(range_start,25) == 1 ){
				mpz_init_set(prime, range_start);
			}
			else
				mpz_nextprime(prime, range_start);	

			//gets next prime, calculates gap
			mpz_nextprime(nextprime, prime);
			while ( mpz_cmp(range_end, nextprime) >= 0){
				mpz_sub(primegap, nextprime, prime);
	
				if ( mpz_cmp(primegap,max_primegap) > 0){
					mpz_set(max_primegap, primegap);
				}
				mpz_set(prime, nextprime);
				mpz_nextprime(nextprime, prime);
			}	
			if ( rank < p-1 ){
				mpz_sub(primegap, nextprime, prime);
				if ( mpz_cmp(primegap,max_primegap) > 0){
					mpz_set(max_primegap, primegap);
				}
			}
		}
	}	

	//converts local max_gap
	max_gap = mpz_get_ui(max_primegap);	
	
	double t2 = MPI_Wtime();

	//sends local max_gap to process 0	
	MPI_Reduce(&max_gap, &final_max_gap, 1, MPI_UNSIGNED_LONG, MPI_MAX, 0, MPI_COMM_WORLD);
               
	double t3 = t2 - t1;

	//prints timer each processes timer results for [0,1000000000]
	if ( p == 32 && input1 == 0 && input2 == 1000000000)
		printf("%d: timer: %f\n", rank, t3 );  	

	double t4;
	if ( rank == 0){
		t4 = MPI_Wtime();
		t4 = t4 - t1;
	}

	//sends process times to process 0	
	MPI_Reduce(&t3,&fastestTime, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
	MPI_Reduce(&t3,&slowestTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	
	//prints run info	 
	if (rank == 0){
		gmp_printf("[%Zd , %Zd]\t", in1, in2);
		printf("Gap: %li\tp: %d\tPrts: %d\tShrt: %.3f\tLng: %.3f\tTot: %.3f\n",final_max_gap,p,p*p, fastestTime, slowestTime, t4);
	}

	MPI_Finalize();

  return 0;
} 

