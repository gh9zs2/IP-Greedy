#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <random>
#include <time.h>


// dataset identifier
unsigned int dataset_id = 3;//0:netflix, 1:amazon_M, 2:amazon_K, 3:MovieLens

// data dimensionality
unsigned int dimensionality = 200;

// dataset sampling rate
double sampling_rate = 1.0;

// top-k
unsigned int k =10;//5,10,15,20

double lamda =0.5;

// #threads
unsigned int thread_num = 1;

// system parameter
unsigned int dimensionality_part = dimensionality / 2;
const unsigned int coord_dimensionality = 5;


// get current time
void get_current_time() {

	time_t t = time(NULL);
	printf(" %s\n", ctime(&t));
}

