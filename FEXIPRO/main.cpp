#include "fexipro.hpp"

mt19937 mt(1);
int main() {

    //parameter
    k = input_k();
    lamda = input_lamda();

	// file input
	input_matrix();
    input_matrix_i2v();

    // print the current time
    get_current_time();

    std::cout << " --------------------\n";
	std::cout << " dataset id: " << dataset_id << "\n";
	std::cout << " query-set cardinality: " << query_set.size() << "\n";
    std::cout << " item-set cardinality: " << item_set.size() << "\n";
	std::cout << " dataset dimensionality: " << dimensionality << "\n";
    std::cout << " sampling rate: " << sampling_rate << "\n";
    std::cout << " k: " << k << "\n";
    std::cout << " number of threads: " << thread_num << "\n";
	std::cout << " --------------------\n\n";


    // pre-processing
	pre_processing();
	std::cout << "Pre-processing time: " << time_pre_processing / 1000 << "[millisec]\n\n";

    // MIP-join
    std::uniform_int_distribution<> rand_query_data(0, query_set.size()-1); 
    vector<int> forcus_id_vec;
    for(int i=0; i <100; i++){
        int user_id = rand_query_data(mt);
        forcus_id_vec.push_back(user_id);
    }
    mip_join(forcus_id_vec);



    re_output_result(forcus_id_vec);


    return 0;
}