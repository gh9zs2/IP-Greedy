#include "data.hpp"
#include <chrono>
#include <unistd.h>


// variable for time measure
std::chrono::system_clock::time_point start, end, init;

// computation time
double time_pre_processing = 0;
double time_join_processing = 0;
unsigned long access_count = 0;
unsigned long ip_count = 0;


// result
std::map<double, std::pair<unsigned int, unsigned int>, std::greater<double>> topk;

double normfun(vector<double> &vec){
    double norm_tmp=0;
    for(int i=0;i<vec.size();i++){
        norm_tmp=norm_tmp + vec[i]*vec[i];
    }
    double norm =sqrt(norm_tmp);
    return norm;
}

double compute_scale(){
    vector<double> max_scalar_vec = item_set[0].i2v_vec;
    vector<double> min_scalar_vec = item_set[0].i2v_vec;

    for(int i=1; i< item_set.size(); i++){
        for(int j=0; j< item_set[i].i2v_vec.size(); j++){
            if(max_scalar_vec[j] < item_set[i].i2v_vec[j]){
                max_scalar_vec[j] = item_set[i].i2v_vec[j];
            }

            if(min_scalar_vec[j] > item_set[i].i2v_vec[j]){
                min_scalar_vec[j] = item_set[i].i2v_vec[j];
            }
        }
    }
    vector<double> tmp;
    for(int i=0; i < max_scalar_vec.size(); i++){
        tmp.push_back( max_scalar_vec[i] - min_scalar_vec[i]);
    }
    double scale = normfun(tmp);
    return scale;   
}

double dist_fun(vector<double> &a, vector<double> &b){
    double tmp=0;
    for(int i=0; i< a.size(); i++){
        tmp += pow((a[i]-b[i]) ,2);
    }
    return sqrt(tmp);
}

double ip_fun(vector<double> &a, vector<double> &b){
	double tmp =0;
	for(int i=0; i< a.size(); i++){
		tmp += a[i] * b[i];
	}
	return tmp;
}


void div_check(double scale,  int forcus_id){

	vector<int> answer_id_vec;
	for(auto itr = query_set[forcus_id].topk.begin(); itr != query_set[forcus_id].topk.end(); ++itr){
		int tmp_id = itr->second;
		answer_id_vec.push_back(tmp_id);
	}
	//ip_sum
	for(int i=0; i< answer_id_vec.size(); i++){
		int ans_id = answer_id_vec[i];
		query_set[forcus_id].ip_sum += ip_fun(query_set[forcus_id].vec, item_set[ans_id].vec);
	}

	//min_dist
	double min_dist = 10000000;
	for(int i=0; i < answer_id_vec.size(); i++){
		int left =answer_id_vec[i];
		for(int j=i+1; j < answer_id_vec.size(); j++){
			int right = answer_id_vec[j];
			double dist_tmp = dist_fun(item_set[left].i2v_vec, item_set[right].i2v_vec);

			if(dist_tmp < min_dist){
				min_dist = dist_tmp;
			}
		}
	}
	query_set[forcus_id].min_dist = min_dist;
	query_set[forcus_id].score = (double)(lamda * query_set[forcus_id].ip_sum)/k + scale * (1-lamda) * query_set[forcus_id].min_dist;



}





void re_output_result(vector<int> &forcus_id_vec) {

    // file-name
	std::string f_name ;
	if (dataset_id == 0) f_name = "FEXIPRO_netflix_";
	if (dataset_id == 1) f_name = "FEXIPRO_amazon_M_";
	if (dataset_id == 2) f_name = "FEXIPRO_amazon_K_";
	if (dataset_id == 3) f_name = "FEXIPRO_MovieLens_";

	f_name += "k-" + std::to_string(k) + "_lam-" + std::to_string(lamda) +  ".csv";
	std::ofstream file(f_name);

	double scale = compute_scale();
	cout << scale << endl;

	for(int ii=0; ii < forcus_id_vec.size(); ii++){
		int i = forcus_id_vec[ii];
		div_check((double)5/scale, i);
		file << i << "," <<query_set[i].score <<"," << query_set[i].min_dist
		<<","<< (double)query_set[i].time_processing /1000 << ",";
		int tmp=0;
		

		for(auto itr = query_set[i].topk.begin(); itr != query_set[i].topk.end(); ++itr){
			int tmp_id = itr->second;
			if(tmp != k-1){
				file << item_set[tmp_id].data_id;
				file << ",";
			}else{
				file << item_set[tmp_id].data_id;
			}
			tmp++;
		}
		file << "\n"; 
		
	}
	file.close();
	std::cout << "finish" << std::endl;	


}
