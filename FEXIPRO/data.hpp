#include "file_input.hpp"
#include <vector>
#include <unordered_map>
#include <map>
using namespace std;


// definition of data
class data {
public:

	unsigned int identifier = 0;
	std::vector<double> vec;
	std::vector<double> vec_svd;
	std::vector<double> vec_svd_norm;
	std::vector<double> vec_svd_int;

	std::vector<double> i2v_vec;
	double norm = 0;
	double norm_part = 0;
	
	std::pair<double, unsigned int> result;

	std::map<double, unsigned int, std::greater<double>> topk;

	
	std::map<double, unsigned int, std::greater<double>> re_topk;

	double time_processing = 0;
	
	double threshold = 0;
	int data_id;


	double min_dist=0;
	double score =0;
	double ip_sum =0;

	// constructor
	data() { }

	// norm computation for item
	void norm_computation_item(const unsigned int dimensionality_part) {

		for (unsigned int i = 0; i < vec.size(); ++i) {
			norm += vec[i] * vec[i];
			if (i >= dimensionality_part) norm_part += vec_svd[i] * vec_svd[i];
		}

		norm = sqrt(norm);
		norm_part = sqrt(norm_part);
	}

	// norm computation for query
	void norm_computation_query(const unsigned int dimensionality_part) {

		for (unsigned int i = 0; i < vec.size(); ++i) norm += vec[i] * vec[i];
		norm = sqrt(norm);
	}

	// svd vec normalization
	void vec_svd_normalization(double max_val) {

		double temp = 0;

		for (unsigned int i = 0; i < vec_svd.size(); ++i) {

			temp = 100 * vec_svd[i] / max_val;
			vec_svd_norm.push_back(temp);
			vec_svd_int.push_back(floor(temp));
		}
	}

	// k-MIPS update
	
	void update_topk(const double ip, const unsigned int id, const unsigned int k) {

		if (ip > threshold) {
			topk.insert({ip,id});
		}

		if (topk.size() > k) {
			auto it = topk.end();
			--it;
			topk.erase(it);
		}

		if (topk.size() == k) {
			auto it = topk.end();
			--it;
			threshold = it->first;
		}
	}
	

	// for norm-based sorting
	bool operator <(const data &d) const { return norm < d.norm; }
	bool operator >(const data &d) const { return norm > d.norm; }
};


// set of query and item vectors
std::vector<data> query_set, item_set;


// definition of sorted list
struct sorted_list {
	unsigned int identifier;
	double coordinate_value;

	bool operator >(const sorted_list &s) const {
		return coordinate_value > s.coordinate_value;
	}
};

std::vector<std::vector<sorted_list>> sorted_list_set;


// input matrix
std::vector<string> split(string& input, char delimiter){
    istringstream stream(input);
    string field;
    std::vector<string> result;
    while (getline(stream, field, delimiter)) {
        result.push_back(field);
    }
    return result;
}



void input_matrix() {

	std::mt19937_64 mt(1);
	std::uniform_real_distribution<> uni_rnd(0,1);

	// file-name
	std::string f_name;

	// variables for query and item vec
	data query, item;

	if (dataset_id == 0 || dataset_id == 1 || dataset_id == 2 || dataset_id == 3) {

		// MF data
		if (dataset_id == 0) f_name = "../dataset/dataset_MF200/netflix_mf-200.txt";
		if (dataset_id == 1) f_name = "../dataset/dataset_MF200/amazon_Movies_and_TV_mf-200.txt";
		if (dataset_id == 2) f_name = "../dataset/dataset_MF200/amazon_Kindle_Store_mf-200.txt";
		if (dataset_id == 3) f_name = "../dataset/dataset_MF200/MovieLens_mf-200.txt";
	}

	std::ifstream ifs_file(f_name);
	if (ifs_file.fail()) {
		std::cout << " file does not exist." << std::endl;
		std::exit(0);
	}
	string line;
    while (getline(ifs_file, line)) {
        std::vector<string> strvec = split(line,' ');
        if(strvec[1]=="F"){
            continue;
        }

        string tmp_str = strvec[0];
        if(tmp_str[0] == 'p'){
			int tmp_id = std::stoi(tmp_str.substr(1));
            strvec.erase(strvec.begin());
            strvec.erase(strvec.begin());

            vector<double> tmp;
            for(int j=0;j<strvec.size();j++){
				double d = std::stod(strvec.at(j));
                tmp.push_back(d);
            }
            query.vec=tmp;
			query.data_id = tmp_id;
			query_set.push_back(query);

			++query.identifier;
        }

        if(tmp_str[0] == 'q'){
			int tmp_id = std::stoi(tmp_str.substr(1));
            strvec.erase(strvec.begin());
            strvec.erase(strvec.begin());

            vector<double> tmp;
            for(int j=0; j<strvec.size();j++){
				double d = std::stof(strvec.at(j));
                tmp.push_back(d);
            }
            item.vec=tmp;
			item.data_id = tmp_id;
			item_set.push_back(item);

			++item.identifier;
            
        }        
    }
}

void input_matrix_i2v(){
	// file-name
	std::string f_name;
	if (dataset_id == 0 || dataset_id == 1 || dataset_id == 2 || dataset_id == 3) {

		// MF data
		if (dataset_id == 0) f_name = "../dataset/dataset_item2vec/netflix_item2vec_d-200.txt";
		if (dataset_id == 1) f_name = "../dataset/dataset_item2vec/amazon_Movie_item2vec_d-200.txt";
		if (dataset_id == 2) f_name = "../dataset/dataset_item2vec/amazon_Kindle_item2vec_d-200.txt";
		if (dataset_id == 3) f_name = "../dataset/dataset_item2vec/MovieLens_item2vec_d-200.txt";
	}
    ifstream ifs(f_name);
    if(!ifs){
        cout<<"Error! File can not be opened"<<endl;
    }

    string line;
	int count=0;
    while (getline(ifs, line)) {
        vector<string> strvec = split(line,' ');
        strvec.erase(strvec.begin());

        vector<double> tmp;
        for(int i=0; i<strvec.size(); i++){
            double d = stod(strvec.at(i));
            tmp.push_back(d);
        }
		item_set[count].i2v_vec =tmp;
        count++;
    }

}

//input parameter
int input_k(){
    ifstream ifs("../parameter/k.txt");
    if(!ifs){
        cout<<"Error! File can not be opened"<<endl;
    }
    string line;
    int k;
    while (getline(ifs, line)) {
        vector<string> strvec = split(line,' ');
        k = stoi(strvec[0]);
    }
    return k;
}

double input_lamda(){
    ifstream ifs("../parameter/lamda.txt");
    if(!ifs){
        cout<<"Error! File can not be opened"<<endl;
    }
    string line;
    double lamda;
    while (getline(ifs, line)) {
        vector<string> strvec = split(line,' ');
        lamda = stod(strvec[0]);
    }
    return lamda;
}
