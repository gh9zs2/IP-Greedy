#include "Datainput.h"


double normfun(vector<double> &vec){
    double norm_tmp=0;
    for(int i=0;i<vec.size();i++){
        norm_tmp=norm_tmp + vec[i]*vec[i];
    }
    double norm =sqrt(norm_tmp);
    return norm;
}

double ip_fun(vector<double> &a, vector<double> &b){
    double ip =0;
    for(int i=0; i< a.size(); i++){
        ip = ip + a[i]*b[i];
    }
    return ip;
}

double dist_fun(vector<double> &a, vector<double> &b){
    double tmp=0;
    for(int i=0; i< a.size(); i++){
        tmp += pow((a[i]-b[i]) ,2);
    }
    return sqrt(tmp);
}

void Init_dist_ralations(int data_number, vector<unordered_map<int, double> > &dist_relataions){
    for(int i=0; i<data_number; i++){
        unordered_map<int, double> tmp;
        dist_relataions.push_back(tmp);
    }
}

double get_scale(vector<Item> &items){
    vector<double> max_scalar_vec = items[0].i2v_vec;
    vector<double> min_scalar_vec = items[0].i2v_vec;

    for(int i=1; i< items.size(); i++){
        for(int j=0; j< items[i].i2v_vec.size(); j++){
            if(max_scalar_vec[j] < items[i].i2v_vec[j]){
                max_scalar_vec[j] = items[i].i2v_vec[j];
            }

            if(min_scalar_vec[j] > items[i].i2v_vec[j]){
                min_scalar_vec[j] = items[i].i2v_vec[j];
            }
        }
    }
    vector<double> tmp;
    for(int i=0; i < max_scalar_vec.size(); i++){
        tmp.push_back( max_scalar_vec[i] - min_scalar_vec[i]);
    }
    double s = normfun(tmp);
    return s;
}

void compute_items_norm_mf(vector<Item> &items){
    for(int i=0; i< items.size(); i++){
        double norm = normfun(items[i].mf_vec);
        items[i].mf_norm = norm;
    }
}


void compute_items_norm_i2v(vector<Item> &items){
    for(int i=0; i< items.size(); i++){
        double norm = normfun(items[i].i2v_vec);
        items[i].i2v_norm = norm;
    }
}


void Init_items_rev(vector<Item> &items){
    for(int i=0; i < items.size(); i++){
        items[i].ip = 0;
        items[i].flag = 0;
        items[i].min_dist= 10000000;
    }
}


double compute_time(chrono::system_clock::time_point start, chrono::system_clock::time_point end){
    auto time = end - start;
    double msec = std::chrono::duration_cast<chrono::microseconds>(time).count();
    return (double)msec/1000;
}


void greedy_search_rev(vector<double> &query, int k, double lamda, double s, vector<Item> &items, Result_data &result){
    chrono::system_clock::time_point  start1, end1, start2, end2, start, end;

    start1 = chrono::system_clock::now();
    start2 = chrono::system_clock::now();

    for(int i=0; i < items.size(); i++){
        items[i].ip = ip_fun(query, items[i].mf_vec);
    }
    end1 = chrono::system_clock::now();
    result.msec1 = compute_time(start1, end1);

    start = chrono::system_clock::now();
    sort(items.begin(), items.end(), ip_comp);
    
    int start_id = 0;
    result.answer_id.push_back(start_id);
    items[start_id].flag=1;
    result.ip_sum += items[start_id].ip;

    end2 = chrono::system_clock::now();
    result.time_vec.push_back(compute_time(start2, end2));

   
    for(int i=1; i<k; i++){
        Result_tmp result_tmp;
        start2 = chrono::system_clock::now();
        
        for(int j=0; j<items.size(); j++){
            if(items[j].flag ==1){
                continue;
            }

            //dist
            double min_dist = result.dist_min;

            int partner_id = result.answer_id[i-1];
            double a = dist_fun(items[j].i2v_vec, items[partner_id].i2v_vec);

            if(a < items[j].min_dist){
                items[j].min_dist = a;
            }

            if(items[j].min_dist < min_dist){
                min_dist = items[j].min_dist;
            }
        
            double score_tmp = lamda * items[j].ip + s*(1-lamda) * min_dist;
            if(result_tmp.score < score_tmp){
                result_tmp.score = score_tmp;
                result_tmp.id = j;
                result_tmp.dist_min = min_dist;
                result_tmp.ip = items[j].ip;
            }  
        }
        //update
        result.answer_id.push_back(result_tmp.id);
        items[result_tmp.id].flag = 1;
        result.ip_sum += result_tmp.ip;
        result.dist_min = result_tmp.dist_min;

        end2 = chrono::system_clock::now();
        result.time_vec.push_back(compute_time(start2, end2));
    }

    end = chrono::system_clock::now();
    result.msec2 = compute_time(start, end);

    result.final_score = (double)(lamda * result.ip_sum )/k  + s * (1-lamda) * result.dist_min;
    
}


void record_fun_double(int user_id, vector<vector<double> > &record, vector<double> &record_tmp){
    vector<double> tmp;
    tmp.push_back(user_id);
    for(int i=0; i < record_tmp.size(); i++){
        tmp.push_back(record_tmp[i]);
    }
    record.push_back(tmp);
}
void record_fun_int(int user_id, vector<vector<double> > &record, vector<int> &record_tmp){
    vector<double> tmp;
    tmp.push_back(user_id);
    for(int i=0; i < record_tmp.size(); i++){
        tmp.push_back(record_tmp[i]);
    }
    record.push_back(tmp);
}

void compute_file(string output_file_name, vector<vector<double> > &record){
    ofstream file(output_file_name);
    for(int i=0; i< record.size();i++){
        for(int j=0; j<record[i].size()-1;j++){
            file<< record[i][j];
            file<<',';
        }
        file<< record[i][ record[i].size() -1];
        file <<'\n';
    }
    file.close();
}
