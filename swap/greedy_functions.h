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
        items[i].dist_relations.clear();
    }
}


double compute_time(chrono::system_clock::time_point start, chrono::system_clock::time_point end){
    auto time = end - start;
    double msec = std::chrono::duration_cast<chrono::microseconds>(time).count();
    return (double)msec/1000;
}

double search_min_dist(vector<int> &id_vec, vector<Item> &items){
    double min_dist = 10000000;
    for(int i=0; i < id_vec.size(); i++){
        int left = id_vec[i];
        for(int j= i+1; j < id_vec.size(); j++){
            int right = id_vec[j];
            double tmp = items[left].dist_relations[right];

            if(tmp ==0){
                double dist = dist_fun(items[left].i2v_vec, items[right].i2v_vec);
                items[left].dist_relations[right] = dist;
                items[right].dist_relations[left] = dist;
                tmp = dist;
            }

            if(tmp < min_dist){
                min_dist = tmp;
            }
        }
    }
    return min_dist;
}


double search_min_dist_skip(vector<int> &id_vec, vector<Item> &items, int skip_id){
    double min_dist = 10000000;
    for(int i=0; i < id_vec.size(); i++){
        if(i == skip_id){
            continue;
        }
        int left = id_vec[i];
        for(int j= i+1; j < id_vec.size(); j++){
            if(j == skip_id){
                continue;
            }
            int right = id_vec[j];
            double tmp = items[left].dist_relations[right];

            if(tmp ==0){
                double dist = dist_fun(items[left].i2v_vec, items[right].i2v_vec);
                items[left].dist_relations[right] = dist;
                items[right].dist_relations[left] = dist;
                tmp = dist;
            }

            if(tmp < min_dist){
                min_dist = tmp;
            }
        }
    }
    return min_dist;
}


double search_mindist_DataToSet(int forcus, vector<int> &id_vec, vector<Item> &items, int skip_id){
    double dist_min = 10000000;
    for(int i =0; i< id_vec.size();i++){
        if(i == skip_id){
            continue;
        }
        int right = id_vec[i];
        double tmp = items[forcus].dist_relations[right];
        if(tmp ==0){
            tmp = dist_fun(items[forcus].i2v_vec, items[right].i2v_vec);
            items[forcus].dist_relations[right] = tmp;
            items[right].dist_relations[forcus] =tmp;
        }
        if(tmp < dist_min){
            dist_min = tmp;
        }
    }
    return dist_min;
}

void swap_data(vector<double> &query, int k, double lamda, double s, vector<Item> &items, Result_data &result){
    chrono::system_clock::time_point  start1, end1, start2, end2;

    start1 = chrono::system_clock::now();
    
    for(int i=0; i<items.size(); i++){
        items[i].ip = ip_fun(query, items[i].mf_vec);
    }
    
    for(int i=0; i < k; i++){
        items[i].flag =1;
        result.answer_id.push_back(i);
        result.ip_sum += items[i].ip;
    }
    
    result.dist_min = search_min_dist(result.answer_id, items);
    result.final_score = (double)(lamda * result.ip_sum )/k  + s * (1-lamda) * result.dist_min;

    end1 = chrono::system_clock::now();
    result.time_vec.push_back(compute_time(start1, end1));

    for(int i=0; i < k-1; i++){
        start2 = chrono::system_clock::now();
        vector<double> remain_ip_sum;
        vector<double> remain_dist_min;
        for(int j=0; j < k; j++){
            int tmp_id = result.answer_id[j];
            remain_ip_sum.push_back(result.ip_sum - items[tmp_id].ip);
            remain_dist_min.push_back(search_min_dist_skip(result.answer_id, items, j));
        }
        
        Best_pair best_pair;
        best_pair.score = result.final_score;
        int flag=1;
        for(int j=0; j<items.size(); j++){
            if(items[j].flag ==1){
                continue;
            }

            for(int l=0; l <k; l++){
                double dist_min = min(remain_dist_min[l], search_mindist_DataToSet(j, result.answer_id, items, l));//
                double ip_sum = remain_ip_sum[l] + items[j].ip;
                double new_score = (double)(lamda * ip_sum)/k  + s * (1-lamda) * dist_min;
                if(best_pair.score < new_score){
                    best_pair.score = new_score;
                    best_pair.ip_sum = ip_sum;
                    best_pair.dist_min = dist_min;
                    best_pair.in_to_out = make_pair(j,l);
                    flag =0;
                }
            }
        }
        if(flag ==1){
            //cout << "early break " << endl;
            end2 = chrono::system_clock::now();
            result.time_vec.push_back(compute_time(start2, end2));
            break;
        }

        //update
        result.final_score = best_pair.score;
        result.dist_min = best_pair.dist_min;
        result.ip_sum = best_pair.ip_sum;
        int out_id = result.answer_id[best_pair.in_to_out.second];
        items[out_id].flag =0;
        result.answer_id[best_pair.in_to_out.second] = best_pair.in_to_out.first;
        items[best_pair.in_to_out.first].flag =1;

        end2 = chrono::system_clock::now();
        result.time_vec.push_back(compute_time(start2, end2));
    }

    
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
    ofstream file(output_file_name);//
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