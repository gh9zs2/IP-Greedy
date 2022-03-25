#include "file_output.hpp"
#include <omp.h>
#include <algorithm>
#include <functional>
#include <random>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/SVD>


Eigen::MatrixXf svd_matrix_for_query;

// maximum scalar
double max_p = 0;


// pre-processing
void pre_processing() {

	std::chrono::system_clock::time_point start = std::chrono::system_clock::now();

	// matrix generation
	Eigen::MatrixXf item_matrix = Eigen::MatrixXf::Zero(dimensionality, item_set.size());
	for (unsigned int i = 0; i < item_set.size(); ++i) {
		for (unsigned int j = 0; j < dimensionality; ++j) item_matrix(j, i) = item_set[i].vec[j];
	}

	// SVD
	Eigen::JacobiSVD<Eigen::MatrixXf> svd(item_matrix, Eigen::ComputeThinU | Eigen::ComputeThinV);
	Eigen::MatrixXf D = svd.singularValues().asDiagonal();
	Eigen::MatrixXf U = svd.matrixU();
	Eigen::MatrixXf V = svd.matrixV();
	for (unsigned int i = 0; i < item_set.size(); ++i) {
		for (unsigned int j = 0; j < dimensionality; ++j) item_set[i].vec_svd.push_back(V(i, j));
	}

	// w decision
	Eigen::MatrixXf DD = Eigen::MatrixXf::Zero(dimensionality, dimensionality);
	double sum = 0;
	double cum_sum = 0;
	for (unsigned int i = 0; i < dimensionality; ++i) {
		cum_sum += D(i, i);
		DD(i, i) = D(i, i);
	}
	for (unsigned int i = 0; i < dimensionality; ++i) {
		sum += D(i, i);
		if (sum / cum_sum >= 0.7) {
			dimensionality_part = i + 1;
			break;
		}
	}

	svd_matrix_for_query = DD * U.transpose();

	// norm computation
	for (unsigned int i = 0; i < query_set.size(); ++i) query_set[i].norm_computation_query(dimensionality_part);
	for (unsigned int i = 0; i < item_set.size(); ++i) item_set[i].norm_computation_item(dimensionality_part);

	// sort by norm in descending order
	std::sort(item_set.begin(), item_set.end(), std::greater<data>());

	// get maximum scalar of P
	for (unsigned int i = 0; i < item_set.size(); ++i) {
		for (unsigned int j = 0; j < dimensionality; ++j) {
			if (item_set[i].vec_svd[j] > max_p) max_p = item_set[i].vec_svd[j];
		}
	}

	// vec_svd normalization
	for (unsigned int i = 0; i < item_set.size(); ++i) item_set[i].vec_svd_normalization(max_p);


	std::chrono::system_clock::time_point end = std::chrono::system_clock::now();

	time_pre_processing = (double)std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
}

// MIP join processing


void mip_join(vector<int> &forcus_id_vec) {

	std::chrono::system_clock::time_point start = std::chrono::system_clock::now();

	#pragma omp parallel num_threads(thread_num)
	{
		#pragma omp for schedule(static) reduction(+:access_count)
		for (unsigned int ii = 0; ii < forcus_id_vec.size(); ++ii) {//query_set.size()
			int i = forcus_id_vec[ii];
			
			std::chrono::system_clock::time_point start1 = std::chrono::system_clock::now();
			///


			// query transformation
			Eigen::MatrixXf q = Eigen::MatrixXf::Zero(dimensionality, 1);
			for (unsigned int j = 0; j < dimensionality; ++j) q(j, 0) = query_set[i].vec[j];
			Eigen::MatrixXf q_svd = svd_matrix_for_query * q;

			// partial norm computation
			double norm_part = 0;
			for (unsigned int j = dimensionality_part; j < dimensionality; ++j) norm_part += q_svd(j, 0) * q_svd(j, 0);
			norm_part = sqrt(norm_part);

			// maximum scalar
			double max_q = 0;

			// vec_svd normalization
			std::vector<double> q_svd_norm;
			for (unsigned int j = 0; j < dimensionality; ++j) if (q_svd(j, 0) > max_q) max_q = q_svd(j, 0);
			for (unsigned int j = 0; j < dimensionality; ++j) q_svd_norm.push_back(100 * q_svd(j, 0) / max_q);

			// integer scaling
			std::vector<double> q_svd_int;
			for (unsigned int j = 0; j < dimensionality; ++j) q_svd_int.push_back(floor(q_svd_norm[j]));

			//double ip_max = 0;
			//unsigned int cnt = 0;


			for (unsigned int j = 0; j < item_set.size(); ++j) {

				bool flag = 1;

				// cauchy schwaltz inequality
				if (query_set[i].norm * item_set[j].norm < query_set[i].threshold) break;

				++access_count;

				double ip_svd = 0;

				// integer-based bounding
				double bound_1 = 0;
				for (unsigned int d = 0; d < dimensionality_part; ++d) bound_1 += (q_svd_int[d] * item_set[j].vec_svd_int[d]) + abs(q_svd_int[d]) + abs(item_set[j].vec_svd_int[d]) + 1;
				bound_1 *= max_p * max_q / 10000;

				if (bound_1 + norm_part * item_set[j].norm_part < query_set[i].threshold) flag = 0;

				if (flag) {

					double bound_2 = 0;
					for (unsigned int d = dimensionality_part; d < dimensionality; ++d) bound_2 += (q_svd_int[d] * item_set[j].vec_svd_int[d]) + abs(q_svd_int[d]) + abs(item_set[j].vec_svd_int[d]) + 1;
					bound_2 *= max_p * max_q / 10000;

					if (bound_1 + bound_2 < query_set[i].threshold) flag = 0;
				}
				
				if (flag) {

					// ip computation with incr
					for (unsigned int d = 0; d < dimensionality; ++d) {

						// incr bounding
						
						if (d == dimensionality_part) {
							if (ip_svd + norm_part * item_set[j].norm_part < query_set[i].threshold) {
								flag = 0;
								break;
							}
						}
						
						ip_svd += q_svd(d, 0) * item_set[j].vec_svd[d];
						
						
					}

					if (flag) ++ip_count;
				}
				

				// result update
				query_set[i].update_topk(ip_svd, j, k);
				
			}
			std::chrono::system_clock::time_point end1 = std::chrono::system_clock::now();
			query_set[i].time_processing = (double)std::chrono::duration_cast<std::chrono::microseconds>(end1 - start1).count();

		}
	}

	std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
	time_join_processing = (double)std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
}
