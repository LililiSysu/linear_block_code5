/*****************************************************************//**
 * \file	MRIP_PSM.h
 * \brief	use prestored Matrix to lower the complexity of Gaussian elimination, and parallized it.
 *			the file contains different versions for one solution
 * 
 * \author	26259
 * \date	June 2023
 *********************************************************************/

#pragma once
#include"./reprocess_common.h"
using namespace std;

// the following function is evolved, suitable for writing in article

// equal distance partition
class PreStored_Matrix_red_equal_dist {

public:
	int n;
	int k;
	int red;
	Matrix<Matrix<int>> permute_set;
	Matrix<Matrix<int>> sorted_info_set;				// for select_G_set, to be removed
	Matrix<int> k_MRP_sorted;							// for select_G_set
	Matrix<int> k_distribution_in_n;					// for select_G_set_continous_count
	Matrix<Matrix<info_bin_continous_record>> continous_info_bin;	// for select_G_set_continous_count
	Matrix<int> common_pos_num;
	Matrix<Matrix<GF2>> G_set;							// to be changed into parity-check part

	Basis_update BU;

	// fetch it through G_taregt, by permute the columns of G_target with columns of redundency set, then do Gaussian elemination
	// return: G_target, permute_target and the iterations that BU has processed
	int get_MRIP_sys_G(const Matrix<int>& reliability_to_col_ind, Matrix<GF2>& G_target, Matrix<int>& permute_target) {
		// in this function we start from redundant set and find the least reliability number of the column index
		// this is of less hardware complexity, but computation cost is larger than 'get_MRIP_sys_G'

		//select_G_set(reliability_to_col_ind);
		int target_ind = select_G_set_continous_count(reliability_to_col_ind);

		//cout << "common_pos_num(target_ind) = " << common_pos_num(target_ind) << endl;
		G_target = G_set(target_ind);
		permute_target = permute_set(target_ind);			// initialize the answer

		//cout << "G_target" << G_target;
		//cout << "permute_target" << permute_target;

		return BU.solve(reliability_to_col_ind, G_target, permute_target);
	}

	// select G that systematic position contains most position in sorted_pos
	// sorted_pos should indicates positions with ordered reliability from most to least
	int select_G_set_accurate(const Matrix<int>& reliability_to_col_ind) {
		// first sort the reliability_to_col_ind

		//cout << "reliability_to_col_ind" << reliability_to_col_ind;

		reliability_to_col_ind.get_part(0, 0, 0, k - 1, k_MRP_sorted);
		k_MRP_sorted.sort('<');

		//cout << "k_MRP_sorted" << k_MRP_sorted;

		// count the match position for each permute set, that is same as sorted_info_set
		common_pos_num.reset(0);
		for (int i = 0, imax = sorted_info_set.size(); i < imax; ++i) {
			for (int sorted_info_set_ind = 0, k_MRP_sorted_ind = 0; sorted_info_set_ind < k && k_MRP_sorted_ind < k;) {
				if (sorted_info_set(i)(sorted_info_set_ind) < k_MRP_sorted(k_MRP_sorted_ind)) {
					sorted_info_set_ind++;
				}
				else if (sorted_info_set(i)(sorted_info_set_ind) > k_MRP_sorted(k_MRP_sorted_ind)) {
					k_MRP_sorted_ind++;
				}
				else {
					sorted_info_set_ind++;
					k_MRP_sorted_ind++;
					common_pos_num(i)++;
				}
			}
		}
		//cout << "common_pos_num" << common_pos_num;

		// set the G_target to the G_set(index), index = arg max (common_pos_num)
		int target_ind = common_pos_num.max_ele_ind();
		return target_ind;
	}

	// select G that systematic position contains most position in sorted_pos
	// sorted_pos should indicates positions with ordered reliability from most to least
	int select_G_set_continous_count(const Matrix<int>& reliability_to_col_ind) {

		// set the distribution of MRP
		k_distribution_in_n.reset(0);		// the first position is sure to be 0
		for (int i = 0; i < k; ++i) {
			k_distribution_in_n(reliability_to_col_ind(i) + 1) = 1;
		}
		//cout << "k_distribution_in_n" << k_distribution_in_n;

		// accumulate the MRP
		for (int i = 1; i <= n; ++i) {
			k_distribution_in_n(i) += k_distribution_in_n(i - 1);		// maybe need to think of parallel implementation
		}
		//cout << "k_distribution_in_n (accumulate)" << k_distribution_in_n;

		// count the common info number for each info set in  permute_set, not absolutely accurate although
		common_pos_num.reset(0);
		for (int i = 0, imax = continous_info_bin.size(); i < imax; i++) {
			for (int j = 0, jmax = continous_info_bin(i).size(); j < jmax; ++j) {
				// this will save a lot complexity for selecting G 
				common_pos_num(i) += \
					k_distribution_in_n(continous_info_bin(i)(j).end_ind + 1) - k_distribution_in_n(continous_info_bin(i)(j).start_ind);
				// anding end_ind by one for accumulated the MRP
			}
		}

		// first sort the reliability_to_col_ind

		//cout << "common_pos_num" << common_pos_num;

		// set the G_target to the G_set(index), index = arg max (common_pos_num)
		int target_ind = common_pos_num.max_ele_ind();
		return target_ind;
	}

	/**
	 * .initialize G_set by a generator Matrix of the code, this is offline computation, once for a code, complexity doesn't matter
	 *  but storage does matter. we only consider G such that its row number is no greater than half of its column number
	 *
	 * \param G: generator matrix or parity matrix of a code, we base on this to generate systematic matrix associated with 'MRB'
	 * \param info_partition_level: e.g., 1 means info set not partitioned, 2 means that info set is partitioned into large_bin_num parts,
	 *								3 menas that info set is partitioned into large_bin_num^2 parts, etc.
	 */
	PreStored_Matrix_red_equal_dist(const Matrix<GF2>& G, int partition_patterns_num) : n(G.col()), k(G.row()), red(n - k), BU(n, k) {

		cout << "k = " << k << endl;
		cout << "n - k = " << n - k << endl;

		// we increase the info_partition_num from 1 until the bin_len = 1
		Matrix<int> info_bin_lb(1, k, 'v');		// info bins' length lower bound
		Matrix<int> info_bin_res(1, k, 'v');	// info bins' residual
		Matrix<int> distance_num(1, k, 'v');	// which is partition number
		info_bin_lb.push_back(k);				// firstly push back no partition info lenght, that is k
		info_bin_res.push_back(0);
		distance_num.push_back(1);
		int partition_num = 2;
		while (partition_num <= k && partition_num <= partition_patterns_num) {
			int info_len_lb_tmp = k / partition_num;
			int info_len_res_tmp = k % partition_num;
			if (info_len_lb_tmp == info_bin_lb.back()) {
				// the info length not changed, switch to info length decreasing mode
				info_bin_lb.pop_back();
				info_bin_res.pop_back();
				distance_num.pop_back(); 
				partition_num--;
				break;
			}
			else {
				// the info length has changed, store the info length lower bound
				info_bin_lb.push_back(info_len_lb_tmp);
				info_bin_res.push_back(info_len_res_tmp);
				distance_num.push_back(partition_num);		// partition number is totally equals to distance number
			}
			partition_num++;
		}
		
		// continue to the info length decreasing step
		for (int i = info_bin_lb.back() - 1; i >= 1 && partition_num <= partition_patterns_num; --i) {
			int info_len_res_tmp = k % i;		// this is magic
			int distance_num_tmp = k / i;
			if (info_len_res_tmp == distance_num_tmp) {
				continue;
			}
			info_bin_lb.push_back(i);
			info_bin_res.push_back(info_len_res_tmp);
			distance_num.push_back(distance_num_tmp);
			partition_num++;
		}
		cout << "info_bin_lb" << info_bin_lb;
		cout << "info_bin_res" << info_bin_res;
		cout << "distance_num" << distance_num;

		int set_num = info_bin_lb.size();
		cout << "set_num = " << set_num << endl;

		Matrix<int> distance_lb(1, set_num, 'v');
		Matrix<int> distance_res(1, set_num, 'v');
		for (int i = 0; i < set_num; ++i) {
			distance_lb.push_back((n - k) / distance_num(i));
			distance_res.push_back((n - k) % distance_num(i));
		}
		cout << "distance_lb" << distance_lb;
		cout << "distance_res" << distance_res;

		// this is simpler, the simpler the more efficient

		continous_info_bin.resize(1, set_num, false);
		permute_set.resize(1, set_num, false);
		permute_set.reset(Matrix<int>(1, n, 'v'));

		cout << "----------------" << endl;
		// decide each permute_set and continuous info bin
		Matrix<int> info_bin_specific(1, k);
		Matrix<int> distance_specific(1, k);
		for (int i = 0; i < set_num; ++i) {
			// construct info_bin_specific
			info_bin_specific.resize(1, 0, false);
			for (int j = 0, jmax = distance_num(i) - info_bin_res(i); j < jmax; ++j) {
				info_bin_specific.push_back(info_bin_lb(i));
			}
			for (int j = 0; j < info_bin_res(i); ++j) {
				info_bin_specific.push_back(info_bin_lb(i) + 1);
			}
			cout << "info_bin_specific" << info_bin_specific;

			// construct distance_specific
			distance_specific.resize(1, 0, false);
			for (int j = 0, jmax = distance_num(i) - distance_res(i); j < jmax; ++j) {
				distance_specific.push_back(distance_lb(i));
			}
			for (int j = 0; j < distance_res(i); ++j) {
				distance_specific.push_back(distance_lb(i) + 1);
			}
			cout << "distance_specific" << distance_specific;

			// construct continous_info_bin
			continous_info_bin(i).resize(1, distance_num(i), false);
			continous_info_bin(i)(0) = info_bin_continous_record(0, info_bin_specific(0) - 1);
			for (int j = 1; j < distance_num(i); ++j) {
				continous_info_bin(i)(j) = info_bin_continous_record(\
					continous_info_bin(i)(j - 1).end_ind + distance_specific(j - 1) + 1, \
					continous_info_bin(i)(j - 1).end_ind + distance_specific(j - 1) + info_bin_specific(j));
			}

			// construct permute_set from continous_info_bin
			for (int j = 0, jmax = continous_info_bin(i).size(); j < jmax; ++j) {
				for (int p = continous_info_bin(i)(j).start_ind; p <= continous_info_bin(i)(j).end_ind; ++p) {
					permute_set(i).push_back(p);
				}			
			}
		}

		cout << "continous_info_bin" << continous_info_bin;
		cout << "permute_set" << permute_set;

		// add redundancy positions
		for (int i = 0; i < set_num; ++i) {
			int permute_set_i_ind = 0;
			int info_bin_res_collect_i_ind = 0;
			for (int j = 0; j < n; ++j) {
				if (permute_set_i_ind < k && j == permute_set(i)(permute_set_i_ind)) {
					++permute_set_i_ind;
				}
				else {
					permute_set(i).push_back(j);
				}
			}
		}

		cout << "permute_set (add red)" << permute_set;

		// generate set_num systematic matrices
		G_set.resize(1, set_num, false);
		for (int i = 0; i < set_num; ++i) {
			G_set(i) = G;
			G_set(i).permute_col(permute_set(i));
			//cout << " ---------- " << endl;
			//cout << "permute_set(i)" << permute_set(i);
			//cout << "G_set(i)" << G_set(i);

			/*G_set(i).row_transformation_to_up_triangle();
			Matrix<int> permute_record = G_set(i).col_permute_to_full_rank_on_left();
			G_set(i).row_transformation_left_up_triangle_to_identity();*/

			Matrix<int> permute_record = G_set(i).GE_left_identity_4_GF2();
			permute_set(i).permute(permute_record);

			//cout << "permute_set(i) after" << permute_set(i);
			//cout << "G_set(i) after" << G_set(i);
		}

		sorted_info_set.resize(1, set_num, false);
		for (int i = 0; i < set_num; ++i) {

			// store the first k value of permute_set into sorted_info_set, where the latter is sorted index for choosing best G_set		
			sorted_info_set(i) = permute_set(i).get_part(0, 0, 0, k - 1);
			sorted_info_set(i).sort('<');
		}

		k_MRP_sorted.resize(1, k, false);
		k_distribution_in_n.resize(1, n + 1, false);
		common_pos_num.resize(1, set_num, false);
	}
};

// equal distance partition, filling all codeword length
class PreStored_Matrix_red_equal_dist_fill {		// todo: add the accurate count strategy back and simulate

public:
	int n;
	int k;
	int red;
	Matrix<Matrix<int>> permute_set;
	Matrix<int> k_distribution_in_n;					// for select_G_set_continous_count
	Matrix<Matrix<info_bin_continous_record>> continous_info_bin;	// for select_G_set_continous_count
	int max_overflow;
	Matrix<int> common_pos_num;
	Matrix<Matrix<GF2>> G_set;							// to be changed into parity-check part

	Basis_update BU;

	// fetch it through G_taregt, by permute the columns of G_target with columns of redundency set, then do Gaussian elemination
	// return: G_target, permute_target and the max common on G_target elements on start
	int get_MRIP_sys_G(const Matrix<int>& reliability_to_col_ind, Matrix<GF2>& G_target, Matrix<int>& permute_target) {
		// in this function we start from redundant set and find the least reliability number of the column index
		// this is of less hardware complexity, but computation cost is larger than 'get_MRIP_sys_G'

		//select_G_set(reliability_to_col_ind);
		int target_ind = select_G_set_continous_count(reliability_to_col_ind);

		//cout << "common_pos_num(target_ind) = " << common_pos_num(target_ind) << endl;
		G_target = G_set(target_ind);
		permute_target = permute_set(target_ind);			// initialize the answer


		//cout << "G_target" << G_target;
		//cout << "permute_target" << permute_target;

		BU.solve(reliability_to_col_ind, G_target, permute_target);

		return common_pos_num(target_ind);
	}

	// select G that systematic position contains most position in sorted_pos
	// sorted_pos should indicates positions with ordered reliability from most to least
	int select_G_set_continous_count(const Matrix<int>& reliability_to_col_ind) {

		// set the distribution of MRP
		k_distribution_in_n.reset(0);		// the first position is sure to be 0
		for (int i = 0; i < k; ++i) {
			k_distribution_in_n(reliability_to_col_ind(i) + 1) = 1;
		}
		// set the post fix of k_distribution_in_n with length k
		for (int i = 1; i <= max_overflow; ++i) {
			k_distribution_in_n(n + i) = k_distribution_in_n(i);
		}
		//cout << "k_distribution_in_n" << k_distribution_in_n;

		// accumulate the MRP
		for (int i = 1, imax = n + max_overflow; i <= imax; ++i) {
			k_distribution_in_n(i) += k_distribution_in_n(i - 1);		// maybe need to think of parallel implementation
		}
		//cout << "k_distribution_in_n (accumulate)" << k_distribution_in_n;

		// count the common info number for each info set in  permute_set, not absolutely accurate although
		common_pos_num.reset(0);
		for (int i = 0, imax = continous_info_bin.size(); i < imax; i++) {
			for (int j = 0, jmax = continous_info_bin(i).size(); j < jmax; ++j) {
				// this will save a lot complexity for selecting G 
				common_pos_num(i) += \
					k_distribution_in_n(continous_info_bin(i)(j).end_ind + 1) - k_distribution_in_n(continous_info_bin(i)(j).start_ind);
				// anding end_ind by one for accumulated the MRP
			}
		}

		// first sort the reliability_to_col_ind

		//cout << "common_pos_num" << common_pos_num;

		// set the G_target to the G_set(index), index = arg max (common_pos_num)
		int target_ind = common_pos_num.max_ele_ind();
		return target_ind;
	}

	/**
	 * .initialize G_set by a generator Matrix of the code, this is offline computation, once for a code, complexity doesn't matter
	 *  but storage does matter. we only consider G such that its row number is no greater than half of its column number
	 *
	 * \param G: generator matrix or parity matrix of a code, we base on this to generate systematic matrix associated with 'MRB'
	 * \param info_partition_level: e.g., 1 means info set not partitioned, 2 means that info set is partitioned into large_bin_num parts,
	 *								3 menas that info set is partitioned into large_bin_num^2 parts, etc.
	 */
	PreStored_Matrix_red_equal_dist_fill(const Matrix<GF2>& G, int partition_patterns_num) : n(G.col()), k(G.row()), red(n - k), BU(n, k) {

		//cout << "k = " << k << endl;
		//cout << "n - k = " << n - k << endl;

		// we increase the info_partition_num from 1 until the bin_len = 1
		Matrix<int> info_bin_lb(1, k, 'v');		// info bins' length lower bound
		Matrix<int> info_bin_res(1, k, 'v');	// info bins' residual
		Matrix<int> distance_num(1, k, 'v');	// which is partition number
		info_bin_lb.push_back(k);				// firstly push back no partition info lenght, that is k
		info_bin_res.push_back(0);
		distance_num.push_back(1);
		int partition_num = 2;
		while (partition_num <= k && partition_num <= partition_patterns_num) {
			int info_len_lb_tmp = k / partition_num;
			int info_len_res_tmp = k % partition_num;
			if (info_len_lb_tmp == info_bin_lb.back()) {
				// the info length not changed, switch to info length decreasing mode
				info_bin_lb.pop_back();
				info_bin_res.pop_back();
				distance_num.pop_back();
				partition_num--;
				break;
			}
			else {
				// the info length has changed, store the info length lower bound
				info_bin_lb.push_back(info_len_lb_tmp);
				info_bin_res.push_back(info_len_res_tmp);
				distance_num.push_back(partition_num);		// partition number is totally equals to distance number
			}
			partition_num++;
		}

		// continue to the info length decreasing step
		for (int i = info_bin_lb.back() - 1; i >= 1 && partition_num <= partition_patterns_num; --i) {
			int info_len_res_tmp = k % i;		// this is magic
			int distance_num_tmp = k / i;
			if (info_len_res_tmp == distance_num_tmp) {
				continue;
			}
			info_bin_lb.push_back(i);
			info_bin_res.push_back(info_len_res_tmp);
			distance_num.push_back(distance_num_tmp);
			partition_num++;
		}
		//cout << "info_bin_lb" << info_bin_lb;
		//cout << "info_bin_res" << info_bin_res;
		//cout << "distance_num" << distance_num;

		int partition_patterns_num_real = info_bin_lb.size();

		Matrix<int> distance_lb(1, partition_patterns_num_real, 'v');
		Matrix<int> distance_res(1, partition_patterns_num_real, 'v');
		for (int i = 0; i < partition_patterns_num_real; ++i) {
			distance_lb.push_back((n - k) / distance_num(i));
			distance_res.push_back((n - k) % distance_num(i));
		}
		//cout << "distance_lb" << distance_lb;
		//cout << "distance_res" << distance_res;

		// this is simpler, the simpler the more efficient

		int n_over_k_ub = (int)ceil(double(n) / k);
		int set_num = partition_patterns_num_real * n_over_k_ub;
		//cout << "set_num = " << set_num << endl;

		continous_info_bin.resize(1, set_num, false);
		permute_set.resize(1, set_num, false);
		permute_set.reset(Matrix<int>(1, n, 'v'));

		//cout << "----------------" << endl;
		// decide each permute_set and continuous info bin
		Matrix<int> info_bin_specific(1, k);
		Matrix<int> distance_specific(1, k);
		max_overflow = 0;
		for (int i = 0; i < partition_patterns_num_real; ++i) {
			// construct info_bin_specific
			info_bin_specific.resize(1, 0, false);
			for (int j = 0, jmax = distance_num(i) - info_bin_res(i); j < jmax; ++j) {
				info_bin_specific.push_back(info_bin_lb(i));
			}
			for (int j = 0; j < info_bin_res(i); ++j) {
				info_bin_specific.push_back(info_bin_lb(i) + 1);
			}
			//cout << "info_bin_specific" << info_bin_specific;

			// construct distance_specific
			distance_specific.resize(1, 0, false);
			for (int j = 0, jmax = distance_num(i) - distance_res(i); j < jmax; ++j) {
				distance_specific.push_back(distance_lb(i));
			}
			for (int j = 0; j < distance_res(i); ++j) {
				distance_specific.push_back(distance_lb(i) + 1);
			}
			//cout << "distance_specific" << distance_specific;

			// construct continous_info_bin
			for (int p = 0; p < n_over_k_ub; ++p) {
				int continous_info_bin_ind = i * n_over_k_ub + p;
				continous_info_bin(continous_info_bin_ind).resize(1, distance_num(i), false);
				continous_info_bin(continous_info_bin_ind)(0) = \
					info_bin_continous_record(p * info_bin_specific(0), p * info_bin_specific(0) + info_bin_specific(0) - 1);
				// continous_info_bin is shifted with p as p * info_bin_specific(0)

				for (int j = 1; j < distance_num(i); ++j) {
					continous_info_bin(continous_info_bin_ind)(j) = info_bin_continous_record(\
						continous_info_bin(continous_info_bin_ind)(j - 1).end_ind + distance_specific(j - 1) + 1, \
						continous_info_bin(continous_info_bin_ind)(j - 1).end_ind + distance_specific(j - 1) + info_bin_specific(j));
				}
				max_overflow = my::max(max_overflow, continous_info_bin(continous_info_bin_ind).back().end_ind);
			}
			//cout << "continous_info_bin" << continous_info_bin;
		}
		max_overflow++;
		max_overflow -= n;

		for (int i = 0; i < set_num; ++i) {
			// construct permute_set from continous_info_bin

			// first check the last continous_info_bin, if it overflow, put its position to the front
			int front_push = continous_info_bin(i).back().end_ind - n + 1;
			front_push = my::max(0, front_push);
			for (int p = 0; p < front_push; ++p) {
				permute_set(i).push_back(p);
			}

			// add the middle continous_info_bin
			for (int j = 0, jmax = continous_info_bin(i).size() - 1; j < jmax; ++j) {

				for (int p = continous_info_bin(i)(j).start_ind, pmax = continous_info_bin(i)(j).end_ind; p <= pmax; ++p) {
					permute_set(i).push_back(p);
				}
			}

			// add the last continous_info_bin
			for (int p = continous_info_bin(i).back().start_ind, pmax = continous_info_bin(i).back().end_ind - front_push; p <= pmax; ++p) {
				permute_set(i).push_back(p);
			}
		}

		//cout << "permute_set" << permute_set;

		// add redundancy positions
		for (int i = 0; i < set_num; ++i) {
			int permute_set_i_ind = 0;
			int info_bin_res_collect_i_ind = 0;
			for (int j = 0; j < n; ++j) {
				if (permute_set_i_ind < k && j == permute_set(i)(permute_set_i_ind)) {
					++permute_set_i_ind;
				}
				else {
					permute_set(i).push_back(j);
				}
			}
		}

		//cout << "permute_set (add red)" << permute_set;

		// generate set_num systematic matrices
		G_set.resize(1, set_num, false);
		for (int i = 0; i < set_num; ++i) {
			G_set(i) = G;
			G_set(i).permute_col(permute_set(i));
			//cout << " ---------- " << endl;
			//cout << "permute_set(i)" << permute_set(i);
			//cout << "G_set(i)" << G_set(i);

			/*G_set(i).row_transformation_to_up_triangle();
			Matrix<int> permute_record = G_set(i).col_permute_to_full_rank_on_left();
			G_set(i).row_transformation_left_up_triangle_to_identity();*/

			Matrix<int> permute_record = G_set(i).GE_left_identity_4_GF2();
			permute_set(i).permute(permute_record);

			//cout << "permute_set(i) after" << permute_set(i);
			//cout << "G_set(i) after" << G_set(i);
		}

		k_distribution_in_n.resize(1, n + max_overflow + 1, false);
		common_pos_num.resize(1, set_num, false);
	}
};

// this is for cycle codes, where generator matrix or parity matrix can shift, 
// e.g., BCH(127,64,21), process on its dual code since rate > 1/2
class PreStored_Matrix_red_equal_dist_cycle {

public:
	int n;
	int k;
	int red;
	Matrix<Matrix<int>> permute_set;
	Matrix<int> k_distribution_in_n;								// for select_G_set_continous_count
	Matrix<Matrix<info_bin_continous_record>> continous_info_bin;	// for select_G_set_continous_count
	Matrix<int> shift_space;										// for select_G_set_continous_count
	Matrix<int> common_pos_num;
	Matrix<int> common_pos_max_shift;
	Matrix<Matrix<GF2>> G_set;										// to be changed into parity-check part, the storage is different
	Matrix<int> natural;

	Basis_update BU;

	// fetch it through G_taregt, by permute the columns of G_target with columns of redundency set, then do Gaussian elemination
	// return: G_target, permute_target and the max common on G_target elements on start
	int get_MRIP_sys_G(const Matrix<int>& reliability_to_col_ind, Matrix<GF2>& G_target, Matrix<int>& permute_target) {
		// in this function we start from redundant set and find the least reliability number of the column index
		// this is of less hardware complexity, but computation cost is larger than 'get_MRIP_sys_G'

		//select_G_set(reliability_to_col_ind);
		target_ind_and_shift opt_info = select_G_set_continous_count(reliability_to_col_ind);
		int target_ind = opt_info.target_ind;
		int shift = opt_info.shift;

		//cout << "target_ind = " << target_ind << endl;
		//cout << "shift = " << shift << endl;

		//cout << "G_set(target_ind)" << G_set(target_ind);
		//cout << "permute_set(target_ind)" << permute_set(target_ind);

		//cout << "common_pos_num(target_ind) = " << common_pos_num(target_ind) << endl;
		G_target = G_set(target_ind);
		natural.col_shift_right_cir(n - shift, permute_target);			// initialize the answer
		permute_target.permute(permute_set(target_ind));

		//cout << "G_target" << G_target;
		//cout << "permute_target" << permute_target;

		BU.solve(reliability_to_col_ind, G_target, permute_target);

		return common_pos_num(target_ind);
	}

	// select G that systematic position contains most position in sorted_pos
	// sorted_pos should indicates positions with ordered reliability from most to least
	target_ind_and_shift select_G_set_continous_count(const Matrix<int>& reliability_to_col_ind) {

		// set the distribution of MRP
		k_distribution_in_n.reset(0);		// the first position is sure to be 0
		for (int i = 0; i < k; ++i) {
			k_distribution_in_n(reliability_to_col_ind(i) + 1) = 1;
		}

		// set the post fix of k_distribution_in_n with length k
		for (int i = 1; i <= k; ++i) {
			k_distribution_in_n(n + i) = k_distribution_in_n(i);
		}
		//cout << "k_distribution_in_n" << k_distribution_in_n;

		// accumulate the MRP
		for (int i = 1, imax = n + k; i <= imax; ++i) {
			k_distribution_in_n(i) += k_distribution_in_n(i - 1);		// maybe need to think of parallel implementation
		}
		//cout << "k_distribution_in_n (accumulate)" << k_distribution_in_n;

		// count the common info number for each info set in  permute_set, not absolutely accurate although
		common_pos_num.reset(0);
		common_pos_max_shift.reset(0);

		//int computation_cost = 0;
		for (int i = 0, imax = continous_info_bin.size(); i < imax; i++) {
			for (int p = 0; p < shift_space(i); ++p) {
				int common_pos_num_tmp = 0;
				for (int j = 0, jmax = continous_info_bin(i).size(); j < jmax; ++j) {
					// this will save a lot complexity for selecting G 
					common_pos_num_tmp += k_distribution_in_n(p + continous_info_bin(i)(j).end_ind + 1) \
						- k_distribution_in_n(p + continous_info_bin(i)(j).start_ind);

					// anding end_ind by one for accumulated the MRP

					//computation_cost++;
				}
				bool update = common_pos_num(i) < common_pos_num_tmp;
				common_pos_num(i) = update ? common_pos_num_tmp : common_pos_num(i);
				common_pos_max_shift(i) = update ? p : common_pos_max_shift(i);
			}
		}

		// first sort the reliability_to_col_ind

		//cout << "common_pos_num" << common_pos_num;
		//cout << "total computation in select_G_set_continous_count = " << computation_cost << endl;

		// set the G_target to the G_set(index), index = arg max (common_pos_num)
		int target_ind = common_pos_num.max_ele_ind();
		int shift_val = common_pos_max_shift(target_ind);
		return target_ind_and_shift(target_ind, shift_val);
	}

	/**
	 * .initialize G_set by a generator Matrix of the code, this is offline computation, once for a code, complexity doesn't matter
	 *  but storage does matter. we only consider G such that its row number is no greater than half of its column number
	 *
	 * \param G: generator matrix or parity matrix of a code, we base on this to generate systematic matrix associated with 'MRB'
	 * \param info_partition_level: e.g., 1 means info set not partitioned, 2 means that info set is partitioned into large_bin_num parts,
	 *								3 menas that info set is partitioned into large_bin_num^2 parts, etc.
	 */
	PreStored_Matrix_red_equal_dist_cycle(const Matrix<GF2>& G, int partition_patterns_num) : n(G.col()), k(G.row()), red(n - k), BU(n, k) {

		//cout << "k = " << k << endl;
		//cout << "n - k = " << n - k << endl;
		int ratio_of_cr = n / k;

		// we increase the info_partition_num from 1 until the bin_len = 1
		Matrix<int> info_bin_lb(1, k, 'v');		// info bins' length lower bound
		Matrix<int> info_bin_res(1, k, 'v');	// info bins' residual
		Matrix<int> distance_num(1, k, 'v');	// which is partition number
		info_bin_lb.push_back(k);				// firstly push back no partition info lenght, that is k
		info_bin_res.push_back(0);
		distance_num.push_back(1);
		int partition_num = 2;
		while (partition_num <= k && partition_num <= partition_patterns_num) {
			int info_len_lb_tmp = k / partition_num;
			int info_len_res_tmp = k % partition_num;
			if (info_len_lb_tmp == info_bin_lb.back()) {
				// the info length not changed, switch to info length decreasing mode
				info_bin_lb.pop_back();
				info_bin_res.pop_back();
				distance_num.pop_back();
				partition_num--;
				break;
			}
			else {
				// the info length has changed, store the info length lower bound
				info_bin_lb.push_back(info_len_lb_tmp);
				info_bin_res.push_back(info_len_res_tmp);
				distance_num.push_back(partition_num);		// partition number is totally equals to distance number
			}
			partition_num++;
		}

		// continue to the info length decreasing step
		for (int i = info_bin_lb.back() - 1; i >= 1 && partition_num <= partition_patterns_num; --i) {
			int info_len_res_tmp = k % i;		// this is magic
			int distance_num_tmp = k / i;
			if (info_len_res_tmp == distance_num_tmp) {
				continue;
			}
			info_bin_lb.push_back(i);
			info_bin_res.push_back(info_len_res_tmp);
			distance_num.push_back(distance_num_tmp);
			partition_num++;
		}
		//cout << "info_bin_lb" << info_bin_lb;
		//cout << "info_bin_res" << info_bin_res;
		//cout << "distance_num" << distance_num;

		int set_num = info_bin_lb.size();
		cout << "set_num = " << set_num << endl;

		Matrix<int> distance_lb(1, set_num, 'v');
		Matrix<int> distance_res(1, set_num, 'v');
		for (int i = 0; i < set_num; ++i) {
			distance_lb.push_back((n - k) / distance_num(i));
			distance_res.push_back((n - k) % distance_num(i));
		}
		//cout << "distance_lb" << distance_lb;
		//cout << "distance_res" << distance_res;

		// this is simpler, the simpler the more efficient

		continous_info_bin.resize(1, set_num, false);
		permute_set.resize(1, set_num, false);

		permute_set.reset(Matrix<int>(1, n, 'v'));

		//cout << "----------------" << endl;
		// decide each permute_set and continuous info bin
		Matrix<int> info_bin_specific(1, k);
		Matrix<int> distance_specific(1, k);
		for (int i = 0; i < set_num; ++i) {
			// construct info_bin_specific
			info_bin_specific.resize(1, 0, false);
			for (int j = 0, jmax = distance_num(i) - info_bin_res(i); j < jmax; ++j) {
				info_bin_specific.push_back(info_bin_lb(i));
			}
			for (int j = 0; j < info_bin_res(i); ++j) {
				info_bin_specific.push_back(info_bin_lb(i) + 1);
			}
			//cout << "info_bin_specific" << info_bin_specific;

			// construct distance_specific
			distance_specific.resize(1, 0, false);
			for (int j = 0, jmax = distance_num(i) - distance_res(i); j < jmax; ++j) {
				distance_specific.push_back(distance_lb(i));
			}
			for (int j = 0; j < distance_res(i); ++j) {
				distance_specific.push_back(distance_lb(i) + 1);
			}
			//cout << "distance_specific" << distance_specific;

			// construct continous_info_bin
			continous_info_bin(i).resize(1, distance_num(i), false);
			continous_info_bin(i)(0) = info_bin_continous_record(0, info_bin_specific(0) - 1);
			for (int j = 1; j < distance_num(i); ++j) {
				continous_info_bin(i)(j) = info_bin_continous_record(\
					continous_info_bin(i)(j - 1).end_ind + distance_specific(j - 1) + 1, \
					continous_info_bin(i)(j - 1).end_ind + distance_specific(j - 1) + info_bin_specific(j));
			}

			// construct permute_set
			for (int j = 0, jmax = continous_info_bin(i).size(); j < jmax; ++j) {
				for (int p = continous_info_bin(i)(j).start_ind; p <= continous_info_bin(i)(j).end_ind; ++p) {
					permute_set(i).push_back(p);
				}
			}

		}

		//cout << "continous_info_bin" << continous_info_bin;
		//cout << "permute_set" << permute_set;

		// construct shift_space
		shift_space.resize(1, set_num, false);
		for (int i = 0; i < set_num; ++i) {
			shift_space(i) = distance_lb(i) + info_bin_lb(i);
		}
		//cout << "shift_space" << shift_space;

		// add redundancy positions, remember that permuted cycle code is not cycle
		for (int i = 0; i < set_num; ++i) {
			int permute_set_i_ind = 0;
			int info_bin_res_collect_i_ind = 0;
			for (int j = 0; j < n; ++j) {
				if (permute_set_i_ind < k && j == permute_set(i)(permute_set_i_ind)) {
					++permute_set_i_ind;
				}
				else {
					permute_set(i).push_back(j);
				}
			}
		}

		//cout << "permute_set (add red)" << permute_set;

		// generate set_num systematic matrices
		G_set.resize(1, set_num, false);
		for (int i = 0; i < set_num; ++i) {
			G_set(i) = G;
			G_set(i).permute_col(permute_set(i));
			//cout << " ---------- " << endl;
			//cout << "permute_set(i)" << permute_set(i);
			//cout << "G_set(i)" << G_set(i);

			/*G_set(i).row_transformation_to_up_triangle();
			Matrix<int> permute_record = G_set(i).col_permute_to_full_rank_on_left();
			G_set(i).row_transformation_left_up_triangle_to_identity();*/

			Matrix<int> permute_record = G_set(i).GE_left_identity_4_GF2();
			permute_set(i).permute(permute_record);

			//cout << "permute_set(i) after" << permute_set(i);
			//cout << "G_set(i) after" << G_set(i);
		}
		natural = Matrix<int>(1, n, 'N');

		k_distribution_in_n.resize(1, n + k + 1, false);		// consider storing cyclic postfix
		common_pos_num.resize(1, set_num, false);
		common_pos_max_shift.resize(1, set_num, false);
	}
};

// this is for generator of extended cycle codes, where the generator matrix can shift without the last column
// e.g., BCH(127,64,21), process on its dual code since rate > 1/2 will have lower complexity
// the last column must be excluded in the Prestored matrix
class PreStored_Matrix_red_equal_dist_extend_cycle {

public:
	int n;
	int k;
	int red;
	Matrix<Matrix<int>> permute_set;
	Matrix<int> k_distribution_in_n;								// for select_G_set_continous_count
	Matrix<Matrix<info_bin_continous_record>> continous_info_bin;	// for select_G_set_continous_count
	Matrix<int> shift_space;										// for select_G_set_continous_count
	Matrix<int> common_pos_num;
	Matrix<int> common_pos_max_shift;
	Matrix<Matrix<GF2>> G_set;										// to be changed into parity-check part, the storage is different
	Matrix<int> natural;

	Basis_update BU;

	// fetch it through G_taregt, by permute the columns of G_target with columns of redundency set, then do Gaussian elemination
	// return: G_target, permute_target and the max common on G_target elements on start
	int get_MRIP_sys_G(const Matrix<int>& reliability_to_col_ind, Matrix<GF2>& G_target, Matrix<int>& permute_target) {
		// in this function we start from redundant set and find the least reliability number of the column index
		// this is of less hardware complexity, but computation cost is larger than 'get_MRIP_sys_G'

		//select_G_set(reliability_to_col_ind);
		target_ind_and_shift opt_info = select_G_set_continous_count(reliability_to_col_ind);
		int target_ind = opt_info.target_ind;
		int shift = opt_info.shift;

		//cout << "target_ind = " << target_ind << endl;
		//cout << "shift = " << shift << endl;

		//cout << "G_set(target_ind)" << G_set(target_ind);
		//cout << "permute_set(target_ind)" << permute_set(target_ind);

		//cout << "common_pos_num(target_ind) = " << common_pos_num(target_ind) << endl;
		G_target = G_set(target_ind);
		natural.col_shift_right_cir(n - 1 - shift, permute_target);			// initialize the answer
		//cout << "permute_target (shift from narual)" << permute_target;
		permute_target.push_back(n - 1);
		//cout << "permute_target (push back last position)" << permute_target;
		permute_target.permute(permute_set(target_ind));
		
		//cout << "G_target" << G_target;
		//cout << "permute_target" << permute_target;

		BU.solve(reliability_to_col_ind, G_target, permute_target);

		return common_pos_num(target_ind);
	}

	// select G that systematic position contains most position in sorted_pos
	// sorted_pos should indicates positions with ordered reliability from most to least
	target_ind_and_shift select_G_set_continous_count(const Matrix<int>& reliability_to_col_ind) {
				
		// set the distribution of MRP
		k_distribution_in_n.reset(0);		// the first position is sure to be 0
		for (int i = 0; i < k; ++i) {
			k_distribution_in_n(reliability_to_col_ind(i) + 1) = 1;
		}

		// set the post fix of k_distribution_in_n with length k
		for (int i = 1; i <= k; ++i) {
			k_distribution_in_n(n - 1 + i) = k_distribution_in_n(i);	// use n-1, ignoring the last column, be careful
		}
		//cout << "k_distribution_in_n" << k_distribution_in_n;

		// accumulate the MRP
		for (int i = 1, imax = n - 1 + k; i <= imax; ++i) {
			k_distribution_in_n(i) += k_distribution_in_n(i - 1);		// maybe need to think of parallel implementation
		}
		//cout << "k_distribution_in_n (accumulate)" << k_distribution_in_n;

		// count the common info number for each info set in  permute_set, not absolutely accurate although
		common_pos_num.reset(0);
		common_pos_max_shift.reset(0);

		//int computation_cost = 0;
		for (int i = 0, imax = continous_info_bin.size(); i < imax; i++) {
			for (int p = 0; p < shift_space(i); ++p) {
				int common_pos_num_tmp = 0;
				for (int j = 0, jmax = continous_info_bin(i).size(); j < jmax; ++j) {
					// this will save a lot complexity for selecting G 
					common_pos_num_tmp += k_distribution_in_n(p + continous_info_bin(i)(j).end_ind + 1) \
						- k_distribution_in_n(p + continous_info_bin(i)(j).start_ind);

					// anding end_ind by one for accumulated the MRP

					//computation_cost++;
				}
				bool update = common_pos_num(i) < common_pos_num_tmp;
				common_pos_num(i) = update ? common_pos_num_tmp : common_pos_num(i);
				common_pos_max_shift(i) = update ? p : common_pos_max_shift(i);
			}
		}

		// first sort the reliability_to_col_ind

		//cout << "common_pos_num" << common_pos_num;
		//cout << "total computation in select_G_set_continous_count = " << computation_cost << endl;

		// set the G_target to the G_set(index), index = arg max (common_pos_num)
		int target_ind = common_pos_num.max_ele_ind();
		int shift_val = common_pos_max_shift(target_ind);
		return target_ind_and_shift(target_ind, shift_val);
	}

	/**
	 * .initialize G_set by a generator Matrix of the code, this is offline computation, once for a code, complexity doesn't matter
	 *  but storage does matter. we only consider G such that its row number is no greater than half of its column number
	 *
	 * \param G: generator matrix or parity matrix of a code, we base on this to generate systematic matrix associated with 'MRB'
	 * \param info_partition_level: e.g., 1 means info set not partitioned, 2 means that info set is partitioned into large_bin_num parts,
	 *								3 menas that info set is partitioned into large_bin_num^2 parts, etc.
	 */
	PreStored_Matrix_red_equal_dist_extend_cycle(const Matrix<GF2>& G, bool is_generator, int partition_patterns_num) 		
		: n(G.col()), k(G.row()), red(n - k), BU(n, k) {
		
		if (is_generator == false) {
			k--;		// k should be shrink for the last position is sure in the info set
		}

		//cout << "k = " << k << endl;
		//cout << "n - k = " << n - k << endl;

		// we increase the info_partition_num from 1 until the bin_len = 1
		Matrix<int> info_bin_lb(1, k, 'v');		// info bins' length lower bound
		Matrix<int> info_bin_res(1, k, 'v');	// info bins' residual
		Matrix<int> distance_num(1, k, 'v');	// which is partition number
		info_bin_lb.push_back(k);				// firstly push back no partition info lenght, that is k
		info_bin_res.push_back(0);
		distance_num.push_back(1);
		int partition_num = 2;
		while (partition_num <= k && partition_num <= partition_patterns_num) {
			int info_len_lb_tmp = k / partition_num;
			int info_len_res_tmp = k % partition_num;
			if (info_len_lb_tmp == info_bin_lb.back()) {
				// the info length not changed, switch to info length decreasing mode
				info_bin_lb.pop_back();
				info_bin_res.pop_back();
				distance_num.pop_back();
				partition_num--;
				break;
			}
			else {
				// the info length has changed, store the info length lower bound
				info_bin_lb.push_back(info_len_lb_tmp);
				info_bin_res.push_back(info_len_res_tmp);
				distance_num.push_back(partition_num);		// partition number is totally equals to distance number
			}
			partition_num++;
		}

		// continue to the info length decreasing step
		for (int i = info_bin_lb.back() - 1; i >= 1 && partition_num <= partition_patterns_num; --i) {
			int info_len_res_tmp = k % i;		// this is magic
			int distance_num_tmp = k / i;
			if (info_len_res_tmp == distance_num_tmp) {
				continue;
			}
			info_bin_lb.push_back(i);
			info_bin_res.push_back(info_len_res_tmp);
			distance_num.push_back(distance_num_tmp);
			partition_num++;
		}
		//cout << "info_bin_lb" << info_bin_lb;
		//cout << "info_bin_res" << info_bin_res;
		//cout << "distance_num" << distance_num;

		int set_num = info_bin_lb.size();
		//cout << "set_num = " << set_num << endl;

		Matrix<int> distance_lb(1, set_num, 'v');
		Matrix<int> distance_res(1, set_num, 'v');
		for (int i = 0; i < set_num; ++i) {
			distance_lb.push_back((n - 1 - k) / distance_num(i));		// the n here is changed to n - 1
			distance_res.push_back((n - 1 - k) % distance_num(i));
		}
		//cout << "distance_lb" << distance_lb;
		//cout << "distance_res" << distance_res;

		// this is simpler, the simpler the more efficient

		continous_info_bin.resize(1, set_num, false);
		permute_set.resize(1, set_num, false);
		permute_set.reset(Matrix<int>(1, n, 'v'));

		//cout << "----------------" << endl;
		// decide each permute_set and continuous info bin
		Matrix<int> info_bin_specific(1, k);
		Matrix<int> distance_specific(1, k);
		for (int i = 0; i < set_num; ++i) {
			// construct info_bin_specific
			info_bin_specific.resize(1, 0, false);
			for (int j = 0, jmax = distance_num(i) - info_bin_res(i); j < jmax; ++j) {
				info_bin_specific.push_back(info_bin_lb(i));
			}
			for (int j = 0; j < info_bin_res(i); ++j) {
				info_bin_specific.push_back(info_bin_lb(i) + 1);
			}
			//cout << "info_bin_specific" << info_bin_specific;

			// construct distance_specific
			distance_specific.resize(1, 0, false);
			for (int j = 0, jmax = distance_num(i) - distance_res(i); j < jmax; ++j) {
				distance_specific.push_back(distance_lb(i));
			}
			for (int j = 0; j < distance_res(i); ++j) {
				distance_specific.push_back(distance_lb(i) + 1);
			}
			//cout << "distance_specific" << distance_specific;

			// construct continous_info_bin
			continous_info_bin(i).resize(1, distance_num(i), false);
			continous_info_bin(i)(0) = info_bin_continous_record(0, info_bin_specific(0) - 1);
			for (int j = 1; j < distance_num(i); ++j) {
				continous_info_bin(i)(j) = info_bin_continous_record(\
					continous_info_bin(i)(j - 1).end_ind + distance_specific(j - 1) + 1, \
					continous_info_bin(i)(j - 1).end_ind + distance_specific(j - 1) + info_bin_specific(j));
			}

			// construct permute_set
			for (int j = 0, jmax = continous_info_bin(i).size(); j < jmax; ++j) {
				for (int p = continous_info_bin(i)(j).start_ind; p <= continous_info_bin(i)(j).end_ind; ++p) {
					permute_set(i).push_back(p);
				}
			}
		}

		//cout << "continous_info_bin" << continous_info_bin;

		// construct shift_space
		shift_space.resize(1, set_num, false);
		for (int i = 0; i < set_num; ++i) {
			shift_space(i) = distance_lb(i) + info_bin_lb(i);
			//shift_space(i + set_num) = shift_space(i);	// the shift space is unchanged for permute_set index of i and i + set_num
		}
		//cout << "shift_space" << shift_space;	

		// add redundancy positions, remember that permuted cycle code is not cycle
		for (int i = 0; i < set_num; ++i) {
			int permute_set_i_ind = 0;
			int info_bin_res_collect_i_ind = 0;
			for (int j = 0; j < n; ++j) {
				if (permute_set_i_ind < k && j == permute_set(i)(permute_set_i_ind)) {
					++permute_set_i_ind;
				}
				else {
					permute_set(i).push_back(j);
				}
			}
		}
		//cout << "permute_set (add red)" << permute_set;

		// generate set_num systematic matrices
		int real_k = is_generator ? k : k + 1;
		G_set.resize(1, set_num, false);	

		for (int i = 0; i < set_num; ++i) {
			G_set(i).resize(real_k, n);
			G.get_part(0, 0, k - 1, n - 1, G_set(i));
			G_set(i).permute_col(permute_set(i));

			/*G_set(i).row_transformation_to_up_triangle();
			Matrix<int> permute_record = G_set(i).col_permute_to_full_rank_on_left();
			G_set(i).row_transformation_left_up_triangle_to_identity();*/

			Matrix<int> permute_record = G_set(i).GE_left_identity_4_GF2();
			permute_set(i).permute(permute_record);

			//cout << "permute_set(i) after" << permute_set(i);
			//cout << "G_set(i) after" << G_set(i);
		}

		if(is_generator == false){
			// add the last row to G_set(i)

			//cout << "permute_set (GE)" << permute_set;
			//cout << "G_set" << G_set;

			// for each shift, store the last column in last_row_4_shifts
			Matrix<GF2> last_row_in_G(1, n);
			for (int i = 0; i < set_num; ++i) {
				permute_set(i).pop_back();							// to adopt the size of permute_set, will add back later

				last_row_in_G = G.get_part(k, 0, k, n - 2);
				last_row_in_G.permute(permute_set(i));					// no shifting is needed for the last column
				//cout << "last_row_in_G" << last_row_in_G;

				// make sure the info set column has 0, by doing GE

				for (int p = 0; p < k; ++p) {
					if (last_row_in_G(p) == 1) {
						// add the pth row to last_row_in_G to make sure column p being 0
						last_row_in_G += G_set(i).get_part(p, 0, p, n - 2);
					}
				}

				last_row_in_G.push_back(1);

				// switch the last column to the k^th position
				G_set(i).resize(k + 1, n);

				// store the last row
				for (int j = 0; j < n; ++j) {
					G_set(i)(k, j) = last_row_in_G(j);
				}
				// switch the last column to the k^th position
				G_set(i).switch_col(k, n - 1);
				permute_set(i).push_back(n - 1);
				permute_set(i).switch_ele(k, n - 1);
			}
			k++;		// adding k back, the MRB is of length k
		}
		
		//cout << "G_set" << G_set;
		//cout << "permute_set (final)" << permute_set;

		natural = Matrix<int>(1, n - 1, 'N');
		k_distribution_in_n.resize(1, n + k + 1, false);		// consider storing cyclic postfix, actually one more space unused
		common_pos_num.resize(1, set_num, false);
		common_pos_max_shift.resize(1, set_num, false);
	}
};
