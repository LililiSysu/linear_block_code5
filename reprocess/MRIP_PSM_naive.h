/*****************************************************************//**
 * \file   MRIP_PSM_naive.h
 * \brief  class to be discarded, it is the implementation of the very first idea
 * 
 * \author 26259
 * \date   August 2023
 *********************************************************************/

#pragma once
#include"./reprocess_common.h"
using namespace std;

class PreStored_Matrix {
public:
	int n;
	int k;
	int r;
	Matrix<Matrix<int>> permute_set;
	Matrix<Matrix<int>> sorted_info_set;
	Matrix<int> k_MRP_sorted;
	Matrix<int> match_pos_num;
	Matrix<Matrix<GF2>> G_set;

	// the variable during get_MRIP_sys_G
	Matrix<int> col_ind_to_reliability;		// the smaller the reliability number, indicating that its order is more forward
	Matrix<GF2> G_target;
	Matrix<int> permute_target;
	Matrix<int> permute_target_col_ind_to_reliability;
	Matrix<int> permute_target_reliability_to_col_ind;
	Matrix<int> info_set_reliability;
	Matrix<int> redundancy_set_reliability;
	Matrix<bool> redundancy_erase_flag;
	int pass_redundancy;

	// fetch it through G_taregt, by permute the columns of G_target with columns of redundency set, then do Gaussian elemination
	void get_MRIP_sys_G(const Matrix<int>& reliability_to_col_ind) {
		// this is done for all i want for this function

		select_G_set(reliability_to_col_ind);
		//cout << "G_target" << G_target;
		//cout << "permute_target" << permute_target;

		for (int i = 0; i < n; ++i) {
			col_ind_to_reliability(reliability_to_col_ind(i)) = i;
		}
		//cout << "col_ind_to_reliability" << col_ind_to_reliability;

		for (int i = 0; i < n; ++i) {
			permute_target_col_ind_to_reliability(i) = col_ind_to_reliability(permute_target(i));
		}
		//cout << "permute_target_col_ind_to_reliability" << permute_target_col_ind_to_reliability;

		for (int i = 0; i < n; ++i) {
			permute_target_reliability_to_col_ind(permute_target_col_ind_to_reliability(i)) = i;
		}
		//cout << "permute_target_reliability_to_col_ind" << permute_target_reliability_to_col_ind;

		for (int i = 0; i < k; ++i) {
			info_set_reliability(i) = permute_target_col_ind_to_reliability(i);
		}
		info_set_reliability.sort('<');
		//cout << "info_set_reliability" << info_set_reliability;

		for (int i = k; i < n; ++i) {
			redundancy_set_reliability(i - k) = permute_target_col_ind_to_reliability(i);
		}
		redundancy_set_reliability.sort('<');
		//cout << "redundancy_set_reliability" << redundancy_set_reliability;		// has problem here, since we donot delet the element

		redundancy_erase_flag.reset(false);
		pass_redundancy = 0;
		// start from the least reliable position of info_set, permute the column to the most reliable position of redundancy set
		for (int w = k - 1; w >= 0 && info_set_reliability(w) > k; w--) {
			// iterate towards the reliable position

			// get the column number w.r.t. reliability position
			int throw_col_ind = permute_target_reliability_to_col_ind(info_set_reliability(w));
			int accept_col_ind_can = throw_col_ind;
			bool have_some_col_to_permute = false;
			for (int q = pass_redundancy; q < r && redundancy_set_reliability(q) < info_set_reliability(w); ++q) {
				if (redundancy_erase_flag(q));
				else {
					accept_col_ind_can = permute_target_reliability_to_col_ind(redundancy_set_reliability(q));
					if (G_target(throw_col_ind, accept_col_ind_can) == 1) {
						have_some_col_to_permute = true;
						redundancy_erase_flag(q) = true;
						while (pass_redundancy < r && redundancy_erase_flag(pass_redundancy) == true) {
							pass_redundancy++;
						}
						break;		// this is the column to permute
					}
				}
			}

			if (have_some_col_to_permute) {

				// update info_set_reliability
				// we erase the reliability out of the for loop, since this loop is related to info_set_reliability
				//info_set_reliability.erase(permute_target_col_ind_to_reliability(throw_col_ind));
				//info_set_reliability.insert(permute_target_col_ind_to_reliability(accept_col_ind_can));// to record max info reliability

				// update redundancy_set_reliability
				//redundancy_set_reliability.erase(permute_target_col_ind_to_reliability(accept_col_ind_can));
				//redundancy_set_reliability.insert(permute_target_col_ind_to_reliability(throw_col_ind));// do not need, this column is done

				// switch columns of throw_col_ind and accept_col_ind_can
				G_target.switch_col(throw_col_ind, accept_col_ind_can);

				// update permute_target, permute_target_col_ind_to_reliability, permute_target_reliability_to_col_ind
				permute_target.switch_ele(throw_col_ind, accept_col_ind_can);
				permute_target_col_ind_to_reliability.switch_ele(throw_col_ind, accept_col_ind_can);
				permute_target_reliability_to_col_ind(permute_target_col_ind_to_reliability(throw_col_ind)) = throw_col_ind;
				permute_target_reliability_to_col_ind(permute_target_col_ind_to_reliability(accept_col_ind_can)) = accept_col_ind_can;

				//cout << "permute_target" << permute_target;
				//cout << "G_target" << G_target;

				// do row transformation on G_target to keep identity on left
				for (int i = 0; i < k; ++i) {
					if (G_target(i, throw_col_ind) == 0 || i == throw_col_ind);
					else {
						G_target(i, throw_col_ind) = 0;
						for (int j = k; j < n; ++j) {
							G_target(i, j) += G_target(throw_col_ind, j);		// this can be done in parallel and save computation time
						}
					}
				}
				//cout << "G_target (2)" << G_target;
			}
		}

		//cout << "G_target" << G_target;

		//cout << "permute_target" << permute_target;

		//cout << "redundancy_erase_flag" << redundancy_erase_flag;

		//cout << "pass_redundancy = " << pass_redundancy << endl;

		// note that this is not the MRB yet

		// what if some MRP are out? we can continue permuting the redundancy position over the info position

		// check from the MRP of redundancy set, and try to insert it over the info set

		// some new idea should be presented here.	

		// we only consider redundancy_set with reliability <= k, since info_set with reliability > k has finished is decision
		for (int q = pass_redundancy; q < r && redundancy_set_reliability(q) < k; ++q) {

			if (redundancy_erase_flag(q)) {
				continue;
			}

			// test if can switch the iter_red to the info_set

			// count the 1 position of column w.r.t. iter_red, find the max reliability
			int accept_col_ind_can = permute_target_reliability_to_col_ind(redundancy_set_reliability(q));
			int max_possible_switch_reliability = -1;
			int throw_col_ind = -1;
			for (int i = 0; i < k; ++i) {
				if (G_target(i, accept_col_ind_can) == 1) {
					bool update = max_possible_switch_reliability < permute_target_col_ind_to_reliability(i);
					max_possible_switch_reliability = update ? permute_target_col_ind_to_reliability(i) : max_possible_switch_reliability;
					throw_col_ind = update ? i : throw_col_ind;
				}
			}

			if (max_possible_switch_reliability > redundancy_set_reliability(q)) {
				// we can switch the redundant set, to get a more reliable base

				// switch columns of throw_col_ind and accept_col_ind_can
				G_target.switch_col(throw_col_ind, accept_col_ind_can);

				//// update info_set_reliability, redundancy_set_reliability
				//redundancy_set_reliability.erase(permute_target_col_ind_to_reliability(accept_col_ind_can));
				//redundancy_set_reliability.insert(permute_target_col_ind_to_reliability(throw_col_ind));

				// update permute_target, permute_target_col_ind_to_reliability, permute_target_reliability_to_col_ind
				permute_target.switch_ele(throw_col_ind, accept_col_ind_can);
				permute_target_col_ind_to_reliability.switch_ele(throw_col_ind, accept_col_ind_can);
				permute_target_reliability_to_col_ind(permute_target_col_ind_to_reliability(throw_col_ind)) = throw_col_ind;
				permute_target_reliability_to_col_ind(permute_target_col_ind_to_reliability(accept_col_ind_can)) = accept_col_ind_can;

				//cout << "permute_target" << permute_target;
				//cout << "G_target" << G_target;

				// do row transformation on G_target to keep identity on left
				for (int i = 0; i < k; ++i) {
					if (G_target(i, throw_col_ind) == 0 || i == throw_col_ind);
					else {
						G_target(i, throw_col_ind) = 0;
						for (int j = k; j < n; ++j) {
							G_target(i, j) += G_target(throw_col_ind, j);		// this can be done in parallel and save computation time
						}
					}
				}
				//cout << "G_target (2)" << G_target;
			}
		}
	}

	// select G that systematic position contains most position in sorted_pos
	// sorted_pos should indicates positions with ordered reliability from most to least
	void select_G_set(const Matrix<int>& reliability_to_col_ind) {
		// first sort the reliability_to_col_ind

		//cout << "reliability_to_col_ind" << reliability_to_col_ind;

		reliability_to_col_ind.get_part(0, 0, 0, k - 1, k_MRP_sorted);
		k_MRP_sorted.sort('<');

		//cout << "k_MRP_sorted" << k_MRP_sorted;

		// count the match position for each permute set, that is same as sorted_info_set
		match_pos_num.reset(0);
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
					match_pos_num(i)++;
				}
			}
		}
		//cout << "match_pos_num" << match_pos_num;

		// set the G_target to the G_set(index), index = arg max (match_pos_num)
		int target_ind = match_pos_num.max_ele_ind();
		G_target = G_set(target_ind);
		permute_target = permute_set(target_ind);
	}

	// initialize G_set by a generator Matrix of the code, this is offline computation, once for a code, complexity doesn't matter
	PreStored_Matrix(const Matrix<GF2>& G, int G_partition_num) {
		n = G.col();
		k = G.row();
		r = n - k;

		// divide k into G_partition_num_part as average as possible
		int group_size_lower_bound = k / G_partition_num;
		int group_size_residual = k % G_partition_num;
		int group_size_upper_bound = group_size_lower_bound + (group_size_residual != 0);
		Matrix<int> group_size(1, G_partition_num);
		group_size.reset(group_size_lower_bound);
		for (int i = 0; i < group_size_residual; ++i) {
			group_size(i)++;	// so that the sum of group size equals k
		}
		//cout << "k = " << k << endl;
		//cout << "group_size" << group_size;
		//cout << "group_size_upper_bound = " << group_size_upper_bound << endl;
		int minority_num = my::min(group_size_residual, G_partition_num - group_size_residual);
		//cout << "minority_num = " << minority_num << endl;

		int group_place_num = (int)ceil(float(n) / group_size_upper_bound);
		//cout << "group_place_num = " << group_place_num << endl;
		int set_num = my::n_choose_k(group_place_num, G_partition_num);
		//cout << "set_num = " << set_num << endl;

		// generate group_place_num choose G_partition_num pattern
		Matrix<int> patterns = Matrix_common::generating_all_n_choose_k_pattern(group_place_num, G_partition_num);
		//cout << "patterns" << patterns;

		Matrix<int> natrual(1, n, 'N');
		Matrix<int> info_set(1, k, 'v');
		permute_set.resize(1, set_num, false);
		for (int i = 0; i < set_num; ++i) {
			info_set.resize(1, 0, false);

			// scan position of group_size to switch the size of group_size
			for (int j = my::rand_int(0, G_partition_num - 1), j1 = 0, j2 = 0; j1 != minority_num || j2 != minority_num; ++j) {
				j = j % G_partition_num;
				if (j1 != minority_num && group_size(j) == group_size_upper_bound) {
					group_size(j)--;
					j1++;
				}
				else if (j2 != minority_num && group_size(j) == group_size_lower_bound) {
					group_size(j)++;
					j2++;
				}
			}

			//group_size.permute_rand();
			//cout << "group_size" << group_size << endl;		// all done, this is very like to uniform distribution
			for (int j = 0; j < G_partition_num; ++j) {
				int group_place_pos_ind = patterns(i, j) * group_size_upper_bound;
				for (int b = 0, bmax = group_size(j); b < bmax; ++b) {
					int push_back_pos = group_place_pos_ind + b;
					if (push_back_pos < n) {
						info_set.push_back(push_back_pos);
					}
					else {
						// randomly access the unused position to push back in the info_set
						int adding_size = k - info_set.size();
						Matrix<int> natrual_unused_position = natrual.erase_cols(info_set, true);
						Matrix<int> random_position_left = natrual_unused_position.get_random_element(adding_size);

						// add the random position into info set
						for (int p = 0; p < adding_size; ++p) {
							info_set.push_back(random_position_left(p));
						}
						break;
					}
				}
			}
			//cout << "info_set" << info_set;

			// add back the info_set to form a whole set
			Matrix<int> natrual_unused_position = natrual.erase_cols(info_set, false);

			// rand permute the info_set and natrual_unused_position
			info_set.permute_rand();
			natrual_unused_position.permute_rand();
			for (int p = 0, pmax = n - k; p < pmax; ++p) {
				info_set.push_back(natrual_unused_position(p));
			}

			// store into info_set_matrix
			permute_set(i) = info_set;
		}

		//cout << "permute_set" << permute_set;

		// generate set_num systematic generator matrix
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


			//cout << "info_set_matrix(i)" << permute_set(i);
			//cout << "G_set(i)" << G_set(i);
		}

		sorted_info_set.resize(1, set_num, false);
		for (int i = 0; i < set_num; ++i) {

			// store the first k value of permute_set into sorted_info_set, where the latter is sorted index for choosing best G_set		
			sorted_info_set(i) = permute_set(i).get_part(0, 0, 0, k - 1);
			sorted_info_set(i).sort('<');
		}

		k_MRP_sorted.resize(1, k, false);
		match_pos_num.resize(1, set_num, false);
		G_target.resize(k, n, false);
		permute_target.resize(1, n, false);
		col_ind_to_reliability.resize(1, n, false);
		permute_target_col_ind_to_reliability.resize(1, n, false);
		permute_target_reliability_to_col_ind.resize(1, n, false);
		info_set_reliability.resize(1, k, false);
		redundancy_set_reliability.resize(1, r, false);
		redundancy_erase_flag.resize(1, r, false);
		pass_redundancy = 0;
	}
};

class PreStored_Matrix_red {

public:
	int n;
	int k;
	int r;
	Matrix<Matrix<int>> permute_set;
	Matrix<Matrix<int>> sorted_info_set;
	Matrix<int> k_MRP_sorted;
	Matrix<int> match_pos_num;
	Matrix<Matrix<GF2>> G_set;

	// the variable during get_MRIP_sys_G
	Matrix<int> col_ind_to_reliability;		// the smaller the reliability number, indicating that its order is more forward
	Matrix<GF2> G_target;
	Matrix<int> permute_target;
	Matrix<int> permute_target_col_ind_to_reliability;
	Matrix<int> permute_target_reliability_to_col_ind;
	Matrix<int> redundancy_set_reliability;
	Heap_max<int> info_set_reliability_heap;

	// fetch it through G_taregt, by permute the columns of G_target with columns of redundency set, then do Gaussian elemination
	void get_MRIP_sys_G(const Matrix<int>& reliability_to_col_ind) {
		// in this function we start from redundant set and find the least reliability number of the column index
		// this is of less hardware complexity, but computation cost is larger than 'get_MRIP_sys_G'

		select_G_set(reliability_to_col_ind);
		//cout << "G_target" << G_target;
		//cout << "permute_target" << permute_target;

		for (int i = 0; i < n; ++i) {
			col_ind_to_reliability(reliability_to_col_ind(i)) = i;
		}
		//cout << "col_ind_to_reliability" << col_ind_to_reliability;

		for (int i = 0; i < n; ++i) {
			permute_target_col_ind_to_reliability(i) = col_ind_to_reliability(permute_target(i));
		}
		//cout << "permute_target_col_ind_to_reliability" << permute_target_col_ind_to_reliability;

		for (int i = 0; i < n; ++i) {
			permute_target_reliability_to_col_ind(permute_target_col_ind_to_reliability(i)) = i;
		}
		//cout << "permute_target_reliability_to_col_ind" << permute_target_reliability_to_col_ind;


		// how to mantain the max reliability number of info set? may be using priority queue, or heap of myself
		for (int i = 0; i < k; ++i) {
			info_set_reliability_heap(i) = permute_target_col_ind_to_reliability(i);
		}
		//cout << "info_set_reliability_heap (before build)" << info_set_reliability_heap;
		info_set_reliability_heap.build();
		//cout << "info_set_reliability_heap (after build)" << info_set_reliability_heap;

		for (int i = k; i < n; ++i) {
			redundancy_set_reliability(i - k) = permute_target_col_ind_to_reliability(i);
		}
		redundancy_set_reliability.sort('<');
		//cout << "redundancy_set_reliability" << redundancy_set_reliability;		// has problem here, since we donot delet the element

		// note that this is not the MRB yet

		// what if some MRP are out? we can continue permuting the redundancy position over the info position

		// check from the MRP of redundancy set, and try to insert it over the info set

		// some new idea should be presented here.	

		// we only consider redundancy_set with reliability <= k, since info_set with reliability > k has finished is decision
		for (int q = 0; q < r && redundancy_set_reliability(q) < info_set_reliability_heap.top(); ++q) {

			// test if can switch the iter_red to the info_set

			// count the 1 position of column w.r.t. iter_red, find the max reliability
			int accept_col_ind_can = permute_target_reliability_to_col_ind(redundancy_set_reliability(q));
			int max_possible_switch_reliability = -1;
			int throw_col_ind = -1;
			for (int i = 0; i < k; ++i) {
				if (G_target(i, accept_col_ind_can) == 1) {	// this is for software implementation, some optimization to be done for hardware
					bool update = max_possible_switch_reliability < permute_target_col_ind_to_reliability(i);
					max_possible_switch_reliability = update ? permute_target_col_ind_to_reliability(i) : max_possible_switch_reliability;
					throw_col_ind = update ? i : throw_col_ind;
				}
			}

			if (max_possible_switch_reliability > redundancy_set_reliability(q)) {
				// we can switch the redundant set, to get a more reliable base

				// switch columns of throw_col_ind and accept_col_ind_can
				G_target.switch_col(throw_col_ind, accept_col_ind_can);

				// mantian the info_set_reliability_heap

				// find the throw_col_ind through iteration over info_set_reliability_heap, we must find it, by increasing the index
				int info_set_reliability_heap_throw_ind;
				for (info_set_reliability_heap_throw_ind = 0; info_set_reliability_heap_throw_ind < k && \
					info_set_reliability_heap(info_set_reliability_heap_throw_ind) != max_possible_switch_reliability; \
					++info_set_reliability_heap_throw_ind);
				// since we replace the column with as most reliability number as possible, this for loop will end shortly

				info_set_reliability_heap.replace(info_set_reliability_heap_throw_ind, \
					permute_target_col_ind_to_reliability(accept_col_ind_can));

				// update permute_target, permute_target_col_ind_to_reliability, permute_target_reliability_to_col_ind
				permute_target.switch_ele(throw_col_ind, accept_col_ind_can);
				permute_target_col_ind_to_reliability.switch_ele(throw_col_ind, accept_col_ind_can);
				permute_target_reliability_to_col_ind(permute_target_col_ind_to_reliability(throw_col_ind)) = throw_col_ind;
				permute_target_reliability_to_col_ind(permute_target_col_ind_to_reliability(accept_col_ind_can)) = accept_col_ind_can;

				//cout << "permute_target" << permute_target;
				//cout << "G_target" << G_target;

				// do row transformation on G_target to keep identity on left
				for (int i = 0; i < k; ++i) {
					if (G_target(i, throw_col_ind) == 1 && i != throw_col_ind) {
						G_target(i, throw_col_ind) = 0;
						for (int j = k; j < n; ++j) {
							G_target(i, j) += G_target(throw_col_ind, j);		// this can be done in parallel and save computation time
						}
					}
				}
				//cout << "G_target (2)" << G_target;
			}
		}
	}

	// select G that systematic position contains most position in sorted_pos
	// sorted_pos should indicates positions with ordered reliability from most to least
	void select_G_set(const Matrix<int>& reliability_to_col_ind) {
		// first sort the reliability_to_col_ind

		//cout << "reliability_to_col_ind" << reliability_to_col_ind;

		reliability_to_col_ind.get_part(0, 0, 0, k - 1, k_MRP_sorted);
		k_MRP_sorted.sort('<');

		//cout << "k_MRP_sorted" << k_MRP_sorted;

		// count the match position for each permute set, that is same as sorted_info_set
		match_pos_num.reset(0);
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
					match_pos_num(i)++;
				}
			}
		}
		//cout << "match_pos_num" << match_pos_num;

		// set the G_target to the G_set(index), index = arg max (match_pos_num)
		int target_ind = match_pos_num.max_ele_ind();
		G_target = G_set(target_ind);
		permute_target = permute_set(target_ind);
	}

	// initialize G_set by a generator Matrix of the code, this is offline computation, once for a code, complexity doesn't matter
	PreStored_Matrix_red(const Matrix<GF2>& G, int G_partition_num) {
		n = G.col();
		k = G.row();
		r = n - k;

		// divide k into G_partition_num_part as average as possible
		int group_size_lower_bound = k / G_partition_num;
		int group_size_residual = k % G_partition_num;
		int group_size_upper_bound = group_size_lower_bound + (group_size_residual != 0);
		Matrix<int> group_size(1, G_partition_num);
		group_size.reset(group_size_lower_bound);
		for (int i = 0; i < group_size_residual; ++i) {
			group_size(i)++;	// so that the sum of group size equals k
		}
		//cout << "k = " << k << endl;
		//cout << "group_size" << group_size;
		//cout << "group_size_upper_bound = " << group_size_upper_bound << endl;
		int minority_num = my::min(group_size_residual, G_partition_num - group_size_residual);
		//cout << "minority_num = " << minority_num << endl;

		int group_place_num = (int)ceil(float(n) / group_size_upper_bound);
		//cout << "group_place_num = " << group_place_num << endl;
		int set_num = my::n_choose_k(group_place_num, G_partition_num);
		//cout << "set_num = " << set_num << endl;

		// generate group_place_num choose G_partition_num pattern
		Matrix<int> patterns = Matrix_common::generating_all_n_choose_k_pattern(group_place_num, G_partition_num);
		//cout << "patterns" << patterns;

		Matrix<int> natrual(1, n, 'N');
		Matrix<int> info_set(1, k, 'v');
		permute_set.resize(1, set_num, false);
		for (int i = 0; i < set_num; ++i) {
			info_set.resize(1, 0, false);

			// scan position of group_size to switch the size of group_size
			for (int j = my::rand_int(0, G_partition_num - 1), j1 = 0, j2 = 0; j1 != minority_num || j2 != minority_num; ++j) {
				j = j % G_partition_num;
				if (j1 != minority_num && group_size(j) == group_size_upper_bound) {
					group_size(j)--;
					j1++;
				}
				else if (j2 != minority_num && group_size(j) == group_size_lower_bound) {
					group_size(j)++;
					j2++;
				}
			}

			//group_size.permute_rand();
			//cout << "group_size" << group_size << endl;		// all done, this is very like to uniform distribution
			for (int j = 0; j < G_partition_num; ++j) {
				int group_place_pos_ind = patterns(i, j) * group_size_upper_bound;
				for (int b = 0, bmax = group_size(j); b < bmax; ++b) {
					int push_back_pos = group_place_pos_ind + b;
					if (push_back_pos < n) {
						info_set.push_back(push_back_pos);
					}
					else {
						// randomly access the unused position to push back in the info_set
						int adding_size = k - info_set.size();
						Matrix<int> natrual_unused_position = natrual.erase_cols(info_set, true);
						Matrix<int> random_position_left = natrual_unused_position.get_random_element(adding_size);

						// add the random position into info set
						for (int p = 0; p < adding_size; ++p) {
							info_set.push_back(random_position_left(p));
						}
						break;
					}
				}
			}
			//cout << "info_set" << info_set;

			// add back the info_set to form a whole set
			Matrix<int> natrual_unused_position = natrual.erase_cols(info_set, false);

			// rand permute the info_set and natrual_unused_position
			info_set.permute_rand();
			natrual_unused_position.permute_rand();
			for (int p = 0, pmax = n - k; p < pmax; ++p) {
				info_set.push_back(natrual_unused_position(p));
			}

			// store into info_set_matrix
			permute_set(i) = info_set;
		}

		//cout << "permute_set" << permute_set;

		// generate set_num systematic generator matrix
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

			//cout << "info_set_matrix(i)" << permute_set(i);
			//cout << "G_set(i)" << G_set(i);
		}

		sorted_info_set.resize(1, set_num, false);
		for (int i = 0; i < set_num; ++i) {

			// store the first k value of permute_set into sorted_info_set, where the latter is sorted index for choosing best G_set		
			sorted_info_set(i) = permute_set(i).get_part(0, 0, 0, k - 1);
			sorted_info_set(i).sort('<');
		}

		k_MRP_sorted.resize(1, k, false);
		match_pos_num.resize(1, set_num, false);
		G_target.resize(k, n, false);
		permute_target.resize(1, n, false);
		col_ind_to_reliability.resize(1, n, false);
		permute_target_col_ind_to_reliability.resize(1, n, false);
		permute_target_reliability_to_col_ind.resize(1, n, false);
		redundancy_set_reliability.resize(1, r, false);
		info_set_reliability_heap.resize(1, k, false);
	}
};

// this is for general code
class PreStored_Matrix_red_recursive_partition {

public:
	int n;
	int k;
	int red;
	Matrix<Matrix<int>> permute_set;
	Matrix<Matrix<int>> sorted_info_set;				// for select_G_set
	Matrix<int> k_MRP_sorted;							// for select_G_set
	Matrix<int> k_distribution_in_n;					// for select_G_set_continous_count
	Matrix<Matrix<info_bin_continous_record>> continous_info_bin;	// for select_G_set_continous_count
	Matrix<int> match_pos_num;
	Matrix<Matrix<GF2>> G_set;

	// the variable during get_MRIP_sys_G
	Matrix<int> col_ind_to_reliability;					// the smaller the reliability number, indicating that its order is more forward
	Matrix<GF2> G_target;
	Matrix<int> permute_target;
	Matrix<int> permute_target_col_ind_to_reliability;
	Matrix<int> permute_target_reliability_to_col_ind;
	Matrix<int> redundancy_set_reliability;
	Heap_max<int> info_set_reliability_heap;

	// fetch it through G_taregt, by permute the columns of G_target with columns of redundency set, then do Gaussian elemination
	void get_MRIP_sys_G(const Matrix<int>& reliability_to_col_ind) {
		// in this function we start from redundant set and find the least reliability number of the column index
		// this is of less hardware complexity, but computation cost is larger than 'get_MRIP_sys_G'

		//select_G_set(reliability_to_col_ind);
		select_G_set_continous_count(reliability_to_col_ind);
		//cout << "G_target" << G_target;
		//cout << "permute_target" << permute_target;

		for (int i = 0; i < n; ++i) {
			col_ind_to_reliability(reliability_to_col_ind(i)) = i;
		}
		//cout << "col_ind_to_reliability" << col_ind_to_reliability;

		for (int i = 0; i < n; ++i) {
			permute_target_col_ind_to_reliability(i) = col_ind_to_reliability(permute_target(i));
		}
		//cout << "permute_target_col_ind_to_reliability" << permute_target_col_ind_to_reliability;

		for (int i = 0; i < n; ++i) {
			permute_target_reliability_to_col_ind(permute_target_col_ind_to_reliability(i)) = i;
		}
		//cout << "permute_target_reliability_to_col_ind" << permute_target_reliability_to_col_ind;


		// how to mantain the max reliability number of info set? may be using priority queue, or heap of myself
		for (int i = 0; i < k; ++i) {
			info_set_reliability_heap(i) = permute_target_col_ind_to_reliability(i);
		}
		//cout << "info_set_reliability_heap (before build)" << info_set_reliability_heap;
		info_set_reliability_heap.build();
		//cout << "info_set_reliability_heap (after build)" << info_set_reliability_heap;

		for (int i = k; i < n; ++i) {
			redundancy_set_reliability(i - k) = permute_target_col_ind_to_reliability(i);
		}
		redundancy_set_reliability.sort('<');
		//cout << "redundancy_set_reliability" << redundancy_set_reliability;		// has problem here, since we donot delet the element

		// note that this is not the MRB yet

		// what if some MRP are out? we can continue permuting the redundancy position over the info position

		// check from the MRP of redundancy set, and try to insert it over the info set

		// some new idea should be presented here.	

		// we only consider redundancy_set with reliability <= k, since info_set with reliability > k has finished is decision
		for (int q = 0; q < red && redundancy_set_reliability(q) < info_set_reliability_heap.top(); ++q) {

			// test if can switch the iter_red to the info_set

			// count the 1 position of column w.r.t. iter_red, find the max reliability
			int accept_col_ind_can = permute_target_reliability_to_col_ind(redundancy_set_reliability(q));
			int max_possible_switch_reliability = -1;
			int throw_col_ind = -1;
			for (int i = 0; i < k; ++i) {
				if (G_target(i, accept_col_ind_can) == 1) {	// this is for software implementation, some optimization to be done for hardware
					bool update = max_possible_switch_reliability < permute_target_col_ind_to_reliability(i);
					max_possible_switch_reliability = update ? permute_target_col_ind_to_reliability(i) : max_possible_switch_reliability;
					throw_col_ind = update ? i : throw_col_ind;
				}
			}

			if (max_possible_switch_reliability > redundancy_set_reliability(q)) {
				// we can switch the redundant set, to get a more reliable base

				// switch columns of throw_col_ind and accept_col_ind_can
				G_target.switch_col(throw_col_ind, accept_col_ind_can);

				// mantian the info_set_reliability_heap

				// find the throw_col_ind through iteration over info_set_reliability_heap, we must find it, by increasing the index
				int info_set_reliability_heap_throw_ind;
				for (info_set_reliability_heap_throw_ind = 0; info_set_reliability_heap_throw_ind < k && \
					info_set_reliability_heap(info_set_reliability_heap_throw_ind) != max_possible_switch_reliability; \
					++info_set_reliability_heap_throw_ind);
				// since we replace the column with as most reliability number as possible, this for loop will end shortly

				info_set_reliability_heap.replace(info_set_reliability_heap_throw_ind, \
					permute_target_col_ind_to_reliability(accept_col_ind_can));

				// update permute_target, permute_target_col_ind_to_reliability, permute_target_reliability_to_col_ind
				permute_target.switch_ele(throw_col_ind, accept_col_ind_can);
				permute_target_col_ind_to_reliability.switch_ele(throw_col_ind, accept_col_ind_can);
				permute_target_reliability_to_col_ind(permute_target_col_ind_to_reliability(throw_col_ind)) = throw_col_ind;
				permute_target_reliability_to_col_ind(permute_target_col_ind_to_reliability(accept_col_ind_can)) = accept_col_ind_can;

				//cout << "permute_target" << permute_target;
				//cout << "G_target" << G_target;

				// do row transformation on G_target to keep identity on left
				for (int i = 0; i < k; ++i) {
					if (G_target(i, throw_col_ind) == 1 && i != throw_col_ind) {
						G_target(i, throw_col_ind) = 0;
						for (int j = k; j < n; ++j) {
							G_target(i, j) += G_target(throw_col_ind, j);		// this can be done in parallel and save computation time
						}
					}
				}
				//cout << "G_target (2)" << G_target;
			}
		}
	}

	// select G that systematic position contains most position in sorted_pos
	// sorted_pos should indicates positions with ordered reliability from most to least
	void select_G_set_accurate(const Matrix<int>& reliability_to_col_ind) {
		// first sort the reliability_to_col_ind

		//cout << "reliability_to_col_ind" << reliability_to_col_ind;

		reliability_to_col_ind.get_part(0, 0, 0, k - 1, k_MRP_sorted);
		k_MRP_sorted.sort('<');

		//cout << "k_MRP_sorted" << k_MRP_sorted;

		// count the match position for each permute set, that is same as sorted_info_set
		match_pos_num.reset(0);
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
					match_pos_num(i)++;
				}
			}
		}
		//cout << "match_pos_num" << match_pos_num;

		// set the G_target to the G_set(index), index = arg max (match_pos_num)
		int target_ind = match_pos_num.max_ele_ind();
		G_target = G_set(target_ind);
		permute_target = permute_set(target_ind);
	}

	// select G that systematic position contains most position in sorted_pos
	// sorted_pos should indicates positions with ordered reliability from most to least
	void select_G_set_continous_count(const Matrix<int>& reliability_to_col_ind) {

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
		match_pos_num.reset(0);
		for (int i = 0, imax = continous_info_bin.size(); i < imax; i++) {
			for (int j = 0, jmax = continous_info_bin(i).size(); j < jmax; ++j) {
				// this will save a lot complexity for selecting G 
				match_pos_num(i) += \
					k_distribution_in_n(continous_info_bin(i)(j).end_ind + 1) - k_distribution_in_n(continous_info_bin(i)(j).start_ind);
				// anding end_ind by one for accumulated the MRP
			}
		}

		// first sort the reliability_to_col_ind

		//cout << "match_pos_num" << match_pos_num;

		// set the G_target to the G_set(index), index = arg max (match_pos_num)
		int target_ind = match_pos_num.max_ele_ind();
		G_target = G_set(target_ind);
		permute_target = permute_set(target_ind);
	}

	/**
	 * .initialize G_set by a generator Matrix of the code, this is offline computation, once for a code, complexity doesn't matter
	 *  but storage does matter. we only consider G such that its row number is no greater than half of its column number
	 *
	 * \param G: generator matrix or parity matrix of a code, we base on this to generate systematic matrix associated with 'MRB'
	 * \param info_partition_level: e.g., 1 means info set not partitioned, 2 means that info set is partitioned into large_bin_num parts,
	 *								3 menas that info set is partitioned into large_bin_num^2 parts, etc.
	 */
	PreStored_Matrix_red_recursive_partition(const Matrix<GF2>& G, int info_partition_level) {
		n = G.col();
		k = G.row();
		red = n - k;

		if (k > red) {
			cout << "k = " << k << endl;
			cout << "n - k = " << n - k << endl;
			cout << "k > n - k, not allowed in [PreStored_Matrix_red_recursive_partition]" << endl;
			return;
		}

		int div_num = (int)floor((double)n / k);
		//cout << "div_num = " << div_num << endl;			// this number should greater than 1

		int set_num = 0;
		for (int i = 0; i < info_partition_level; ++i) {
			set_num += (int)pow(div_num, pow(div_num, i));	// this is exponential of exponential, be careful
		}
		//cout << "set_num = " << set_num << endl;

		permute_set.resize(1, set_num, false);
		permute_set.reset(Matrix<int>(1, n, 'v'));
		//cout << "permute_set" << permute_set;

		int permute_set_start = 0;

		Matrix<int> nbin(1, (int)round((1 - pow(div_num, info_partition_level + 1)) / (1 - div_num)), 'v');	// generate bin len recursively
		Matrix<int> kbin(1, (int)round((1 - pow(div_num, info_partition_level)) / (1 - div_num)), 'v');		// sum of geometric progression
		int nbin_start = 0;
		int kbin_start = 0;
		nbin.push_back(n);
		kbin.push_back(k);
		continous_info_bin.resize(1, set_num, false);

		for (int i = 0; i < info_partition_level; ++i) {
			int qs = (int)pow(div_num, i);
			for (int j = 0; j < qs; ++j) {
				int len = nbin(nbin_start);
				nbin_start++;				// this is identitcal to pop() in queue

				// partition len into div_num parts, with small bin at left
				int lb = len / div_num;		// lower bound of divided bin
				int ub = lb + 1;			// upper bound of divided bin
				int res = len % div_num;
				for (int w = 0, wmax = div_num - res; w < wmax; ++w) {
					nbin.push_back(lb);		// push lb at frist
				}
				for (int w = 0; w < res; ++w) {
					nbin.push_back(ub);
				}
			}
			//cout << "nbin" << nbin;

			//cout << "kbin" << kbin;

			int kbin_size = kbin.size() - kbin_start;		// the father node number, div_number is the branch number
			//int nbin_size = kbin_size * div_num;

			Matrix<int> my_scan = Matrix_common::n_ary_incremental_scan(div_num, kbin_size);
			//cout << "my_scan" << my_scan;

			Matrix<int> permute_set_bias(kbin_size, div_num);
			permute_set_bias(0) = 0;
			for (int g = 0, gmax = kbin_size * div_num - 1; g < gmax; ++g) {
				permute_set_bias(g + 1) = permute_set_bias(g) + nbin(nbin_start + g);
			}
			//cout << "permute_set_bias" << permute_set_bias;

			// diliver all number in kbin to nbin, forming info_set
			for (int b = 0, bmax = my_scan.row(); b < bmax; ++b) {
				continous_info_bin(permute_set_start).resize(1, kbin_size, false);

				// iterate through all the permute_set to deliver one bin in kbin
				for (int q = 0; q < kbin_size; ++q) {

					int codeword_bin = nbin(nbin_start + q * div_num + my_scan(b, q));
					int k_bin = kbin(kbin_start + q);
					int info_pos_bias = permute_set_bias(q, my_scan(b, q));

					// set continous_info_bin, in usual codeword_bin >= k_bin
					continous_info_bin(permute_set_start)(q) = \
						info_bin_continous_record(info_pos_bias, info_pos_bias + my::max(codeword_bin, k_bin) - 1);

					Matrix<int> random_deliver(1, k_bin, 'N');
					if (codeword_bin > k_bin) {
						random_deliver = Matrix_common::n_randomly_choose_k(codeword_bin, k_bin);
						random_deliver.sort('<');
					}

					for (int z = 0, zmax = random_deliver.size(); z < zmax; ++z) {
						permute_set(permute_set_start).push_back(random_deliver(z) + info_pos_bias);
					}
				}
				permute_set_start++;
			}

			if (i + 1 < info_partition_level) {
				// update k_bin
				for (int qs = 0; qs < kbin_size; ++qs) {
					int len = kbin(kbin_start);
					kbin_start++;				// this is identitcal to pop() in queue

					// partition len into div_num parts, with small bin at left
					int lb = len / div_num;		// lower bound of divided bin
					int ub = lb + 1;			// upper bound of divided bin
					int res = len % div_num;
					for (int w = 0, wmax = div_num - res; w < wmax; ++w) {
						kbin.push_back(lb);		// push lb at frist
					}
					for (int w = 0; w < res; ++w) {
						kbin.push_back(ub);
					}
				}
			}
		}

		//cout << "continous_info_bin" << continous_info_bin;

		//cout << "permute_set" << permute_set;

		// to do below: adding redundancy position, GE and store

		// add redundancy positions
		for (int i = 0; i < set_num; ++i) {
			int permute_set_i_ind = 0;
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

		// generate set_num systematic generator matrix
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
		match_pos_num.resize(1, set_num, false);
		G_target.resize(k, n, false);
		permute_target.resize(1, n, false);
		col_ind_to_reliability.resize(1, n, false);
		permute_target_col_ind_to_reliability.resize(1, n, false);
		permute_target_reliability_to_col_ind.resize(1, n, false);
		redundancy_set_reliability.resize(1, red, false);
		info_set_reliability_heap.resize(1, k, false);
	}
};

// this is a loose bin plan, partition n into loose bins greater or equal to k, then partition these bins into info_partition_num parts
class PreStored_Matrix_red_loose_bin {

public:
	int n;
	int k;
	int red;
	Matrix<Matrix<int>> permute_set;
	Matrix<Matrix<int>> sorted_info_set;				// for select_G_set, to be removed
	Matrix<int> k_MRP_sorted;							// for select_G_set
	Matrix<int> k_distribution_in_n;					// for select_G_set_continous_count
	Matrix<Matrix<info_bin_continous_record>> continous_info_bin;	// for select_G_set_continous_count
	Matrix<int> match_pos_num;
	Matrix<Matrix<GF2>> G_set;							// to be changed into parity-check part

	// the variable during get_MRIP_sys_G
	Matrix<int> col_ind_to_reliability;					// the smaller the reliability number, indicating that its order is more forward
	//Matrix<GF2> G_target;
	//Matrix<int> permute_target;
	Matrix<int> permute_target_col_ind_to_reliability;
	Matrix<int> permute_target_reliability_to_col_ind;
	Matrix<int> redundancy_set_reliability;
	Heap_max<int> info_set_reliability_heap;

	// fetch it through G_taregt, by permute the columns of G_target with columns of redundency set, then do Gaussian elemination
	// return: G_target, permute_target and the max common on G_target elements on start
	int get_MRIP_sys_G(const Matrix<int>& reliability_to_col_ind, Matrix<GF2>& G_target, Matrix<int>& permute_target) {
		// in this function we start from redundant set and find the least reliability number of the column index
		// this is of less hardware complexity, but computation cost is larger than 'get_MRIP_sys_G'

		//select_G_set(reliability_to_col_ind);
		int target_ind = select_G_set_continous_count(reliability_to_col_ind);

		//cout << "match_pos_num(target_ind) = " << match_pos_num(target_ind) << endl;
		G_target = G_set(target_ind);
		permute_target = permute_set(target_ind);			// initialize the answer

		//cout << "G_target" << G_target;
		//cout << "permute_target" << permute_target;

		for (int i = 0; i < n; ++i) {
			col_ind_to_reliability(reliability_to_col_ind(i)) = i;
		}
		//cout << "col_ind_to_reliability" << col_ind_to_reliability;

		for (int i = 0; i < n; ++i) {
			permute_target_col_ind_to_reliability(i) = col_ind_to_reliability(permute_target(i));
		}
		//cout << "permute_target_col_ind_to_reliability" << permute_target_col_ind_to_reliability;

		for (int i = 0; i < n; ++i) {
			permute_target_reliability_to_col_ind(permute_target_col_ind_to_reliability(i)) = i;
		}
		//cout << "permute_target_reliability_to_col_ind" << permute_target_reliability_to_col_ind;


		// how to mantain the max reliability number of info set? may be using priority queue, or heap of myself
		for (int i = 0; i < k; ++i) {
			info_set_reliability_heap(i) = permute_target_col_ind_to_reliability(i);
		}
		//cout << "info_set_reliability_heap (before build)" << info_set_reliability_heap;
		info_set_reliability_heap.build();
		//cout << "info_set_reliability_heap (after build)" << info_set_reliability_heap;

		for (int i = k; i < n; ++i) {
			redundancy_set_reliability(i - k) = permute_target_col_ind_to_reliability(i);
		}
		redundancy_set_reliability.sort('<');
		//cout << "redundancy_set_reliability" << redundancy_set_reliability;		// has problem here, since we donot delet the element

		// note that this is not the MRB yet

		// what if some MRP are out? we can continue permuting the redundancy position over the info position

		// check from the MRP of redundancy set, and try to insert it over the info set

		// some new idea should be presented here.	

		// we only consider redundancy_set with reliability <= k, since info_set with reliability > k has finished is decision
		for (int q = 0; q < red && redundancy_set_reliability(q) < info_set_reliability_heap.top(); ++q) {

			// test if can switch the iter_red to the info_set

			// count the 1 position of column w.r.t. iter_red, find the max reliability
			int accept_col_ind_can = permute_target_reliability_to_col_ind(redundancy_set_reliability(q));
			int max_possible_switch_reliability = -1;
			int throw_col_ind = -1;
			for (int i = 0; i < k; ++i) {
				if (G_target(i, accept_col_ind_can) == 1) {	// this is for software implementation, some optimization to be done for hardware
					bool update = max_possible_switch_reliability < permute_target_col_ind_to_reliability(i);
					max_possible_switch_reliability = update ? permute_target_col_ind_to_reliability(i) : max_possible_switch_reliability;
					throw_col_ind = update ? i : throw_col_ind;
				}
			}

			if (max_possible_switch_reliability > redundancy_set_reliability(q)) {
				// we can switch the redundant set, to get a more reliable base

				// switch columns of throw_col_ind and accept_col_ind_can
				G_target.switch_col(throw_col_ind, accept_col_ind_can);

				// mantian the info_set_reliability_heap

				// find the throw_col_ind through iteration over info_set_reliability_heap, we must find it, by increasing the index
				int info_set_reliability_heap_throw_ind;
				for (info_set_reliability_heap_throw_ind = 0; info_set_reliability_heap_throw_ind < k && \
					info_set_reliability_heap(info_set_reliability_heap_throw_ind) != max_possible_switch_reliability; \
					++info_set_reliability_heap_throw_ind);
				// since we replace the column with as most reliability number as possible, this for loop will end shortly

				info_set_reliability_heap.replace(info_set_reliability_heap_throw_ind, \
					permute_target_col_ind_to_reliability(accept_col_ind_can));

				// update permute_target, permute_target_col_ind_to_reliability, permute_target_reliability_to_col_ind
				permute_target.switch_ele(throw_col_ind, accept_col_ind_can);
				permute_target_col_ind_to_reliability.switch_ele(throw_col_ind, accept_col_ind_can);
				permute_target_reliability_to_col_ind(permute_target_col_ind_to_reliability(throw_col_ind)) = throw_col_ind;
				permute_target_reliability_to_col_ind(permute_target_col_ind_to_reliability(accept_col_ind_can)) = accept_col_ind_can;

				//cout << "permute_target" << permute_target;
				//cout << "G_target" << G_target;

				// do row transformation on G_target to keep identity on left
				for (int i = 0; i < k; ++i) {
					if (G_target(i, throw_col_ind) == 1 && i != throw_col_ind) {
						G_target(i, throw_col_ind) = 0;
						for (int j = k; j < n; ++j) {
							G_target(i, j) += G_target(throw_col_ind, j);		// this can be done in parallel and save computation time
						}
					}
				}
				//cout << "G_target (2)" << G_target;
			}
		}

		return match_pos_num(target_ind);
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
		match_pos_num.reset(0);
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
					match_pos_num(i)++;
				}
			}
		}
		//cout << "match_pos_num" << match_pos_num;

		// set the G_target to the G_set(index), index = arg max (match_pos_num)
		int target_ind = match_pos_num.max_ele_ind();
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
		match_pos_num.reset(0);
		for (int i = 0, imax = continous_info_bin.size(); i < imax; i++) {
			for (int j = 0, jmax = continous_info_bin(i).size(); j < jmax; ++j) {
				// this will save a lot complexity for selecting G 
				match_pos_num(i) += \
					k_distribution_in_n(continous_info_bin(i)(j).end_ind + 1) - k_distribution_in_n(continous_info_bin(i)(j).start_ind);
				// anding end_ind by one for accumulated the MRP
			}
		}

		// first sort the reliability_to_col_ind

		//cout << "match_pos_num" << match_pos_num;

		// set the G_target to the G_set(index), index = arg max (match_pos_num)
		int target_ind = match_pos_num.max_ele_ind();
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
	PreStored_Matrix_red_loose_bin(const Matrix<GF2>& G, int info_partition_num) {
		n = G.col();
		k = G.row();
		red = n - k;

		if (k > red) {
			cout << "k = " << k << endl;
			cout << "n - k = " << n - k << endl;
			cout << "Warning: k > n - k,  [PreStored_Matrix_red_loose_bin] is trivial" << endl;
		}

		cout << "k = " << k << endl;
		cout << "n - k = " << n - k << endl;

		// generate the loose_bin_len 
		int loose_bin_num = n / k;
		Matrix<int> loose_bin_len(1, loose_bin_num, 'v');
		int loose_bin_len_lb = n / loose_bin_num;
		int loose_bin_res = n % loose_bin_num;
		for (int i = 0, imax = loose_bin_num - loose_bin_res; i < imax; ++i) {
			loose_bin_len.push_back(loose_bin_len_lb);
		}
		for (int i = 0; i < loose_bin_res; ++i) {
			loose_bin_len.push_back(loose_bin_len_lb + 1);
		}

		// generate the info_part_size 
		Matrix<int> info_part_size(1, info_partition_num, 'v');
		int info_part_len_lb = k / info_partition_num;
		int info_part_res = k % info_partition_num;
		for (int i = 0, imax = info_partition_num - info_part_res; i < imax; ++i) {
			info_part_size.push_back(info_part_len_lb);
		}
		for (int i = 0; i < info_part_res; ++i) {
			info_part_size.push_back(info_part_len_lb + 1);
		}
		cout << "info_part_size" << info_part_size;

		// generate the codeword_bin_num 
		int codeword_bin_num = info_partition_num * loose_bin_num;
		cout << "codeword_bin_num = " << codeword_bin_num << endl;
		Matrix<int> bin_size(1, codeword_bin_num, 'v');
		for (int i = 0; i < loose_bin_num; ++i) {
			int bin_len_lb = loose_bin_len(i) / info_partition_num;
			int bin_res = loose_bin_len(i) % info_partition_num;
			for (int i = 0, imax = info_partition_num - bin_res; i < imax; ++i) {
				bin_size.push_back(bin_len_lb);
			}
			for (int i = 0; i < bin_res; ++i) {
				bin_size.push_back(bin_len_lb + 1);
			}
			cout << "bin_size" << bin_size;
		}

		int set_num = my::n_choose_k(codeword_bin_num, info_partition_num);
		//cout << "set_num = " << set_num << endl;


		// this is simpler, the simpler the more efficient

		// generete info pos bias
		Matrix<int> info_pos_bias(1, codeword_bin_num + 1, 'v');
		info_pos_bias.push_back(0);
		for (int i = 0; i < codeword_bin_num; ++i) {
			info_pos_bias.push_back(bin_size(i) + info_pos_bias.back());
		}
		cout << "info_pos_bias" << info_pos_bias;

		Matrix<int> placement_parteen = Matrix_common::generating_all_n_choose_k_pattern(codeword_bin_num, info_partition_num);

		permute_set.resize(1, set_num, false);
		permute_set.reset(Matrix<int>(1, n, 'v'));

		continous_info_bin.resize(1, set_num, false);
		continous_info_bin.reset(Matrix<info_bin_continous_record>(1, info_partition_num, 'v'));

		Matrix<Matrix<int>> info_bin_res_collect(1, set_num);
		info_bin_res_collect.reset(Matrix<int>(1, k, 'v'));

		for (int i = 0; i < set_num; ++i) {
			for (int j = 0; j < info_partition_num; ++j) {
				// place the continus part in this loop
				int bin_ind = placement_parteen(i, j);
				int info_bias_now = info_pos_bias(bin_ind);

				Matrix<int> random_deliver(1, info_part_size(j), 'N');

				if (bin_size(bin_ind) > info_part_size(j)) {
					Matrix<int> tmp_info_bin_res_collect;
					Matrix_common::n_randomly_divide_into_k_and_res(bin_size(bin_ind), info_part_size(j), \
						random_deliver, tmp_info_bin_res_collect);

					tmp_info_bin_res_collect.sort('<');
					for (int p = 0, pmax = tmp_info_bin_res_collect.size(); p < pmax; ++p) {
						info_bin_res_collect(i).push_back(tmp_info_bin_res_collect(p) + info_bias_now);
					}

					random_deliver.sort('<');
				}

				for (int p = 0, pmax = random_deliver.size(); p < pmax; ++p) {
					permute_set(i).push_back(random_deliver(p) + info_bias_now);
				}

				continous_info_bin(i).push_back(info_bin_continous_record(info_bias_now, info_bias_now + bin_size(bin_ind) - 1));
			}
		}

		cout << "continous_info_bin" << continous_info_bin;

		// merge the continous info bin

		Matrix<info_bin_continous_record> tmp_record(1, info_partition_num, 'v');
		for (int i = 0; i < set_num; ++i) {
			tmp_record.resize(1, 0, false);
			tmp_record.push_back(continous_info_bin(i)(0));
			for (int j = 1; j < info_partition_num; ++j) {
				// check if we can merge the next bin
				if (tmp_record.back().end_ind + 1 == continous_info_bin(i)(j).start_ind) {
					// merge the bin
					tmp_record.back().end_ind = continous_info_bin(i)(j).end_ind;
				}
				else {
					tmp_record.push_back(continous_info_bin(i)(j));
				}
			}
			continous_info_bin(i) = tmp_record;			// store the result to continous_info_bin
		}

		cout << "continous_info_bin (after merge)" << continous_info_bin;

		cout << "permute_set" << permute_set;

		cout << "info_bin_res_collect" << info_bin_res_collect;

		// add redundancy positions

		// first add the position in info_bin_res_collect

		for (int i = 0; i < set_num; ++i) {
			for (int j = 0, jmax = info_bin_res_collect(i).size(); j < jmax; ++j) {
				permute_set(i).push_back(info_bin_res_collect(i)(j));
			}
		}

		// then add the positions left
		for (int i = 0; i < set_num; ++i) {
			int permute_set_i_ind = 0;
			int info_bin_res_collect_i_ind = 0;
			int info_bin_res_collect_i_size = info_bin_res_collect(i).size();
			for (int j = 0; j < n; ++j) {
				if (permute_set_i_ind < k && j == permute_set(i)(permute_set_i_ind)) {
					++permute_set_i_ind;
				}
				else if (info_bin_res_collect_i_ind < info_bin_res_collect_i_size \
					&& j == info_bin_res_collect(i)(info_bin_res_collect_i_ind)) {

					++info_bin_res_collect_i_ind;		// also skipping info_bin_res_collect since it inserted

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
		match_pos_num.resize(1, set_num, false);
		col_ind_to_reliability.resize(1, n, false);
		permute_target_col_ind_to_reliability.resize(1, n, false);
		permute_target_reliability_to_col_ind.resize(1, n, false);
		redundancy_set_reliability.resize(1, red, false);
		info_set_reliability_heap.resize(1, k, false);
	}
};
