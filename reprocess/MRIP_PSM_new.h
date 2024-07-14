/*****************************************************************//**
 * \file   MRIP_PSM_new.h
 * \brief  upgrade of MRIP_PSM in cycle info set matching
 * 
 * \author 26259
 * \date   August 2023
 *********************************************************************/

#pragma once
#include"./reprocess_common.h"
using namespace std;

// this is for cycle codes, where generator matrix or parity matrix can shift, 
// e.g., BCH(127,64,21), process on its dual code since rate > 1/2
class PreStored_Matrix_red_cycle {

public:
	int n;
	int k;
	int red;
	Matrix<Matrix<int>> permute_set;
	
	// for select_G_set_continous_count
	Matrix<int> k_distribution_in_n;
	Matrix<int> k_distribution_in_n_permuted;
	Matrix<int> k_distribution_in_n_accumulated;
	Matrix<int> common_pos_num;
	Matrix<int> common_pos_max_shift;
	Matrix<GF2> G_set;										// to be changed into parity-check part, the storage is different
	Matrix<int> natural;

	Basis_update BU;

	// fetch it through G_taregt, by permute the columns of G_target with columns of redundency set, then do Gaussian elemination
	// return: G_target, permute_target and the max common on G_target elements on start
	int get_MRIP_sys_G(const Matrix<int>& reliability_to_col_ind, Matrix<GF2>& G_target, Matrix<int>& permute_target) {
		// in this function we start from redundant set and find the least reliability number of the column index
		// this is of less hardware complexity, but computation cost is larger than 'get_MRIP_sys_G'

		//select_G_set(reliability_to_col_ind);
		target_ind_and_shift opt_info = select_G_set_continous_count(reliability_to_col_ind);
		//int target_ind = opt_info.target_ind;					// this will no greater than m, hence overflow shouldn't be worried
		//int shift = opt_info.shift;
		//shift <<= target_ind;									// this is the solution

		//cout << "target_ind = " << target_ind << endl;
		//cout << "shift = " << shift << endl;

		//cout << "permute_set(target_ind)" << permute_set(target_ind);

		//cout << "common_pos_num(target_ind) = " << common_pos_num(target_ind) << endl;
		G_target = G_set;
		//natural.col_shift_right_cir(n - shift, permute_target);	// initialize the answer
		//cout << "permute_target (1)" << permute_target;

		//permute_target.permute(permute_set(target_ind));		// this is right
		//cout << "permute_target (2)" << permute_target;

		//cout << "G_target" << G_target;
		//cout << "permute_target" << permute_target;

		// a more straight forward way

		natural.set_natural();
		natural.permute(permute_set(opt_info.target_ind));
		natural.col_shift_right_cir(n - opt_info.shift, permute_target);

		return BU.solve(reliability_to_col_ind, G_target, permute_target);									// this is the problem
	}

	// select G that systematic position contains most position in sorted_pos
	// sorted_pos should indicates positions with ordered reliability from most to least
	target_ind_and_shift select_G_set_continous_count(const Matrix<int>& reliability_to_col_ind) {

		//cout << "reliability_to_col_ind" << reliability_to_col_ind;

		// count the common info number for each info set in  permute_set, not absolutely accurate although
		common_pos_num.reset(0);
		common_pos_max_shift.reset(0);

		k_distribution_in_n.reset(0);
		for (int i = 0; i < k; ++i) {
			k_distribution_in_n(reliability_to_col_ind(i)) = 1;
		}
		//cout << "k_distribution_in_n" << k_distribution_in_n;

		for (int g = 0, gmax = permute_set.size(); g < gmax; ++g) {	// the outer for loop is parallizable

			k_distribution_in_n_permuted = k_distribution_in_n;
			k_distribution_in_n_permuted.permute(permute_set(g));	// different automorphism corresponding to a permute_set

			// set the distribution of MRP
			k_distribution_in_n_accumulated.resize(1, n + 1);
			k_distribution_in_n_accumulated.reset(0);		// the first position is sure to be 0
			for (int i = 0; i < n; ++i) {
				k_distribution_in_n_accumulated(i + 1) = k_distribution_in_n_permuted(i);
			}
			//cout << "k_distribution_in_n_accumulated" << k_distribution_in_n_accumulated;


			// set the post fix of k_distribution_in_n with length k
			for (int i = 1; i <= k; ++i) {
				k_distribution_in_n_accumulated.push_back(k_distribution_in_n_accumulated(i));
			}
			//cout << "k_distribution_in_n_accumulated (add_cyc_post_fix)" << k_distribution_in_n_accumulated;

			// accumulate the MRP
			for (int i = 1, imax = n + k; i <= imax; ++i) {
				k_distribution_in_n_accumulated(i) += k_distribution_in_n_accumulated(i - 1);
					// maybe need to think of parallel implementation
			}
			//cout << "k_distribution_in_n_accumulated (final)" << k_distribution_in_n_accumulated;

			for (int p = 0; p < n; ++p) {
				// this will save a lot complexity for selecting G 
				int common_pos_num_tmp = k_distribution_in_n_accumulated(p + k) - k_distribution_in_n_accumulated(p);
					// anding end_ind by one for accumulated the MRP

				if (common_pos_num(g) >= common_pos_num_tmp);
				else{
					// update if a better info set matching is found
					common_pos_num(g) = common_pos_num_tmp;
					common_pos_max_shift(g) = p;
				}
			}
		}

		// first sort the reliability_to_col_ind

		//cout << "common_pos_num" << common_pos_num;
		//cout << "common_pos_max_shift" << common_pos_max_shift;

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
	PreStored_Matrix_red_cycle(const Matrix<GF2>& G) : n(G.col()), k(G.row()), red(n - k), BU(n, k) {

		//cout << "k = " << k << endl;
		//cout << "n - k = " << n - k << endl;

		// this is simpler, the simpler the more efficient

		// generate set_num systematic matrices
		G_set = G;

		Matrix<int> permute_record = G_set.GE_left_identity_4_GF2_echelon();

		cout << "G_set" << G_set << endl;

		natural = Matrix<int>(1, n, 'N');

		// determine how many multiply number we can access
		int set_num = 0;

		// only for GF2 code
		while ((1 << (set_num + 1)) % n != 1) {
			set_num++;
		}
		set_num++;						// start from 2^0, which is in count
		cout << "set_num = " << set_num << endl;		// set_num = m for BCH code of length 2^m-1

		// permute_set is exactly the set that can be shifted and generate different info set, that is a magic
		permute_set.resize(1, set_num, false);
		for (int i = 0; i < set_num; ++i) {
			permute_set(i).resize(1, n);
			for (int j = 0; j < n; ++j) {
				permute_set(i)(j) = (j << i) % n;
			}
		}

		cout << "permute_set" << permute_set;

		k_distribution_in_n.resize(1, n, false);
		k_distribution_in_n_permuted.resize(1, n, false);
		k_distribution_in_n_accumulated.resize(1, n + k + 1, false);		// consider storing cyclic postfix
		common_pos_num.resize(1, set_num, false);
		common_pos_max_shift.resize(1, set_num, false);
	}
};