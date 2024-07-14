/*****************************************************************//**
 * \file   Viterbi_naive.h
 * \brief  this file contains naive implementation of Viterbi alogrithm, they have limit on both
 * complexity and type of codes, you have better never use them.
 * 
 * \author 26259
 * \date   April 2023
 *********************************************************************/

#pragma once

#include "Viterbi_common.h"

class Viterbi_cyclic_code_k {		// for fixed cyclic codes, where k < n - k, no list decoding
private:
	int n;
	int k;
	int r;
	Matrix<GF2> PM;

	int max_list_size;
	Matrix<GF2> list_v_hat;

	int pow_2_k;
	int half_pow_2_k;
	Matrix<int> past_state;
	// the state before the current state that connect to this state and has the minimun path metric
	// past_state(i,j) is best state connecting to i at time instant j, current state is j+1

	Matrix<my_double> path_metric;			// path_metric(i,j) is path metric of state i at time instance j+1

	Matrix<int> PM_col_reps;
	// the state correspond to j^{th} received bits, indicated by j^{th} column to a parity matrix PM, now j=0
	// state is indicated by int numner form by add and shifting of binary column of PM, e.g. [1,1,0]^T -> 6

	Matrix<my_double> dist;// initialize correlation distance, considering both case for soft metric and hard metric


	void initialize_parameter() {
		// testing
		/*PM.switch_row(0, 2);
		PM.switch_row(1, 3);
		cout << "PM" << PM;*/

		pow_2_k = 1 << k;
		half_pow_2_k = pow_2_k >> 1;
		dist = Matrix<my_double>(2, n);

		if (r <= k) {
			// under this case we cannot use viterbi algorithm, we use brutule force
			return;
		}
		PM_col_reps = Matrix<int>(1, n);
		past_state = Matrix<int>(pow_2_k, n);
		path_metric = Matrix<my_double>(pow_2_k, n);

		// compute the state correspond to i^{th} received bits, named it as PM column representation
		PM_col_reps(0) = 1 << (k - 1); 	// special for cyclic code since the PM's first column is [1,0,0,...,0]^T
		for (int i = 1; i < k; ++i) {
			// for the first k-1 time instance, not all state is valid, computation can be saved
			PM_col_reps(i) = 0;
			for (int p = 0; p <= i; ++p) {		// only 2^i state valid
				PM_col_reps(i) <<= 1;
				PM_col_reps(i) += (int)PM(p, i);
			}
			PM_col_reps(i) <<= (k - i - 1);		// shift to the correct state
		}

		// compute the state correspond to i^{th} received bits, named it as PM column representation
		int tmp_PM_col_rep = 0;		// this is a static state column for i in [k,n-k]
		for (int p = 1; p <= k; ++p) {
			tmp_PM_col_rep <<= 1;
			tmp_PM_col_rep += (int)PM(p, k);
		}
		for (int i = k; i <= n - k; ++i) {		// designed for cyclic code
			PM_col_reps(i) = tmp_PM_col_rep;
		}
		for (int i = n - k; i < n; ++i) {
			PM_col_reps(i) = 0;
			for (int p = r + i - n; p < r; ++p) {
				PM_col_reps(i) <<= 1;
				PM_col_reps(i) += (int)PM(p, i);
			}
		}
	}

public:
	Viterbi_cyclic_code_k(const Matrix<GF2>& _PM) {
		// get the code parameters, for a linear block code of (n,k) and r = n - k, 
		// with parity check matrix PM, generator check matrix GM

		n = _PM.col();
		r = _PM.row();
		k = n - r;
		PM = _PM;

		initialize_parameter();
	}

	/**
	 * .apply decode method to 'r_or_hdr', under BPSK with 0 -> 1 and 1 -> -1
	 * optimized for case k<n-k, the decoder needs 2^{k} states, for list_size > 1, not developed
	 *
	 * \param r_or_hdr: received codeword, can be soft or hard
	 * \return decoded result, u
	 */
	template<class T>
	Matrix<GF2> decode_v(const Matrix<T>& r_or_hdr) {
#ifdef RUN_MSG
#ifdef count_operation_number
#ifdef use_my_double

		int double_ope_num_before = my_double_auxiliary_storage::operation_number;
		int double_ope_num_after;
#endif // use_my_double

		int GF2_ope_num_before = GF2_auxiliary_storage::operation_number;
		int GF2_ope_num_after;
#endif // count_operation_number
#endif // RUN_MSG

		// design for rate < 1/2, this is not stupid, use brutle force is slightly more complex than this.

		// but, brutle force is easy to extend to list decoding. you may not write list decoding of this type of viterbi. For now.

		// to merge the method of hard and soft decoding, if T is ordered structure, then use soft decoidng, else use hard deocding
		bool is_soft = !(T() < T());
		//cout << "is_soft = " << is_soft << endl;


		//cout << "pow_2_k = " << pow_2_k << endl;

		past_state.reset(-invalid_large_metric);	// -invalid_large_metric means invalid past state, set them all invalid initially
		path_metric.reset(invalid_large_metric);	// invalid_large_metric means invalid correlation distance, large number is perfered

		// compute distance for each symbol to 0 and 1
		for (int i = 0; i < n; ++i) {
			// not key complexity
			T ri = r_or_hdr(i);

			my_double d0 = (double)(is_soft ? (ri - 1) : (ri != 0));
			my_double d1 = (double)(is_soft ? (ri + 1) : (ri == 0));

			dist(0, i) = d0 * d0;
			dist(1, i) = d1 * d1;
		}

		// for i==0, the first received symbol
		path_metric(0, 0) = dist(0, 0);
		path_metric(PM_col_reps(0), 0) = dist(1, 0);

		// initialize past state
		past_state(0, 0) = 0;					// initial and final state is 0 for a linear block code, set the past state be 0
		past_state(PM_col_reps(0), 0) = 0;		// second path derived from initial 0 state

		//ope_num_after = my_double_auxiliary_storage::operation_number;
		//cout << "ope_num = " << ope_num_after - ope_num_before << endl;

		// starting trellises, special for bch code (cyclic code), i is the time instance in [1,r-1]
		for (int i = 1; i < k; ++i) {
			// compute the state correspond to i^{th} received bits, named it as PM column representation

			int jump_state_num = 1 << (k - i);
			// j is the state index in [0,2^r-1], but can be jumped over invalid states
			for (int j = 0; j < pow_2_k; j += jump_state_num) {
				// consider one state in the for lop, i.e. state j, actually no contention happened

				path_metric(j, i) = path_metric(j, i - 1) + dist(0, i);
				past_state(j, i) = j;

				path_metric(j ^ PM_col_reps(i), i) = path_metric(j, i - 1) + dist(1, i);
				past_state(j ^ PM_col_reps(i), i) = j;
			}
		}

		//ope_num_after = my_double_auxiliary_storage::operation_number;
		//cout << "ope_num = " << ope_num_after - ope_num_before << endl;


		// trellises with full state, i is the time instance in [k,n-k], like shift register to reduced state
		for (int i = k; i < n - k; ++i) {

			// j is the state index in [0,2^r-1]
			for (int j = 0; j < half_pow_2_k; ++j) {
				// actually no contention and no divergence
				int j_prime = j << 1;

				path_metric(j_prime, i) = path_metric(j, i - 1) + dist(0, i);
				past_state(j_prime, i) = j;
			}
			for (int j = half_pow_2_k; j < pow_2_k; ++j) {
				// actually no contention and no divergence

				// throw away the heighest bit, and shift right 1 bit
				int j_prime = (j ^ half_pow_2_k) << 1;

				path_metric(j_prime ^ PM_col_reps(i), i) = path_metric(j, i - 1) + dist(1, i);
				past_state(j_prime ^ PM_col_reps(i), i) = j;
			}
		}

		//ope_num_after = my_double_auxiliary_storage::operation_number;
		//cout << "ope_num = " << ope_num_after - ope_num_before << endl;

		// ending trellises, special for bch code (cyclic code), i is the time instance in [n-k,n-1]
		for (int i = n - k; i < n; ++i) {

			int range = pow_2_k >> (i - n + k + 1);
			// j is the state index in [0,range]
			for (int j = 0; j < range; ++j) {
				// consider one state in the for lop, i.e. state j

				// update correlation distance and past state
				my_double can_1 = path_metric(j, i - 1) + dist(0, i);		// accumulate correlation distance for each path to state j
				my_double can_2 = path_metric(j ^ PM_col_reps(i), i - 1) + dist(1, i);

				// choose the best path for each merging state
				bool is_can_1_better = can_1 < can_2;
				path_metric(j, i) = is_can_1_better ? can_1 : can_2;
				past_state(j, i) = is_can_1_better ? j : (j ^ PM_col_reps(i));
			}
		}

		//cout << "path_metric" << path_metric;
		//cout << "past_state" << past_state;

		// get the decoded codeword, final state is 0
		Matrix<GF2> v_hat_viterbi(1, n);		// the result codeword
		int track_state = 0;					// track the current state, for final state, it must be 0
		for (int i = n - 1; i >= n - k; --i) {		// path back tracking from time instance n-1 to 0
			v_hat_viterbi(i) = past_state(track_state, i) != track_state;		// if state changed, the decoded bit is 1 else it is 0
			track_state = past_state(track_state, i);							// update the current state a time instant earlier
		}
		for (int i = n - k - 1; i >= k; --i) {
			v_hat_viterbi(i) = (past_state(track_state, i) & half_pow_2_k) == half_pow_2_k;
			track_state = past_state(track_state, i);
		}
		for (int i = k - 1; i >= 0; --i) {		// path back tracking from time instance n-1 to 0
			v_hat_viterbi(i) = past_state(track_state, i) != track_state;		// if state changed, the decoded bit is 1 else it is 0
			track_state = past_state(track_state, i);							// update the current state a time instant earlier
		}



#ifdef RUN_MSG
		cout << "path_metric(0, n - 1) = " << path_metric(0, n - 1) << endl;
#ifdef count_operation_number
#ifdef use_my_double

		double_ope_num_after = my_double_auxiliary_storage::operation_number;
		cout << "(Viterbi_cyclic_code_k) double_ope_num = " << double_ope_num_after - double_ope_num_before << endl;
#endif // use_my_double

		GF2_ope_num_after = GF2_auxiliary_storage::operation_number;
		cout << "(Viterbi_cyclic_code_k) GF2_ope_num = " << GF2_ope_num_after - GF2_ope_num_before << endl;
#endif // count_operation_number
#endif // RUN_MSG

		return v_hat_viterbi;
	}

	
};

class Viterbi_cyclic_code_r {
	// this is not the very original version, for that, you have to find the code of old version,
	// see the file in  ../_old_vertion/04-04/reprocess/Viterbi.h, class of Viterbi, function decode_v_cyclic_code_r
	// optimized specially for n-k<k, i.e. code rate >= 1/2
private:
	int n;
	int r;
	int pow_2_r;
	Matrix<GF2> PM;

	Matrix<GF2> list_v_hat;
	Matrix<int> PM_col_reps;
	// the state correspond to j^{th} received bits, indicated by j^{th} column to a parity matrix PM, now j=0
	// state is indicated by int numner form by add and shifting of binary column of PM, e.g. [1,1,0]^T -> 6

	Matrix<my_double> dist;// initialize correlation distance, considering both case for soft metric and hard metric

	/* list decoding */
	Matrix<int> last_opt_state;				// store state at time instant [1,n-1], the initial state 0 is not stored
	Matrix<my_double> relative_metric;		// for each list, need a bias at merge point for metric unification
	priority_queue<metric_point> can_opt;	// candidate merge point for next optimum path
	int list_proceed;
	int diverge_time;

	Matrix<past_and_metric> state;			// consider using Matrix of Matrix for space saving, it will be smarter
	Matrix<int> ending_1_ind;				// store the ending 1 index, -1 means not cared
	Matrix<int> starting_1_ind;				// store the starting 1 index, -1 means not cared

	int min_contention_col_ind;
	Matrix<GF2> unused_PM;
	Matrix<int> valid_list_ind;

public:
	Matrix<GF2> test_zero;

	Viterbi_cyclic_code_r(const Matrix<GF2>& _PM) {
		// get the code parameters, for a linear block code of (n,k) and r = n - k, 
		// with parity check matrix PM

		n = _PM.col();
		r = _PM.row();
		PM = _PM;

		/* apply space */

		PM_col_reps = Matrix<int>(1, n);
		dist = Matrix<my_double>(2, n);

		/* list decoding */
		last_opt_state = Matrix<int>(1, n);

		pow_2_r = 1 << r;
		state = Matrix<past_and_metric>(pow_2_r, n);

		for (int i = 0; i < n; ++i) {
			PM_col_reps(i) = 0;
			for (int p = 0; p < r; ++p) {
				PM_col_reps(i) <<= 1;
				PM_col_reps(i) += (int)PM(p, i);
			}
		}
		min_contention_col_ind = r;		// if all divergence, min_contention_col_ind is r

		ending_1_ind = Matrix<int>(1, n);
		ending_1_ind.reset(-1);
		int divergence_ending_1_ind = 0;
		ending_1_ind(0) = divergence_ending_1_ind;
		for (int i = 1; i < n; ++i) {
			if (PM(divergence_ending_1_ind + 1, i) == 1) {
				// divergence
				divergence_ending_1_ind++;
				ending_1_ind(i) = divergence_ending_1_ind;
			}
			else {
				min_contention_col_ind = my::min(min_contention_col_ind, i);
				ending_1_ind(i) = ending_1_ind(i - 1);
			}
			if (divergence_ending_1_ind == r - 1) {
				// it must come here after all divergence
				break;
			}
		}

		starting_1_ind = Matrix<int>(1, n);
		starting_1_ind.reset(-1);
		int shrink_ending_1_ind = r - 1;
		starting_1_ind(n - 1) = shrink_ending_1_ind;
		for (int i = n - 2; i >= 0; --i) {
			if (PM(shrink_ending_1_ind - 1, i) == 1) {
				// shrinking
				shrink_ending_1_ind--;
				starting_1_ind(i) = shrink_ending_1_ind;
			}
			else {
				starting_1_ind(i) = starting_1_ind(i + 1);
			}
			if (shrink_ending_1_ind == 0) {
				// it must come here after all divergence
				break;
			}
		}

		// verify
#ifdef RUN_MSG

		cout << "PM" << PM;
		cout << "ending_1_ind" << ending_1_ind;
		cout << "starting_1_ind" << starting_1_ind;
#endif
	}

	/**
	 * .change PM not changing the size of code, i.e., n, and k not changed.
	 * not initializing state, for decode_v_once especially.
	 */
	void change_unused_PM(const Matrix<GF2>& _unused_PM) {
		unused_PM = _unused_PM;
	}

	/**
	 * . find the next sub-optimum codeword, with parameter set in private area
	 */
	void next_subopt_v() {
		// we have problem here

		//cout << "r=" << r << endl;
		//cout << "min_contention_col_ind=" << min_contention_col_ind << endl;

		// consider the space saving later, i.e., we donot need to store the first r state, since no contention
		for (int i = diverge_time; i >= min_contention_col_ind; --i) {
			int los = last_opt_state(i);
			bool is_best_path_0 = list_v_hat(list_proceed, i) == 0;

			// take the worse path, but we should judge whether the worse path exists
			int former_state = is_best_path_0 ? (los ^ PM_col_reps(i)) : los;
			my_double dist_tmp = is_best_path_0 ? dist(1, i) : dist(0, i);

			my_double mt = state(former_state, i - 1).metric;
			if (mt != invalid_large_metric) {

				can_opt.push(metric_point(former_state, i - 1, mt + dist_tmp \
					- state(los, i).metric + relative_metric(list_proceed), los, list_proceed));

			}

			//finished_state(last_opt_state(i), i) = list_proceed;
		}
		//cout << "finished_state" << finished_state;
		//cout << "last_opt_state" << last_opt_state;
		//cout << can_opt.top() << endl;

		/*while (!second_opt.empty()) {
			cout << second_opt.top();
			second_opt.pop();
		}
		cout << endl;*/

		// now you get the second opt metric path, backtracking it
		list_proceed++;
		//cout << "list_proceed=" << list_proceed << endl;
		diverge_time = can_opt.top().diverge_time;
		diverge_time++;		// adding 1 temporary, this become merge bit,for less computation only

		int diverge_list_ind = can_opt.top().list_belonging;
		//int diverge_list_ind = finished_state(can_opt.top().merge_state, diverge_time);
		// copy the path before diverge
		for (int i = n - 1; i > diverge_time; --i) {
			list_v_hat(list_proceed, i) = list_v_hat(diverge_list_ind, i);
		}
		// flip the merge bit
		list_v_hat(list_proceed, diverge_time) = list_v_hat(diverge_list_ind, diverge_time).flip();
		relative_metric(list_proceed) = can_opt.top().relative_metric;		// right

		diverge_time--;		// minus 1, change the variable back
		int track_state = can_opt.top().diverge_state;
		for (int i = diverge_time; i >= 0; --i) {
			last_opt_state(i) = track_state;
			track_state = state(track_state, i).past;							// update the current state a time instant earlier
			list_v_hat(list_proceed, i) = track_state != last_opt_state(i);		// if state changed, the decoded bit is 1, else 0

			/*last_opt_state(i) = track_state;
			list_v_hat(list_proceed, i) = past_state(track_state, i) != track_state;
			track_state = past_state(track_state, i);	*/
		}
		//cout << "list_v_hat" << list_v_hat;

		can_opt.pop();
		//cout << "second_opt.empty() = " << can_opt.empty() << endl;
	}

	/**
	 * .apply decode method to 'r_or_hdr', under BPSK with 0 -> 1 and 1 -> -1
	 * optimized for saving the storage space, if use the parity matrix repeatedly, use this function
	 *
	 * \param r_or_hdr: received codeword, can be soft or hard
	 * \param list_size: number of decoded v listed by row, which are maximum likelihood
	 * \return decoded result, v
	 */
	template<class T>
	Matrix<GF2> decode_v(const Matrix<T>& r_or_hdr, int list_num = 1) {
#ifdef RUN_MSG
#ifdef count_operation_number
#ifdef use_my_double

		int double_ope_num_before = my_double_auxiliary_storage::operation_number;
		int double_ope_num_after;
#endif // use_my_double

		int GF2_ope_num_before = GF2_auxiliary_storage::operation_number;
		int GF2_ope_num_after;
#endif // count_operation_number
#endif // RUN_MSG

		// need some extension, for divergence and shrinking procees with contention

		// to merge the method of hard and soft decoding, if T is ordered structure, then use soft decoidng, else use hard deocding
		bool is_soft = !(T() < T());
		//cout << "is_soft = " << is_soft << endl;

		//cout << "n = " << n << endl;

		// compute distance for each symbol to 0 and 1
		for (int i = 0; i < n; ++i) {
			// not key complexity
			T ri = r_or_hdr(i);

			my_double d0 = (double)(is_soft ? (ri - 1) : (ri != 0));
			my_double d1 = (double)(is_soft ? (ri + 1) : (ri == 0));

			dist(0, i) = d0 * d0;
			dist(1, i) = d1 * d1;
		}
		//cout << "dist" << dist;

		//cout << "PM_col_reps" << PM_col_reps;

		// for time i, state[i] mapping from state number to state metric and past state pair
		state(0, 0) = past_and_metric(0, dist(0, 0));
		state(PM_col_reps(0), 0) = past_and_metric(0, dist(1, 0));

		for (int i = 1; i < n; ++i) {
			if (ending_1_ind(i) != -1 && starting_1_ind(i) == -1) {		// at the region of divergence
				if (ending_1_ind(i) == ending_1_ind(i - 1) + 1) {
					// for the divergence instance, not all state is valid, computation can be saved
					int jump_state_num = 1 << (r - 1 - ending_1_ind(i - 1));
					for (int j = 0; j < pow_2_r; j += jump_state_num) {
						// j is the state index in [0,2^r-1], but can be jumped over invalid states
						// j stands for the former states, scan over the former state in this for loop
						// if the first r column of PM are independent, actually no contention happened

						// update correlation distance and past state
						// choose the best path for each merging state

						my_double metric_store = state(j, i - 1).metric;
						state(j, i) = past_and_metric(j, metric_store + dist(0, i));
						state(j ^ PM_col_reps(i), i) = past_and_metric(j, metric_store + dist(1, i));
					}
				}
				else {
					// we must have ending_1_ind(i) == ending_1_ind(i - 1)
					// contention only, no divergence, j is the present state index
					int jump_state_num = 1 << (r - 1 - ending_1_ind(i - 1));
					for (int j = 0; j < pow_2_r; j += jump_state_num) {
						// update correlation distance and past state
						my_double can_1 = state(j, i - 1).metric + dist(0, i);	// accumulate correlation distance for each path to state j
						my_double can_2 = state(j ^ PM_col_reps(i), i - 1).metric + dist(1, i);	// key complexity

						// choose the best path for each merging state
						state(j, i) = can_1 < can_2 ? past_and_metric(j, can_1) : past_and_metric(j ^ PM_col_reps(i), can_2);
					}
				}
			}
			else if (ending_1_ind(i) == -1 && starting_1_ind(i) == -1) {
				// trellises with full state, i is the time instance in [r,n-r], to be optimized for k<r
				// compute the state correspond to i^{th} received bits, named it as PM column representation			

				// j is the state index in [0,2^r-1]
				for (int j = 0; j < pow_2_r; ++j) {
					// update correlation distance and past state
					my_double can_1 = state(j, i - 1).metric + dist(0, i);		// accumulate correlation distance for each path to state j
					my_double can_2 = state(j ^ PM_col_reps(i), i - 1).metric + dist(1, i);	// key complexity

					// choose the best path for each merging state
					state(j, i) = can_1 < can_2 ? past_and_metric(j, can_1) : past_and_metric(j ^ PM_col_reps(i), can_2);
				}
			}
			else if (ending_1_ind(i) != -1 && starting_1_ind(i) != -1) {
				// this is the most complex case, consider contention and divergence simultaneously
				cout << "warning, entering unexpected states" << endl;

			}
			else if (ending_1_ind(i) == -1 && starting_1_ind(i) != -1) {
				if (starting_1_ind(i) == starting_1_ind(i - 1) + 1) {
					// case of shrink
					// j is the state index in [0,range]
					int range = pow_2_r >> (starting_1_ind(i) + 1);
					for (int j = 0; j < range; ++j) {
						// consider one state in the for lop, i.e. state j

						// update correlation distance and past state
						my_double can_1 = state(j, i - 1).metric + dist(0, i);	// accumulate correlation distance for each path to state j
						my_double can_2 = state(j ^ PM_col_reps(i), i - 1).metric + dist(1, i);

						// choose the best path for each merging state
						state(j, i) = can_1 < can_2 ? past_and_metric(j, can_1) : past_and_metric(j ^ PM_col_reps(i), can_2);
					}
				}
				else {
					// we must have starting_1_ind(i) == starting_1_ind(i - 1)
					// contention only, no shrink, j is the present state index
					int range = pow_2_r >> starting_1_ind(i);
					for (int j = 0; j < range; ++j) {
						// consider one state in the for lop, i.e. state j

						// update correlation distance and past state
						my_double can_1 = state(j, i - 1).metric + dist(0, i);	// accumulate correlation distance for each path to state j
						my_double can_2 = state(j ^ PM_col_reps(i), i - 1).metric + dist(1, i);

						// choose the best path for each merging state
						state(j, i) = can_1 < can_2 ? past_and_metric(j, can_1) : past_and_metric(j ^ PM_col_reps(i), can_2);
					}
				}

			}
		}


		// this is the original altorithm for cyclic codes, simple but has lots of limitation
		//// starting trellises, special for bch code (cyclic code), i is the time instance in [1,r-1]
		//for (int i = 1; i < r; ++i) {
		//	// compute the state correspond to i^{th} received bits, named it as PM column representation
		//	// for the first r-1 time instance, not all state is valid, computation can be saved
		//	int jump_state_num = 1 << (r - i);
		//	// j is the state index in [0,2^r-1], but can be jumped over invalid states
		//	for (int j = 0; j < pow_2_r; j += jump_state_num) {
		//		// consider one state in the for lop, i.e. state j, actually no contention happened
		//		// if the first r column of PM are independent
		//		// update correlation distance and past state
		//		//my_double can_1 = path_metric(j, i - 1) + dist_to_0;		// accumulate correlation distance for each path to state j
		//		//my_double can_2 = path_metric(j ^ PM_col_rep, i - 1) + dist_to_1;
		//		// choose the best path for each merging state
		//		//bool is_can_1_better = can_1 < can_2;
		//		//path_metric(j, i) = is_can_1_better ? can_1 : can_2;
		//		//past_state(j, i) = is_can_1_better ? j : (j ^ PM_col_rep);
		//		my_double metric_store = state(j, i - 1).metric;
		//		state(j, i) = past_and_metric(j, metric_store + dist(0, i));
		//		state(j ^ PM_col_reps(i), i) = past_and_metric(j, metric_store + dist(1, i));
		//	}
		//}
		////ope_num_after = my_double_auxiliary_storage::operation_number;
		////cout << "ope_num = " << ope_num_after - ope_num_before << endl;
		//// trellises with full state, i is the time instance in [r,n-r], to be optimized for k<r
		//for (int i = r; i < n - r; ++i) {
		//	// compute the state correspond to i^{th} received bits, named it as PM column representation			
		//	// j is the state index in [0,2^r-1]
		//	for (int j = 0; j < pow_2_r; ++j) {
		//		// update correlation distance and past state
		//		my_double can_1 = state(j, i - 1).metric + dist(0, i);		// accumulate correlation distance for each path to state j
		//		my_double can_2 = state(j ^ PM_col_reps(i), i - 1).metric + dist(1, i);	// key complexity
		//		// choose the best path for each merging state
		//		state(j, i) = can_1 < can_2 ? past_and_metric(j, can_1) : past_and_metric(j ^ PM_col_reps(i), can_2);
		//	}
		//}
		////ope_num_after = my_double_auxiliary_storage::operation_number;
		////cout << "ope_num = " << ope_num_after - ope_num_before << endl;
		//// ending trellises, special for bch code (cyclic code), i is the time instance in [n-r,n-1]
		//for (int i = n - r; i < n; ++i) {
		//	int range = pow_2_r >> (i - n + r + 1);
		//	// j is the state index in [0,range]
		//	for (int j = 0; j < range; ++j) {
		//		// consider one state in the for lop, i.e. state j
		//		// update correlation distance and past state
		//		my_double can_1 = state(j, i - 1).metric + dist(0, i);		// accumulate correlation distance for each path to state j
		//		my_double can_2 = state(j ^ PM_col_reps(i), i - 1).metric + dist(1, i);
		//		// choose the best path for each merging state
		//		state(j, i) = can_1 < can_2 ? past_and_metric(j, can_1) : past_and_metric(j ^ PM_col_reps(i), can_2);
		//	}
		//}

		list_v_hat.resize(list_num, n, false);
		relative_metric.resize(1, list_num, false);

		// set for listing

		// get the optimum decoded codeword, final state is 0
		int track_state = 0;					// track the current state, for final state, it must be 0
		for (int i = n - 1; i >= 0; --i) {		// path back tracking from time instance n-1 to 0
			last_opt_state(i) = track_state;
			track_state = state(track_state, i).past;
			list_v_hat(0, i) = track_state != last_opt_state(i);

			//list_v_hat(0, i) = past_state(track_state, i) != track_state;	// if state changed, the decoded bit is 1 else it is 0			
			//track_state = past_state(track_state, i);						// update the current state a time instant earlier
		}
		//cout << "list_v_hat" << list_v_hat;

		valid_list_ind.resize(1, 0, false);		// not releasing the space, just seting row and column
		bool is_unused_PM_exist = unused_PM.size() != 0;

		/* since the parity matrix is determined, we if unused_PM exists, it must be necessary */
		if (!is_unused_PM_exist || (is_unused_PM_exist && unused_PM.check_inner_product(list_v_hat.get_row(0)))) {
			valid_list_ind.push_back(0);
		}

		if (list_num > 1) {
			relative_metric(0) = 0;
			can_opt = priority_queue<metric_point>();
			list_proceed = 0;
			diverge_time = n - 1;

			for (int i = 1; i < list_num; ++i) {
				next_subopt_v();

				/* since the parity matrix is determined, we if unused_PM exists, it must be necessary */
				if (!is_unused_PM_exist || (is_unused_PM_exist && unused_PM.check_inner_product(list_v_hat.get_row(i)))) {
					valid_list_ind.push_back(i);
				}
			}
		}


#ifdef RUN_MSG
		//cout << "relative_metric" << relative_metric;
		
		//while (!can_opt.empty()) {
		//	cout << can_opt.top();
		//	can_opt.pop();
		//}
		//cout << endl;

		//cout << "list_v_hat" << list_v_hat;
#ifdef count_operation_number
#ifdef use_my_double

		double_ope_num_after = my_double_auxiliary_storage::operation_number;
		cout << "(Viterbi_cyclic_code_r) double_ope_num = " << double_ope_num_after - double_ope_num_before << endl;
#endif // use_my_double

		GF2_ope_num_after = GF2_auxiliary_storage::operation_number;
		cout << "(Viterbi_cyclic_code_r) GF2_ope_num = " << GF2_ope_num_after - GF2_ope_num_before << endl;
#endif // count_operation_number
#endif // RUN_MSG

		return list_v_hat.get_rows(valid_list_ind);
	}
};
