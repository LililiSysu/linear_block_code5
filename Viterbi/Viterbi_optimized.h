/*****************************************************************//**
 * \file   Viterbi_optimized.h
 * \brief  this is of the most efficient Viterbi algorithm we developed, it need to find a optimized Parity Matrix
 * first, which should satisfy final condition (see the code), and is with minimun computation estimation
 * the developed algorithm can only deal with the final-condition Parity Matrix, but in our observation, any linear
 * block code, with permutation, can generate a final-condition Parity Matrix to use
 *
 * \author 26259
 * \date   April 2023
 *********************************************************************/

#pragma once

#include"Viterbi_common.h"
#include<set>

 /**
  * .only for the Parity check Matrix satisfying the final condition
  */
class Viterbi_optimized {
private:
	int n;
	int r;
	int pow_2_r;
	Matrix<GF2> PM;

	Matrix<GF2> list_v_hat;
	Matrix<int> PM_col_reps;
	// the state correspond to j^{th} received bits, indicated by j^{th} column to a parity matrix PM, now j=0
	// state is indicated by int numner form by add and shifting of binary column of PM, e.g. [1,1,0]^T -> 6

	Matrix<my_double> dist;					// initialize correlation distance, considering both case for soft metric and hard metric

	/* list decoding */
	Matrix<int> last_opt_state;				// store state at time instant [1,n-1], the initial state 0 is not stored
	Matrix<my_double> relative_metric;		// for each list, need a bias at merge point for metric unification

	Heap_max<metric_point_opt> can_opt_Heap_max;
	Heap_min<metric_point_opt> can_opt_Heap_min;
	int reduncency_in_Heap_max;
	int reduncency_in_Heap_min;

	set<metric_point_opt> can_opt_set;				// please do not use advanced data structure, since it's slow

	Matrix<bool> is_state_1_or_3;					// record state, for simplification, but not easy to read
	Matrix<int> state_2_or_4_time_instance;			// this is for list viterbi

	int list_proceed;
	int diverge_time;

	Matrix_flex_col<past_and_metric> state; // consider using Matrix_flex_col for space saving, it will be smarter	
	Matrix<int> time_instance_state;		// 0 -> unset, 1 -> divergence, 2 -> no divergence, no shrinking, contention, 
											// 3 -> divergence and shrinking, no contention, 4 -> shrinking

	int min_contention_col_ind;
	Matrix<GF2> unused_PM;
	Matrix<int> valid_list_ind;

public:

	Viterbi_optimized(const Matrix<GF2>& _PM) {
		// get the code parameters, for a linear block code of (n,k) and r = n - k, 
		// with parity check matrix PM

		n = _PM.col();
		r = _PM.row();
		PM = _PM;

		/* apply space */
		dist = Matrix<my_double>(2, n);

		/* list decoding */
		last_opt_state = Matrix<int>(1, n);
		pow_2_r = 1 << r;
		min_contention_col_ind = r;				// if all divergence, min_contention_col_ind is r

		Matrix<int> ending_1_ind(1, n);		// store the ending 1 index
		ending_1_ind.reset(r - 1);
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

		Matrix<int> starting_1_ind(1, n);		// store the starting 1 index, -1 means not cared
		starting_1_ind.reset(0);
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

#ifdef RUN_MSG
		// verify
		cout << "PM" << PM;
		cout << "ending_1_ind" << ending_1_ind;
		cout << "starting_1_ind" << starting_1_ind;
#endif // RUN_MSG

		Matrix<int> column_size(1, n);
		for (int i = 0; i < n - 1; ++i) {
			int tmp = starting_1_ind(i) == starting_1_ind(i + 1) ? starting_1_ind(i) - 1 : starting_1_ind(i);
			column_size(i) = 1 << (ending_1_ind(i) - tmp);
		}
		column_size(n - 1) = 1;

#ifdef RUN_MSG
		cout << "column_size" << column_size;
#endif // RUN_MSG


		state.resize(column_size);

		PM_col_reps = Matrix<int>(1, n);		// set the representation suitable for state2
		for (int i = 0; i < n - 1; ++i) {
			PM_col_reps(i) = 0;
			int tmp = starting_1_ind(i) == starting_1_ind(i + 1) ? starting_1_ind(i) : starting_1_ind(i) + 1;
			for (int p = tmp; p <= ending_1_ind(i); ++p) {
				PM_col_reps(i) <<= 1;
				PM_col_reps(i) += (int)PM(p, i);
			}
		}
		PM_col_reps(n - 1) = 0;

#ifdef RUN_MSG
		cout << "PM_col_reps" << PM_col_reps;
#endif // RUN_MSG


		// initialize time instance state
		time_instance_state = Matrix<int>(1, n);
		is_state_1_or_3 = Matrix<bool>(1, n);
		state_2_or_4_time_instance = Matrix<int>(1, n, 'v');
		time_instance_state(0) = 1;
		for (int i = 1; i < n; ++i) {

			int former_state_size = state.col(i - 1);
			int present_state_size = state.col(i);

			if (present_state_size > former_state_size) {
				// divergence
				time_instance_state(i) = 1;
				is_state_1_or_3(i) = true;
			}
			else if (present_state_size == former_state_size && ending_1_ind(i) == ending_1_ind(i - 1)) {
				//  no divergence, no shrinking, contention
				time_instance_state(i) = 2;
				is_state_1_or_3(i) = false;
				state_2_or_4_time_instance.push_back(i);
			}
			else if (present_state_size == former_state_size && ending_1_ind(i) != ending_1_ind(i - 1)) {
				// divergence and shrinking, no contention, the special case
				time_instance_state(i) = 3;
				is_state_1_or_3(i) = true;
			}
			else if (present_state_size < former_state_size) {
				// shrinking
				time_instance_state(i) = 4;
				is_state_1_or_3(i) = false;
				state_2_or_4_time_instance.push_back(i);
			}
		}

#ifdef RUN_MSG
		cout << "time_instance_state" << time_instance_state;
		cout << "is_state_1_or_3" << is_state_1_or_3;
		cout << "state_2_or_4_time_instance" << state_2_or_4_time_instance;
#endif // RUN_MSG

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
	 * using data sturcture of set, the red-black tree, do some simpilification, but not easy to read or modify
	 * but the result is wrong
	 */
	void next_subopt_v(int to_be_found_list_num) {
		// we consider limit the candidate optimum metric point number up to the number of to-be-found list vector
		// consider switch to priority queue, using set is slow, but priority queue is slow too
		// difference of list output is induced by disorder of same relative metric, doesn't matter

		for (int k = state_2_or_4_time_instance.binary_search(diverge_time, true); k >= 0; --k) {
			int i = state_2_or_4_time_instance(k);
			int los = last_opt_state(i);
			bool is_best_path_0 = list_v_hat(list_proceed, i) == 0;
			int former_state = los;
			former_state ^= is_best_path_0 ? PM_col_reps(i) + (time_instance_state(i) == 4 ? state.col(i) : 0) : 0;

			my_double can_relative_metric = state(i - 1, former_state).metric + dist(is_best_path_0, i) \
				- state(i, los).metric + relative_metric(list_proceed);

			if (can_opt_Heap_max.size() - reduncency_in_Heap_max >= to_be_found_list_num) {
				if (can_relative_metric >= can_opt_Heap_max(0).relative_metric);
				else {
					metric_point_opt can_metric_point(former_state, i - 1, list_proceed, can_relative_metric);
					can_opt_Heap_max.replace_top(can_metric_point);
					can_opt_Heap_min.push(can_metric_point);
					reduncency_in_Heap_min++;
				}
			}
			else {
				metric_point_opt can_metric_point(former_state, i - 1, list_proceed, can_relative_metric);
				can_opt_Heap_max.push(can_metric_point);
				can_opt_Heap_min.push(can_metric_point);
			}
		}

		int can_opt_Heap_max_last_ind = can_opt_Heap_max.size() - 1;
		int can_opt_Heap_min_last_ind = can_opt_Heap_min.size() - 1;

		// if the reduncency is less than 50 percent of the tree, keep it, only 5 percent of time reduced
		//if (reduncency_in_Heap_min < 64 || reduncency_in_Heap_max * 2 < to_be_found_list_num);
		//else {
		//	int pop_num = 0;
		//	int pop_ind = can_opt_Heap_max_last_ind;
		//	while (pop_num < reduncency_in_Heap_max) {
		//		if (can_opt_Heap_max(pop_ind) < can_opt_Heap_min(0)) {
		//			can_opt_Heap_max.switch_ele(pop_ind, can_opt_Heap_max_last_ind);
		//			can_opt_Heap_max.pop_back();
		//			can_opt_Heap_max_last_ind--;
		//			//can_opt_Heap_max.pop(pop_ind);
		//			pop_num++;
		//		}
		//		pop_ind--;
		//	}
		//	can_opt_Heap_max.build();
		//	reduncency_in_Heap_max = 0;
		//}

		//if (reduncency_in_Heap_min < 64 || reduncency_in_Heap_min * 2 < to_be_found_list_num);
		//else{
		//	int pop_num = 0;
		//	int pop_ind = can_opt_Heap_min_last_ind;
		//	while (pop_num < reduncency_in_Heap_min) {
		//		if (can_opt_Heap_min(pop_ind) > can_opt_Heap_max(0)) {
		//			can_opt_Heap_min.switch_ele(pop_ind, can_opt_Heap_min_last_ind);
		//			can_opt_Heap_min.pop_back(); 
		//			can_opt_Heap_min_last_ind--;
		//			pop_num++;
		//		}
		//		pop_ind--;
		//	}
		//	can_opt_Heap_min.build();
		//	reduncency_in_Heap_min = 0;
		//}

		// now you get the second opt metric path, backtracking it
		// the main complexity for list viterbi
		list_proceed++;
		//cout << "list_proceed=" << list_proceed << endl;
		diverge_time = can_opt_Heap_min(0).diverge_time;
		diverge_time++;		// adding 1 temporary, this become merge bit,for less computation only

		//cout << "can_opt.top()=" << can_opt.top() << endl;

		int diverge_list_ind = can_opt_Heap_min(0).list_belonging;
		// copy the path before diverge
		for (int i = n - 1; i > diverge_time; --i) {
			list_v_hat(list_proceed, i) = list_v_hat(diverge_list_ind, i);
		}
		// flip the merge bit
		list_v_hat(list_proceed, diverge_time) = list_v_hat(diverge_list_ind, diverge_time).flip();
		relative_metric(list_proceed) = can_opt_Heap_min(0).relative_metric;		// right

		diverge_time--;		// minus 1, change the variable back
		int track_state = can_opt_Heap_min(0).diverge_state;
		for (int i = diverge_time; i >= 0; --i) {		// path back tracking from time instance n-1 to 0
			last_opt_state(i) = track_state;
			track_state = state(i, track_state).past;		// former state
			if (time_instance_state(i) == 1 || time_instance_state(i) == 3) {
				// special case, time_instance_state(i) == 3
				// note that under this case, former sate size == present state size

				// or, time_instance_state(i) == 1
				// the decoded bit is determined only by the present state, i.e. last_opt_state(i)
				list_v_hat(list_proceed, i) = (last_opt_state(i) & 1) == 1;
			}
			else {
				// case of contention or shrinking
				list_v_hat(list_proceed, i) = track_state != last_opt_state(i);	// former state changed, decoded bit is 1
			}
		}

		// pop the minimum
		can_opt_Heap_min.pop();
		reduncency_in_Heap_max++;

		//cout << "can_opt_Heap_max: after" << can_opt_Heap_max;
		//cout << "can_opt_Heap_max.if_heap_gt() = " << can_opt_Heap_max.if_heap_gt() << endl;
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
	Matrix<GF2> decode_v(const Matrix<T>& r_or_hdr, int max_list_num = 1, int valid_list_num = 1) {
		// need some extension, for divergence and shrinking procees with contention

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
		state(0, PM_col_reps(0)) = past_and_metric(0, dist(1, 0));

		for (int i = 1; i < n; ++i) {

			int former_state_size = state.col(i - 1);
			int present_state_size = state.col(i);

			switch (time_instance_state(i)) {
			case(1):
				// for the divergence instance, not all state is valid, computation can be saved
				for (int j = 0; j < former_state_size; j++) {
					// j stands for the former states, scan over the former state in this for loop
					// if the first r column of PM are independent, actually no contention happened

					// update correlation distance and past state
					// choose the best path for each merging state

					my_double metric_store = state(i - 1, j).metric;
					state(i, j << 1) = past_and_metric(j, metric_store + dist(0, i));
					state(i, (j << 1) ^ PM_col_reps(i)) = past_and_metric(j, metric_store + dist(1, i));
				}
				break;
			case(2):
				// contention only, neither divergence nor shrinking, j is the present state index
				for (int j = 0; j < present_state_size; j++) {
					// update correlation distance and past state
					my_double can_1 = state(i - 1, j).metric + dist(0, i);	// accumulate correlation distance for each path to state j
					my_double can_2 = state(i - 1, j ^ PM_col_reps(i)).metric + dist(1, i);	// key complexity

					// choose the best path for each merging state
					state(i, j) = can_1 < can_2 ? past_and_metric(j, can_1) : past_and_metric(j ^ PM_col_reps(i), can_2);
				}
				break;
			case(3):
				// no contention, both divergence and shrinking, j is the former state index, special case (rare)
			{
				int half_former_state_size = former_state_size >> 1;
				for (int j = 0; j < half_former_state_size; j++) {
					// j stands for the former states, scan over the former state in this for loop

					// update correlation distance and past state
					// choose the valid path for each merging state

					my_double metric_store = state(i - 1, j).metric;
					state(i, j << 1) = past_and_metric(j, metric_store + dist(0, i));
				}
				for (int j = half_former_state_size; j < former_state_size; j++) {

					// choose the valid path for each merging state

					my_double metric_store = state(i - 1, j).metric;
					state(i, ((j ^ half_former_state_size) << 1) ^ PM_col_reps(i)) = past_and_metric(j, metric_store + dist(1, i));
				}
				break;
			}

			case(4):
				// contention and shrinking
				// j is the present state
				for (int j = 0; j < present_state_size; j++) {
					// consider one state in the for lop, i.e. state j

					// update correlation distance and past state
					my_double can_1 = state(i - 1, j).metric + dist(0, i);	// accumulate correlation distance for each path to state j
					my_double can_2 = state(i - 1, j ^ (PM_col_reps(i) + present_state_size)).metric + dist(1, i);

					// choose the best path for each merging state
					state(i, j) = can_1 < can_2 ? past_and_metric(j, can_1) : \
						past_and_metric(j ^ (PM_col_reps(i) + present_state_size), can_2);
				}
				break;
			}
		}
		//cout << "state" << state;

		list_v_hat.resize(max_list_num, n, false);
		relative_metric.resize(1, max_list_num, false);
		valid_list_ind.resize(1, valid_list_num, false);

		// set for listing

		// get the optimum decoded codeword, final state is 0
		int track_state = 0;					// track the current state, for final state, it must be 0
		for (int i = n - 1; i >= 0; --i) {		// path back tracking from time instance n-1 to 0
			last_opt_state(i) = track_state;
			track_state = state(i, track_state).past;		// former state
			if (time_instance_state(i) == 2 || time_instance_state(i) == 4) {
				list_v_hat(0, i) = track_state != last_opt_state(i);	// former state changed, decoded bit is 1
			}
			else {
				// special case, time_instance_state(i) == 3
				// note that under this case, former sate size == present state size

				// or, time_instance_state(i) == 1
				// the decoded bit is determined only by the present state, i.e. last_opt_state(i)
				list_v_hat(0, i) = (last_opt_state(i) & 1) == 1;
			}
		}
		//cout << "list_v_hat" << list_v_hat;

		valid_list_ind.resize(1, 0, false);		// not releasing the space, just seting row and column
		bool is_unused_PM_exist = unused_PM.size() != 0;

		/* since the parity matrix is determined, we if unused_PM exists, it must be necessary */
		if (!is_unused_PM_exist || (is_unused_PM_exist && unused_PM.check_inner_product(list_v_hat.get_row(0)))) {
			valid_list_ind.push_back(0);
		}


		int used_list_num = 1;
		if (valid_list_ind.size() < valid_list_num) {

			// using 'next_sub_opt'

			relative_metric(0) = 0;
			list_proceed = 0;
			diverge_time = n - 1;

			can_opt_Heap_max.resize(1, max_list_num + n, false);
			can_opt_Heap_min.resize(1, 2 * max_list_num, false);
			can_opt_Heap_max.resize(1, 0, false);
			can_opt_Heap_min.resize(1, 0, false);

			reduncency_in_Heap_max = 0;
			reduncency_in_Heap_min = 0;

			do {
				next_subopt_v(max_list_num - used_list_num);
				// to be optimized when eleminating invalid list num, a question unsolved

				/* since the parity matrix is determined, if unused_PM exists, it must be necessary */
				if (!is_unused_PM_exist || (is_unused_PM_exist && unused_PM.check_inner_product(list_v_hat.get_row(used_list_num)))) {
					valid_list_ind.push_back(used_list_num);
				}
				used_list_num++;
			} while (used_list_num < max_list_num && valid_list_ind.size() < valid_list_num);


			// using 'next_sub_opt_3'

			//relative_metric(0) = 0;
			//list_proceed = 0;
			//diverge_time = n - 1;
			//can_opt_set.clear();

			//do {
			//	next_subopt_v3(max_list_num - used_list_num);
			//	// to be optimized when eleminating invalid list num, a question unsolved

			//	/* since the parity matrix is determined, if unused_PM exists, it must be necessary */
			//	if (!is_unused_PM_exist || (is_unused_PM_exist && unused_PM.check_inner_product(list_v_hat.get_row(used_list_num)))) {
			//		valid_list_ind.push_back(used_list_num);
			//	}
			//	used_list_num++;
			//} while (used_list_num < max_list_num && valid_list_ind.size() < valid_list_num);
		}


#ifdef RUN_MSG		
		cout << "relative_metric" << relative_metric;
		//cout << "list_v_hat: all" << list_v_hat.get_part(0, 0, used_list_num - 1, -1);
		cout << "best metric = " << relative_metric(valid_list_ind(0)) + state(n - 1, 0).metric << endl;		// final state is zero
		cout << "valid_list_ind" << valid_list_ind;		// see the valid list index, judge if the list size is enough
#ifdef count_operation_number
#ifdef use_my_double

		double_ope_num_after = my_double_auxiliary_storage::operation_number;
		cout << "(Viterbi_optimized) double_ope_num = " << double_ope_num_after - double_ope_num_before << endl;
#endif // use_my_double

		GF2_ope_num_after = GF2_auxiliary_storage::operation_number;
		cout << "(Viterbi_optimized) GF2_ope_num = " << GF2_ope_num_after - GF2_ope_num_before << endl;
#endif // count_operation_number

		if (valid_list_ind.size() == 0) {
			cout << "no valid list found, under themax_list_num = " << max_list_num << endl;
		}
#endif // RUN_MSG

		return list_v_hat.get_rows(valid_list_ind);				// the list contains valid codewords only
	}

	// back up for advanced data structure to find next subopt vector
	
	/**
	 * . find the next sub-optimum codeword, with parameter set in private area
	 * using data sturcture of set, the red-black tree, but the result is wrong
	 */
	void next_subopt_v3(int to_be_found_list_num) {
		// we consider limit the candidate optimum metric point number up to the number of to-be-found list vector

		//cout << "min_contention_col_ind=" << min_contention_col_ind << endl;

		for (int i = diverge_time; i >= min_contention_col_ind; --i) {

			int los = last_opt_state(i);
			bool is_best_path_0 = list_v_hat(list_proceed, i) == 0;
			int former_state;
			switch (time_instance_state(i)) {
			case(2):
				former_state = is_best_path_0 ? (los ^ PM_col_reps(i)) : los;
				break;
			case(4):
				former_state = is_best_path_0 ? (los ^ (PM_col_reps(i) + state.col(i))) : los;
				break;
			default:
				continue;
			}

			// take the worse path, but we should judge whether the worse path exists
			my_double dist_tmp = is_best_path_0 ? dist(1, i) : dist(0, i);

			// add the optional path to candidate list
			my_double can_relative_metric = state(i - 1, former_state).metric + dist_tmp \
				- state(i, los).metric + relative_metric(list_proceed);

			if ((int)can_opt_set.size() >= to_be_found_list_num) {
				auto iter = can_opt_set.end();
				iter--;							// the last element of set
				if (iter->relative_metric <= can_relative_metric);	// the candidate is not added
				else {
					// remove the last element of can_opt_set, and add the candidate
					can_opt_set.erase(iter);
					can_opt_set.emplace(metric_point_opt(former_state, i - 1, list_proceed, can_relative_metric));
				}
			}
			else {
				can_opt_set.emplace(metric_point_opt(former_state, i - 1, list_proceed, can_relative_metric));
			}
		}

		auto iter = can_opt_set.begin();
		// now you get the second opt metric path, backtracking it
		list_proceed++;
		//cout << "list_proceed=" << list_proceed << endl;
		diverge_time = iter->diverge_time;
		diverge_time++;		// adding 1 temporary, this become merge bit,for less computation only
		//cout << "can_opt.top()=" << can_opt.top() << endl;

		int diverge_list_ind = iter->list_belonging;
		// copy the path before diverge
		for (int i = n - 1; i > diverge_time; --i) {
			list_v_hat(list_proceed, i) = list_v_hat(diverge_list_ind, i);
		}
		// flip the merge bit
		list_v_hat(list_proceed, diverge_time) = list_v_hat(diverge_list_ind, diverge_time).flip();
		relative_metric(list_proceed) = iter->relative_metric;		// right

		diverge_time--;		// minus 1, change the variable back
		int track_state = iter->diverge_state;
		for (int i = diverge_time; i >= 0; --i) {		// path back tracking from time instance n-1 to 0
			last_opt_state(i) = track_state;
			track_state = state(i, track_state).past;		// former state
			if (time_instance_state(i) == 1 || time_instance_state(i) == 3) {
				// special case, time_instance_state(i) == 3
				// note that under this case, former sate size == present state size

				// or, time_instance_state(i) == 1
				// the decoded bit is determined only by the present state, i.e. last_opt_state(i)
				list_v_hat(list_proceed, i) = (last_opt_state(i) & 1) == 1;
			}
			else {
				// case of contention or shrinking
				list_v_hat(list_proceed, i) = track_state != last_opt_state(i);	// former state changed, decoded bit is 1
			}

		}

		// pop the minimum
		can_opt_set.erase(iter);
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
	Matrix<GF2> decode_v3(const Matrix<T>& r_or_hdr, int max_list_num = 1, int valid_list_num = 1) {
		// need some extension, for divergence and shrinking procees with contention

#ifdef RUN_MSG
#ifdef count_operation_number
#ifdef use_my_double

		int float_ope_num_before = my_float_auxiliary_storage::operation_number;
		int float_ope_num_after;
#endif // use_my_double

		int GF2_ope_num_before = GF2_auxiliary_storage::operation_number;
		int GF2_ope_num_after;
#endif // count_operation_number
#endif // RUN_MSG

		// to merge the method of hard and soft decoding, if T is ordered structure, then use soft decoidng, else use hard deocding
		bool is_soft = !(T() < T());
		//cout << "is_soft = " << is_soft << endl;

		//cout << "n = " << n << endl;

		// compute distance for each symbol to 0 and 1
		for (int i = 0; i < n; ++i) {
			// not key complexity
			T ri = r_or_hdr(i);

			my_float d0 = (float)(is_soft ? (ri - 1) : (ri != 0));
			my_float d1 = (float)(is_soft ? (ri + 1) : (ri == 0));

			dist(0, i) = d0 * d0;
			dist(1, i) = d1 * d1;
		}
		//cout << "dist" << dist;

		//cout << "PM_col_reps" << PM_col_reps;

		// for time i, state[i] mapping from state number to state metric and past state pair
		state(0, 0) = past_and_metric(0, dist(0, 0));
		state(0, PM_col_reps(0)) = past_and_metric(0, dist(1, 0));

		for (int i = 1; i < n; ++i) {

			int former_state_size = state.col(i - 1);
			int present_state_size = state.col(i);

			switch (time_instance_state(i)) {
			case(1):
				// for the divergence instance, not all state is valid, computation can be saved
				for (int j = 0; j < former_state_size; j++) {
					// j stands for the former states, scan over the former state in this for loop
					// if the first r column of PM are independent, actually no contention happened

					// update correlation distance and past state
					// choose the best path for each merging state

					my_float metric_store = state(i - 1, j).metric;
					state(i, j << 1) = past_and_metric(j, metric_store + dist(0, i));
					state(i, (j << 1) ^ PM_col_reps(i)) = past_and_metric(j, metric_store + dist(1, i));
				}
				break;
			case(2):
				// contention only, neither divergence nor shrinking, j is the present state index
				for (int j = 0; j < present_state_size; j++) {
					// update correlation distance and past state
					my_float can_1 = state(i - 1, j).metric + dist(0, i);	// accumulate correlation distance for each path to state j
					my_float can_2 = state(i - 1, j ^ PM_col_reps(i)).metric + dist(1, i);	// key complexity

					// choose the best path for each merging state
					state(i, j) = can_1 < can_2 ? past_and_metric(j, can_1) : past_and_metric(j ^ PM_col_reps(i), can_2);
				}
				break;
			case(3):
				// no contention, both divergence and shrinking, j is the former state index, special case (rare)
			{
				int half_former_state_size = former_state_size >> 1;
				for (int j = 0; j < half_former_state_size; j++) {
					// j stands for the former states, scan over the former state in this for loop

					// update correlation distance and past state
					// choose the valid path for each merging state

					my_float metric_store = state(i - 1, j).metric;
					state(i, j << 1) = past_and_metric(j, metric_store + dist(0, i));
				}
				for (int j = half_former_state_size; j < former_state_size; j++) {

					// choose the valid path for each merging state

					my_float metric_store = state(i - 1, j).metric;
					state(i, ((j ^ half_former_state_size) << 1) ^ PM_col_reps(i)) = past_and_metric(j, metric_store + dist(1, i));
				}
				break;
			}

			case(4):
				// contention and shrinking
				// j is the present state
				for (int j = 0; j < present_state_size; j++) {
					// consider one state in the for lop, i.e. state j

					// update correlation distance and past state
					my_float can_1 = state(i - 1, j).metric + dist(0, i);	// accumulate correlation distance for each path to state j
					my_float can_2 = state(i - 1, j ^ (PM_col_reps(i) + present_state_size)).metric + dist(1, i);

					// choose the best path for each merging state
					state(i, j) = can_1 < can_2 ? past_and_metric(j, can_1) : \
						past_and_metric(j ^ (PM_col_reps(i) + present_state_size), can_2);
				}
				break;
			}
		}
		//cout << "state" << state;

		list_v_hat.resize(max_list_num, n, false);
		relative_metric.resize(1, max_list_num, false);
		valid_list_ind.resize(1, valid_list_num, false);

		// set for listing

		// get the optimum decoded codeword, final state is 0
		int track_state = 0;					// track the current state, for final state, it must be 0
		for (int i = n - 1; i >= 0; --i) {		// path back tracking from time instance n-1 to 0
			last_opt_state(i) = track_state;
			track_state = state(i, track_state).past;		// former state
			if (time_instance_state(i) == 2 || time_instance_state(i) == 4) {
				list_v_hat(0, i) = track_state != last_opt_state(i);	// former state changed, decoded bit is 1
			}
			else {
				// special case, time_instance_state(i) == 3
				// note that under this case, former sate size == present state size

				// or, time_instance_state(i) == 1
				// the decoded bit is determined only by the present state, i.e. last_opt_state(i)
				list_v_hat(0, i) = (last_opt_state(i) & 1) == 1;
			}
		}
		//cout << "list_v_hat" << list_v_hat;

		valid_list_ind.resize(1, 0);		// not releasing the space, just seting row and column
		bool is_unused_PM_exist = unused_PM.size() != 0;

		/* since the parity matrix is determined, we if unused_PM exists, it must be necessary */
		if (!is_unused_PM_exist || (is_unused_PM_exist && unused_PM.check_inner_product(list_v_hat.get_row(0)))) {
			valid_list_ind.push_back(0);
		}

		int used_list_num = 1;
		if (valid_list_ind.size() < valid_list_num) {

			relative_metric(0) = 0;
			//can_opt = priority_queue<metric_point>();
			list_proceed = 0;
			diverge_time = n - 1;
			can_opt_set.clear();

			do {
				next_subopt_v3(max_list_num - used_list_num);
				// to be optimized when eleminating invalid list num, a question unsolved

				/* since the parity matrix is determined, if unused_PM exists, it must be necessary */
				if (!is_unused_PM_exist || (is_unused_PM_exist && unused_PM.check_inner_product(list_v_hat.get_row(used_list_num)))) {
					valid_list_ind.push_back(used_list_num);
				}
				used_list_num++;
			} while (used_list_num < max_list_num && valid_list_ind.size() < valid_list_num);
		}


#ifdef RUN_MSG		
		cout << "relative_metric" << relative_metric;
		//cout << "list_v_hat: all" << list_v_hat.get_part(0, 0, used_list_num - 1, -1);
		cout << "best metric = " << relative_metric(valid_list_ind(0)) + state(n - 1, 0).metric << endl;		// final state is zero
		cout << "valid_list_ind" << valid_list_ind;		// see the valid list index, judge if the list size is enough
#ifdef count_operation_number
#ifdef use_my_double

		float_ope_num_after = my_float_auxiliary_storage::operation_number;
		cout << "(Viterbi_optimized) float_ope_num = " << float_ope_num_after - float_ope_num_before << endl;
#endif // use_my_double

		GF2_ope_num_after = GF2_auxiliary_storage::operation_number;
		cout << "(Viterbi_optimized) GF2_ope_num = " << GF2_ope_num_after - GF2_ope_num_before << endl;
#endif // count_operation_number

		if (valid_list_ind.size() == 0) {
			cout << "no valid list found, under themax_list_num = " << max_list_num << endl;
		}
#endif // RUN_MSG

		return list_v_hat.get_rows(valid_list_ind);				// the list contains valid codewords only
	}

};