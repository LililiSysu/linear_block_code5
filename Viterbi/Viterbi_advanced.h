/*****************************************************************//**
 * \file   Viterbi_advanced.h
 * \brief  using advanced data structure in Viterbi, which makes Viterbi algorithm be able 
 * to cope with any linear block code, while keep an optimized complexity, but the algorithm is
 * still not best choice for any linear block code, you should use it carefully and only in the
 * testing case, but not simulation/application case
 * 
 * \author 26259
 * \date   April 2023
 *********************************************************************/

#pragma once

#include"Viterbi_common.h"
#include<unordered_map>
#include<set>
using namespace std;

class LC_OSD_r;

template<int m, int t>
class LL_OSD_binary_H_hybrid_Viterbi_v2;

/**
 * .general viterbi algorithm for any linear block code, using the data structure of unordered map
 * the computation cost in this algorithm is minimal given a parity Matrix of a linear block code
 * but the storage is not optimal, due to the advanced data structure, which cost about 2 to 6 times 
 * as optimal data storage need, you can check the instruction of unordered map to know the detail
 * 
 * this is useful for a non-structured parity matrix, you can check the states number, operation number, etc,
 * in decoding with low latency OSD, this is used for avoiding Gaussian Elimination, and unused Parity check
 * Matrix can be added to the class to verify the correctness of the list codeword
 * 
 * i admit that this class is a little complicated, the interface can be simplified, if in the future we need
 */
class Viterbi_unordered_map {
private:
	int n;
	int r;
	Matrix<GF2> PM;
	Matrix<int> PM_permutation;
	int n_total;
	Matrix<my_double> zero_col_dist;
	Matrix<GF2> zero_col_hdr;
	int zero_col_hdr_int;

	Matrix<GF2> list_v_hat;
	Matrix<int> PM_col_reps;
	// the state correspond to j^{th} received bits, indicated by j^{th} column to a parity matrix PM, now j=0
	// state is indicated by int numner form by add and shifting of binary column of PM, e.g. [1,1,0]^T -> 6

	Matrix<my_double> dist;// initialize correlation distance, considering both case for soft metric and hard metric

	/* list decoding */
	Matrix<int> last_opt_state;				// store state at time instant [1,n-1], the initial state 0 is not stored
	Matrix<my_double> relative_metric;		// for each list, need a bias at merge point for metric unification
	priority_queue<metric_point> can_opt;	// candidate merge point for next optimum path
	set<metric_point_opt> can_opt_set;		// using set to reduce the size of can_opt
	my_double absolute_metric;
	my_double absolute_zero_col_mertric;

	int list_proceed;
	int diverge_time;
	int min_contention_col_ind;
	int max_state_num;

	// new variable for list viterbi
	priority_queue<ending_one> fixed_PM_row_end_one_ind;
	// consider using unordered_map for recording state metric
	// replace it by ordered vector, after state is determined, inquire will only take O(log(n)) complexity, to be finished.
	vector<unordered_map<int, past_and_metric>> state;		// unordered_map takes about 6 times storage than vector, which make it slow

	Matrix<char> time_seperation;			// the i^th element indicates the state changing at time i, 
	// 0 means divergence with no contention, 1 means no divergence and no contention, 2 means no divergence and contention
	Matrix<char> no_contention_no_divergence_state;		// for no divergence and no contention time instance, 
	// indicates how to update each state, 0 means that update by path 0, 1 menas that updates by path 1, and 2 means that no updates

	Matrix<int> rank_distribution;
	Matrix<GF2> unused_PM;
	Matrix<int> valid_list_ind;

	Matrix<GF2> can_v_hat;
	Matrix<GF2> can_first_v_hat;
	Matrix<int> non_0_ind_aux;

	void initialize_size() {
		n_total = n;
		/* apply space */

		PM_col_reps = Matrix<int>(1, n);
		dist = Matrix<my_double>(2, n);

		/* list decoding */
		last_opt_state = Matrix<int>(1, n);

		// consider using unordered_map for recording state metric
		state = vector<unordered_map<int, past_and_metric>>(n);

		int old_size = rank_distribution.size();
		rank_distribution.resize(1, r + 1, true);					// confusing, but unimportant
		rank_distribution.set_part(0, old_size, 0, r, 0);

		can_v_hat.resize(1, n_total, false);
		can_first_v_hat.resize(1, n_total, false);
		non_0_ind_aux.resize(1, n_total, false);

		zero_col_hdr_int = 0;
		absolute_metric = 0;
		absolute_zero_col_mertric = 0;
		list_proceed = 0;
	}
	void initialize_PM() {
		n_total = n;
		PM_permutation = find_opt_PM::column_permute(PM);		// although this may not the best
		n = PM.col();		// warning: n changed

		//cout << "PM_preprocessed: inner" << PM;

		for (int i = 0; i < n; ++i) {
			PM_col_reps(i) = 0;
			for (int p = 0; p < r; ++p) {
				PM_col_reps(i) <<= 1;
				PM_col_reps(i) += (int)PM(p, i);
			}
		}
		// in may be erasing states when rank right-part of the PM decreases, 
		// for simplicity we only consider zeros ending that states
		fixed_PM_row_end_one_ind = priority_queue<ending_one>();		// claer out
		for (int i = 0; i < r; ++i) {
			int j = n - 1;
			while (j >= 0 && PM(i, j) == 0) {
				--j;
			}
			fixed_PM_row_end_one_ind.push(ending_one(i, j));
		}
		min_contention_col_ind = -1;

		if (n_total != n) {
			zero_col_dist.resize(1, 1 << (n_total - n), false);
			zero_col_hdr.resize(1, n_total - n, false);
		}
	}

public:
	int log2_max_state_num;
	Matrix<GF2> test_zero;

	Viterbi_unordered_map(const Matrix<GF2>& _PM) {
		// get the code parameters, for a linear block code of (n,k) and r = n - k, 
		// with parity check matrix PM

		n = _PM.col();
		r = _PM.row();
		initialize_size();

		PM = _PM;
		initialize_PM();

		// this is for repeatly use the same PM, initialize state at constructor
		//initialize_state();
	}
	Viterbi_unordered_map(int _n = 0, int _r = 0) {
		n = _n;
		r = _r;
		initialize_size();
	}

	void resize(int _n, int _r) {
		n = _n;
		r = _r;
		initialize_size();
	}

	void initialize_state() {
		// initialize the sturcture of state, this is for using the same parity check matrix repeatedly
		priority_queue<ending_one> PM_row_end_one_ind = fixed_PM_row_end_one_ind;

		int next_ending_one = PM_row_end_one_ind.top().one_pos_ind;
		int ending_row_ind = PM_row_end_one_ind.top().row_ind;
		time_seperation = Matrix<char>(1, n);
		time_seperation.reset(-1);
		no_contention_no_divergence_state = Matrix<char>(1, 0, 'v');

		// for time i, state[i] mapping from state number to state metric and past state pair
		state[0][0] = past_and_metric();	// if exist, means active

		max_state_num = 0;

		// this is only valid if ending_row_ind != 0
		if (next_ending_one != 0) {
			state[0][PM_col_reps(0)] = past_and_metric();
		}
		else {
			while (next_ending_one == 0) {	// clear all rows ending with last non-zero on position 0
				PM_row_end_one_ind.pop();
				if (!PM_row_end_one_ind.empty()) {
					next_ending_one = PM_row_end_one_ind.top().one_pos_ind;
					ending_row_ind = PM_row_end_one_ind.top().row_ind;
				}
				else {
					// this is the ending of decoding
					break;
				}
			}
		}

		max_state_num = max_state_num < (int) state[0].size() ? (int)state[0].size() : max_state_num;

		//Matrix<bool> is_explored(1, pow_2_r, '0');
		//is_explored(0) = true;
		//is_explored(PM_col_reps(0)) = true;

		for (int i = 1; i < n; ++i) {
			int test_new_state = PM_col_reps(i);
			int former = i - 1;

			// copy the former states
			state[i] = state[former];

			if (state[former].find(test_new_state) == state[former].end()) {
				// new state should be opened
				for (auto iter = state[former].begin(); iter != state[former].end(); ++iter) {
					// add the state
					int j = (iter->first) ^ test_new_state;

					// this make sure that state (iter->first) must in key of state[former]
					// if j in key of state[i], then j ^ test_new_state in key of state[former]

					//if (is_explored(j) == false)		// must be false
					state[i][j] = past_and_metric();
				}
			}

			while (i == next_ending_one) {
				// removing the ending state at receiver
				int big_mark = 1 << (r - ending_row_ind - 1);
				auto iter = state[i].begin();
				while (iter != state[i].end()) {
					if (((iter->first) & big_mark) == big_mark) {
						state[i].erase(iter++);
					}
					else {
						iter++;
					}
				}

				PM_row_end_one_ind.pop();
				if (!PM_row_end_one_ind.empty()) {
					next_ending_one = PM_row_end_one_ind.top().one_pos_ind;
					ending_row_ind = PM_row_end_one_ind.top().row_ind;
				}
				else {
					// this is the ending of decoding
					break;
				}
			}

			//cout << "---------- i= " << i << "----------" << endl;
			//for (auto iter = state[i].begin(); iter != state[i].end(); ++iter) {
			//	// current state, last state, and metric
			//	cout << iter->first << ": " << iter->second << endl;
			//}

			if (state[former].find(test_new_state) == state[former].end()) {
				if (state[i].find(test_new_state) != state[i].end()) {
					// divergence with no contentioin
					time_seperation(i) = 0;
					////cout << "i=" << i << ", divergence with no contentioin" << endl;

				}
				else {
					// no divergence, no contention
					time_seperation(i) = 1;
					//cout << "i=" << i << ", no contention" << endl;

					for (auto iter = state[former].begin(); iter != state[former].end(); ++iter) {
						// the following may be embarrasing, the condition of if and program switching may take time

						int j = iter->first;	// starting state

						auto update_state = state[i].find(j);
						if (update_state != state[i].end()) {		// path 0, make sure not creating new state
							
							//update_state->second = past_and_metric();
							no_contention_no_divergence_state.push_back(0);

						}
						else {		// at most one path will exist
							int new_state = j ^ test_new_state;
							update_state = state[i].find(new_state);
							if (update_state != state[i].end()) {		// path 1, whic may be not exists
								
								//update_state->second = past_and_metric();
								no_contention_no_divergence_state.push_back(1);

							}
							else {
								no_contention_no_divergence_state.push_back(2);
							}
						}

					}
				}
			}
			else {
				// contention, the new column is dependent of key of state[former], which forms a linear space in GF2
				time_seperation(i) = 2;

				//cout << "i=" << i << ", contention" << endl;
				min_contention_col_ind = min_contention_col_ind == -1 ? i : min_contention_col_ind;	// first contention ind
			}

			//cout << "---------- after ----------" << endl;
			//for (auto iter = state[i].begin(); iter != state[i].end(); ++iter) {
			//	// current state, last state, and metric
			//	cout << iter->first << ": " << iter->second << endl;
			//}

			max_state_num = max_state_num < (int)state[i].size() ? (int)state[i].size() : max_state_num;
		}
	}
	/**
	 * .change PM not changing the size of code, i.e., n, and k not changed.
	 * not initializing state, for decode_v_once especially.
	 */
	void change_PM(const Matrix<GF2>& _PM) {
		PM = _PM;
		n = PM.col();
		initialize_PM();
	}
	/**
	 * .change PM not changing the size of code, i.e., n, and k not changed.
	 * not initializing state, for decode_v_once especially.
	 */
	void change_unused_PM(const Matrix<GF2>& _unused_PM) {
		if (_unused_PM.size() == 0) {
			return;
		}
		unused_PM = _unused_PM;
		unused_PM.permute_col(PM_permutation);
	}

	void print_state() const{		
		for (int s = 0; s < max_state_num; ++s) {
			for (int i = 0; i < n; ++i) {
				if (state[i].find(s) != state[i].end()) {
					cout << "(";
					my::print_binary(s, my::_log2(max_state_num));
					cout << ")";
					cout << "\t";
				}
				else {
					cout << "\t";
				}
			}
			cout << endl;
		}
	}

	/**
	 * .apply decode method to 'r_or_hdr', under BPSK with 0 -> 1 and 1 -> -1
	 * optimized for saving the storage space, if use the parity matrix once, use this function
	 *
	 * \param r_or_hdr: received codeword, can be soft or hard
	 * \param list_size: number of decoded v listed by row, which are maximum likelihood
	 * \return decoded result, v
	 */
	template<class T>
	void decode_v_once_inner(Matrix<T> r_or_hdr, int list_num = 1) {
		// to merge the method of hard and soft decoding, if T is ordered structure, then use soft decoidng, else use hard deocding
		bool is_soft = !(T() < T());

		//cout << "is_soft = " << is_soft << endl;
		//cout << "n = " << n << endl;

		// throw away symbols w.r.t. zero columns, n_total -> n
		r_or_hdr.permute(PM_permutation);
		if (n == n_total);
		else {
			// make a table to sotre all possible of zero column position
			zero_col_dist.reset(0);
			int zero_position_num = n_total - n;
			int zero_position_possible = 1 << zero_position_num;

			for (int j = 0; j < zero_position_num; ++j) {
				// not key complexity
				T ri = r_or_hdr(n_total - 1 - j);
				my_double d0 = (double)(is_soft ? (ri - 1) : (ri != 0));
				my_double d1 = (double)(is_soft ? (ri + 1) : (ri == 0));
				d0 *= d0;
				d1 *= d1;

				int i = 0;
				int pow_2_j = 1 << j;
				while (i < zero_position_possible) {

					for (int k = 0; k < pow_2_j; ++k) {
						zero_col_dist(i + k) += d0;
						zero_col_dist(i + pow_2_j + k) += d1;
					}
					i += pow_2_j * 2;
				}
				zero_col_hdr(zero_position_num - 1 - j) = d0 > d1;
			}

			/*cout << "r_or_hdr" << r_or_hdr;
			cout << "zero_col_dist" << zero_col_dist;
			cout << "zero_col_hdr" << zero_col_hdr;*/

			zero_col_hdr_int = 0;
			for (int i = 0; i < n_total - n; ++i) {
				zero_col_hdr_int <<= 1;
				zero_col_hdr_int += int(zero_col_hdr(i));
			}

			absolute_zero_col_mertric = zero_col_dist(zero_col_hdr_int);
			for (int i = 0; i < zero_position_possible; ++i) {
				zero_col_dist(i) -= absolute_zero_col_mertric;
			}

			r_or_hdr = r_or_hdr.get_part(0, 0, 0, n - 1);
		}

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

		priority_queue<ending_one> PM_row_end_one_ind = fixed_PM_row_end_one_ind;
		int next_ending_one = PM_row_end_one_ind.top().one_pos_ind;
		int ending_row_ind = PM_row_end_one_ind.top().row_ind;

		// for i==0, the first received symbol (to be discarded)
		//path_metric(0, 0) = dist(0, 0);
		//past_state(0, 0) = 0;				// initial and final state is 0 for a linear block code, set the past state be 0

		// initialize past state
		//path_metric(PM_col_reps(0), 0) = dist(1, 0);
		//past_state(PM_col_reps(0), 0) = 0;		// second path derived from initial 0 state

		//ope_num_after = my_double_auxiliary_storage::operation_number;
		//cout << "ope_num = " << ope_num_after - ope_num_before << endl;

		//cout << "PM_col_reps" << PM_col_reps;

		max_state_num = 0;

		state[0].clear();
		// for time i, state[i] mapping from state number to state metric and past state pair
		state[0][0] = past_and_metric(0, dist(0, 0));	// if exist, means active

		// this is only valid if ending_row_ind != 0
		if (PM_col_reps(0) != 0) {
			if (next_ending_one != 0) {
				state[0][PM_col_reps(0)] = past_and_metric(0, dist(1, 0));
			}
			else {
				while (next_ending_one == 0) {	// clear all rows ending with last non-zero on position 0
					PM_row_end_one_ind.pop();
					if (!PM_row_end_one_ind.empty()) {
						next_ending_one = PM_row_end_one_ind.top().one_pos_ind;
						ending_row_ind = PM_row_end_one_ind.top().row_ind;
					}
					else {
						// this is the ending of decoding
						break;
					}
				}
			}
		}
		else {
			// contention at state[0], this is impossible since zero columns are permuted to the back
			min_contention_col_ind = 0;
			state[0][0].metric = dist(0, 0) < dist(0, 1) ? dist(0, 0) : dist(0, 1);
		}

		max_state_num = max_state_num < (int)state[0].size() ? (int)state[0].size() : max_state_num;

		//Matrix<bool> is_explored(1, pow_2_r, '0');
		//is_explored(0) = true;
		//is_explored(PM_col_reps(0)) = true;

		for (int i = 1; i < n; ++i) {
			int test_new_state = PM_col_reps(i);
			int former = i - 1;

			// copy the former states
			state[i] = state[former];		// this clean out the state, from begin to end

			if (state[former].find(test_new_state) == state[former].end()) {
				// new state should be opened
				for (auto iter = state[former].begin(); iter != state[former].end(); ++iter) {
					// add the state
					int j = (iter->first) ^ test_new_state;

					// this make sure that state (iter->first) must in key of state[former]
					// if j in key of state[i], then j ^ test_new_state in key of state[former]

					//if (is_explored(j) == false)		// must be false
					state[i][j] = past_and_metric(-invalid_large_metric, invalid_large_metric);
				}
			}

			while (i == next_ending_one) {
				// removing the ending state at receiver
				int big_mark = 1 << (r - ending_row_ind - 1);
				auto iter = state[i].begin();
				while (iter != state[i].end()) {
					if (((iter->first) & big_mark) == big_mark) {
						state[i].erase(iter++);
					}
					else {
						iter++;
					}
				}

				PM_row_end_one_ind.pop();
				if (!PM_row_end_one_ind.empty()) {
					next_ending_one = PM_row_end_one_ind.top().one_pos_ind;
					ending_row_ind = PM_row_end_one_ind.top().row_ind;
				}
				else {
					// this is the ending of decoding
					break;
				}
			}

			max_state_num = max_state_num < (int)state[i].size() ? (int)state[i].size() : max_state_num;

			//cout << "---------- i= " << i << "----------" << endl;
			//for (auto iter = state[i].begin(); iter != state[i].end(); ++iter) {
			//	// current state, last state, and metric
			//	cout << iter->first << ": " << iter->second << endl;
			//}

			if (state[former].find(test_new_state) == state[former].end()) {
				if (state[i].find(test_new_state) != state[i].end()) {
					// divergence with no contentioin
					//cout << "i=" << i << ", divergence with no contentioin" << endl;

					for (auto iter = state[former].begin(); iter != state[former].end(); ++iter) {
						int j = iter->first;	// starting state
						my_double last_metric = iter->second.metric;
						state[i][j] = past_and_metric(j, last_metric + dist(0, i));

						/*path_metric(j, i) = path_metric(j, i - 1) + dist(0, i);
						past_state(j, i) = j;*/

						int new_state = j ^ test_new_state;
						//is_explored(new_state) = true;

						state[i][new_state] = past_and_metric(j, last_metric + dist(1, i));

						/*path_metric(new_state, i) = path_metric(j, i - 1) + dist(1, i);
						past_state(new_state, i) = j;*/
					}
				}
				else {
					// no divergence, no contention, this is very rare case
					//cout << "i=" << i << ", no contention" << endl;

					for (auto iter = state[former].begin(); iter != state[former].end(); ++iter) {
						// the following may be embarrasing, the condition of if and program switching may take time

						int j = iter->first;	// starting state
						my_double last_metric = iter->second.metric;

						auto update_state = state[i].find(j);
						if (update_state != state[i].end()) {		// path 0, make sure not creating new state
							update_state->second = past_and_metric(j, last_metric + dist(0, i));

							//path_metric(j, i) = path_metric(j, i - 1) + dist(0, i);
							//past_state(j, i) = j;
						}
						else {		// at most one path will exist
							int new_state = j ^ test_new_state;
							update_state = state[i].find(new_state);
							if (update_state != state[i].end()) {		// path 1, whic may be not exists
								update_state->second = past_and_metric(j, last_metric + dist(1, i));

								/*path_metric(new_state, i) = path_metric(j, i - 1) + dist(1, i);
								past_state(new_state, i) = j;*/
							}
						}

					}
				}
			}
			else {
				// contention, the new column is dependent of key of state[former], which forms a linear space in GF2
				//cout << "i=" << i << ", contention" << endl;

				min_contention_col_ind = min_contention_col_ind == -1 ? i : min_contention_col_ind;	// first contention ind
				for (auto iter = state[i].begin(); iter != state[i].end(); ++iter) {
					int j = iter->first;	// ending state
					int new_state = j ^ test_new_state;

					// if j in state[i], j must be in key of state[former], since it is not erased by ending state
					// the key of state[former] always forms a complete linear space of GF2.
					// also j ^ test_new_state in key of state[former], as we analyzed above

					// update correlation distance and past state
					my_double can_1 = state[former][j].metric + dist(0, i);
					// accumulate correlation distance for each path to state j
					my_double can_2 = state[former][new_state].metric + dist(1, i);	// key complexity

					// choose the best path for each merging state
					state[i][j] = can_1 < can_2 ? past_and_metric(j, can_1) : past_and_metric(new_state, can_2);

					/*bool is_can_1_better = can_1 < can_2;
					path_metric(j, i) = is_can_1_better ? can_1 : can_2;
					past_state(j, i) = is_can_1_better ? j : (j ^ test_new_state);*/
				}
			}

			//cout << "---------- after ----------" << endl;
			//for (auto iter = state[i].begin(); iter != state[i].end(); ++iter) {
			//	// current state, last state, and metric
			//	cout << iter->first << ": " << iter->second << endl;
			//}
		}
		//min_contention_col_ind = min_contention_col_ind == -1 ? 1 : min_contention_col_ind;

		// counting rank distribution
		log2_max_state_num = my::_log2(max_state_num);
		//cout << "log2_max_state_num = " << log2_max_state_num << endl;
		rank_distribution(log2_max_state_num)++;
		//rank_distribution(log2_max_state_num - 1)++;


		//cout << "min_contention_col_ind=" << min_contention_col_ind << endl;

		// memory check
		//cout << "state" << endl;
		//for (int i = 0; i < n; ++i) {
		//	cout << "----- i=" << i << " -----" << endl;
		//	for (auto iter = state[i].begin(); iter != state[i].end(); ++iter) {
		//		// current state, last state, and metric
		//		cout << iter->first << ": " << iter->second <<";";
		//	}
		//	cout << endl << "state[" << i << "].size() = " << state[i].size() << endl;
		//}

		/*cout << "path_metric" << path_metric.get_part(0, 0, 3, n - 1);
		cout << "past_state" << past_state.get_part(0, 0, 3, n - 1);*/

		list_v_hat.resize(list_num, n_total, false);

		// set for listing

		// get the optimum decoded codeword, final state is 0
		int track_state = 0;					// track the current state, for final state, it must be 0
		for (int i = n - 1; i >= 0; --i) {		// path back tracking from time instance n-1 to 0
			last_opt_state(i) = track_state;
			track_state = state[i][track_state].past;
			list_v_hat(0, i) = track_state != last_opt_state(i);

			//list_v_hat(0, i) = past_state(track_state, i) != track_state;	// if state changed, the decoded bit is 1 else it is 0			
			//track_state = past_state(track_state, i);						// update the current state a time instant earlier
		}

		// add the hard decision of zero columns' position at the end of list_v_hat
		for (int i = n; i < n_total; ++i) {
			list_v_hat(0, i) = zero_col_hdr(i - n);
		}

		absolute_metric = state[n - 1][0].metric;
	}

	/**
	 * .apply decode method to 'r_or_hdr', under BPSK with 0 -> 1 and 1 -> -1
	 * optimized for saving the storage space, if use the parity matrix once, use this function
	 *
	 * \param r_or_hdr: received codeword, can be soft or hard
	 * \param list_size: number of decoded v listed by row, which are maximum likelihood
	 * \return decoded result, v
	 */
	template<class T>
	Matrix<GF2> decode_v_once(Matrix<T> r_or_hdr, int list_num = 1) {	

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

		decode_v_once_inner(r_or_hdr, list_num);

		relative_metric.resize(1, list_num, false);	
		valid_list_ind.resize(1, 0, false);		// not releasing the space, just seting row and column

		/* assume that if PM is full rank, then there is no need to use unused_PM */
		can_v_hat = list_v_hat.get_row(0);

		can_first_v_hat.resize(1, 0, false);
		if (unused_PM.check_inner_product_4_GF2(can_v_hat, can_first_v_hat, non_0_ind_aux)) {
			can_first_v_hat = can_v_hat;
			valid_list_ind.push_back(0);
		}

		int used_list_num = 1;
		if (list_num > 1) {
			// using next_subopt_v()
			
			//relative_metric(0) = 0;
			//can_opt = priority_queue<metric_point>();
			//list_proceed = 0;
			//diverge_time = n - 1;
			//for (int i = 1; i < list_num; ++i) {
			//	next_subopt_v();
			//	/* assume that if PM is full rank, then there is no need to use unused_PM*/
			//	if (unused_PM.check_inner_product_4_GF2(list_v_hat.get_row(i), can_first_v_hat, non_0_ind_aux)) {
			//		if (can_first_v_hat.col() == 0) {
			//			can_first_v_hat = list_v_hat.get_row(i);
			//		}
			//		valid_list_ind.push_back(i);
			//	}
			//}

			// using nex_subopt_v_set()

			relative_metric(0) = 0;
			//can_opt = priority_queue<metric_point>();
			list_proceed = 0;
			diverge_time = n - 1;
			can_opt_set.clear();
			add_zero_col_pos_2_can_opt_set(list_num - 1);

			do {
				//cout << "list_proceed = " << list_proceed << endl;
				next_subopt_v_set(list_num - used_list_num);
				/* assume that if PM is full rank, then there is no need to use unused_PM*/
				can_v_hat = list_v_hat.get_row(used_list_num);
				if (unused_PM.check_inner_product_4_GF2(can_v_hat, can_first_v_hat, non_0_ind_aux)) {

					if (can_first_v_hat.col() == 0) {
						can_first_v_hat = can_v_hat;
					}
					valid_list_ind.push_back(used_list_num);
				}
				used_list_num++;
			} while (used_list_num < list_num);
		}
		
		list_v_hat = list_v_hat.get_rows(valid_list_ind);
		//cout << "relative_metric" << relative_metric;

		list_v_hat.permute_col_back(PM_permutation);

		//cout << "zero_col_dist" << zero_col_dist;
		//cout << "relative_metric" << relative_metric;
		//cout << "absolute_metric = " << absolute_metric << endl;

#ifdef RUN_MSG
		//cout << "relative_metric" << relative_metric;
		/*cout << "can_opt.size() = " << can_opt.size() << endl;
		while (!can_opt.empty()) {
			cout << can_opt.top();
			can_opt.pop();
		}
		cout << endl;*/

		//cout << "list_v_hat" << list_v_hat;
		cout << "valid_list_ind" << valid_list_ind;
		cout << "relative_metric" << relative_metric;

#ifdef count_operation_number
#ifdef use_my_double

		double_ope_num_after = my_double_auxiliary_storage::operation_number;
		cout << "(Viterbi_unordered_map) double_ope_num = " << double_ope_num_after - double_ope_num_before << endl;
#endif // use_my_double

		GF2_ope_num_after = GF2_auxiliary_storage::operation_number;
		cout << "(Viterbi_unordered_map) GF2_ope_num = " << GF2_ope_num_after - GF2_ope_num_before << endl;
#endif // count_operation_number

		/*if (valid_list_ind(0) != 0) {
			cout << "select list change the result" << endl;
		}*/
		//cout << "valid_list_ind" << valid_list_ind;
#endif // RUN_MSG
				
		return list_v_hat;
	}

	template<class T>
	Matrix<GF2> decode_v_4_LC_OSD_once(Matrix<T> recv, LC_OSD_r& lc_osd, int max_list_num);

	template<class T, class LL_OSD_type>
	Matrix<GF2> decode_v_4_LL_OSD_once(Matrix<T> recv, LL_OSD_type& ll_osd_Viterbi, int max_list_num);

	template<class T, class LL_OSD_Hybrid_type>
	Matrix<GF2> decode_v_4_LL_OSD_Hybrid_once(Matrix<T> recv, LL_OSD_Hybrid_type& ll_osd_hybrid_Viterbi, int max_list_num);

	template<int m,int c>
	friend class LL_OSD_binary_H_hybrid_Viterbi_v2;

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
			//if (pointer_pm->list_belonging != -1)
			//	continue;		// merge into proceed state, do not need to continue

			/*if (finished_state(last_opt_state(i), i) != -1)
				break;*/

			int former = i - 1;
			bool is_best_path_0 = list_v_hat(list_proceed, i) == 0;
			// take the worse path, but we should judge whether the worse path exists
			int former_state = is_best_path_0 ? (los ^ PM_col_reps(i)) : los;
			my_double dist_tmp = is_best_path_0 ? dist(1, i) : dist(0, i);

			if (former != -1) {
				auto iter = state[former].find(former_state);
				if (iter != state[former].end()) {

					can_opt.push(metric_point(former_state, former, iter->second.metric + dist_tmp \
						- state[i][los].metric + relative_metric(list_proceed), los, list_proceed));

					/*can_opt.push(metric_point(former_state, i - 1, path_metric(former_state, i - 1) + dist_tmp \
						- path_metric(last_opt_state(i), i) + relative_metric(list_proceed), last_opt_state(i)));*/

				}
			}
			else {
				can_opt.push(metric_point(former_state, former, 0 + dist_tmp \
					- state[i][los].metric + relative_metric(list_proceed), los, list_proceed));
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
		//cout << " can_opt.top() = " << can_opt.top() << endl;

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
			track_state = state[i][track_state].past;							
			// update the current state a time instant earlier
			list_v_hat(list_proceed, i) = track_state != last_opt_state(i);		
			// if state changed, the decoded bit is 1, else 0

			/*last_opt_state(i) = track_state;
			list_v_hat(list_proceed, i) = past_state(track_state, i) != track_state;
			track_state = past_state(track_state, i);	*/
		}
		//cout << "list_v_hat" << list_v_hat;

		can_opt.pop();
		//cout << "second_opt.empty() = " << can_opt.empty() << endl;
	}

	/**
	 * . only call this for a new viterbi list, but not during scan of zero col pos
	 */
	void add_zero_col_pos_2_can_opt_set(int to_be_found_list_num) {

		bool is_set_iter_valid = false;
		set<metric_point_opt>::iterator set_iter;

		int zs = zero_col_dist.size();
		for (int i = 0; i < zs; ++i) {
			if (i != zero_col_hdr_int) {
				int former_state = 0;
				int diverge_time_fake = -i - 10;			// remember that this indicates zero columns' bits
				my_double can_relative_metric = zero_col_dist(i) + relative_metric(list_proceed);

				if ((int)can_opt_set.size() >= to_be_found_list_num) {
					if (is_set_iter_valid);
					else {
						set_iter = can_opt_set.end();
						set_iter--;							// the last element of set
						is_set_iter_valid = true;
					}
					if (set_iter->relative_metric <= can_relative_metric);	// the candidate is not added
					else {
						// remove the last element of can_opt_set, and add the candidate
						can_opt_set.erase(set_iter);
						can_opt_set.emplace(metric_point_opt(former_state, diverge_time_fake, list_proceed, can_relative_metric));

						// update the max element in set
						set_iter = can_opt_set.end();
						set_iter--;
					}
				}
				else {
					can_opt_set.emplace(metric_point_opt(former_state, diverge_time_fake, list_proceed, can_relative_metric));
				}
			}
		}
	}

	/**
	 * . find the next sub-optimum codeword, with parameter set in private area
	 */
	void next_subopt_v_set(int to_be_found_list_num) {
		// we consider limit the candidate optimum metric point number up to the number of to-be-found list vector

		//cout << "r=" << r << endl;
		//cout << "min_contention_col_ind=" << min_contention_col_ind << endl;
		// we mush have min_contention_col_ind>=1

		bool is_set_iter_valid = false;
		set<metric_point_opt>::iterator set_iter;

		// consider the space saving later, i.e., we donot need to store the first r state, since no contention
		for (int i = diverge_time; i >= min_contention_col_ind; --i) {

			int los = last_opt_state(i);
			bool is_best_path_0 = list_v_hat(list_proceed, i) == 0;

			// take the worse path, but we should judge whether the worse path exists
			int former_state = is_best_path_0 ? (los ^ PM_col_reps(i)) : los;
			my_double dist_tmp = is_best_path_0 ? dist(1, i) : dist(0, i);

			auto iter = state[i - 1].find(former_state);
			if (iter != state[i - 1].end()) {

				my_double can_relative_metric = iter->second.metric + dist_tmp \
					- state[i][los].metric + relative_metric(list_proceed);

				if ((int)can_opt_set.size() >= to_be_found_list_num) {
					if (is_set_iter_valid);
					else {
						set_iter = can_opt_set.end();
						set_iter--;							// the last element of set
						is_set_iter_valid = true;
					}
					if (set_iter->relative_metric <= can_relative_metric);	// the candidate is not added
					else {
						// remove the last element of can_opt_set, and add the candidate
						can_opt_set.erase(set_iter);
						can_opt_set.emplace(metric_point_opt(former_state, i - 1, list_proceed, can_relative_metric));

						// update the max element in set
						set_iter = can_opt_set.end();
						set_iter--;
					}
				}
				else {
					can_opt_set.emplace(metric_point_opt(former_state, i - 1, list_proceed, can_relative_metric));
				}

				//can_opt.push(metric_point(former_state, i - 1, can_relative_metric, los, list_proceed));

				/*can_opt.push(metric_point(former_state, i - 1, path_metric(former_state, i - 1) + dist_tmp \
					- path_metric(last_opt_state(i), i) + relative_metric(list_proceed), last_opt_state(i)));*/

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

		auto iter = can_opt_set.begin();
		// now you get the second opt metric path, backtracking it
		list_proceed++;
		//cout << "list_proceed=" << list_proceed << endl;
		diverge_time = iter->diverge_time;
		diverge_time++;		// adding 1 temporary, this become merge bit, for less computation only
		//cout << " can_opt.top() = " << can_opt.top() << endl;

		int diverge_list_ind = iter->list_belonging;
		relative_metric(list_proceed) = iter->relative_metric;		// right

		//int diverge_list_ind = finished_state(can_opt.top().merge_state, diverge_time);
		// copy the path before diverge
		for (int i = n - 1; i > diverge_time && i >= 0; --i) {
			list_v_hat(list_proceed, i) = list_v_hat(diverge_list_ind, i);
		}

		if (diverge_time >= 0) {
			// the new list is diverge from Viterbi

			// flip the merge bit
			list_v_hat(list_proceed, diverge_time) = list_v_hat(diverge_list_ind, diverge_time).flip();

			diverge_time--;		// minus 1, change the variable back
			int track_state = iter->diverge_state;
			for (int i = diverge_time; i >= 0; --i) {
				last_opt_state(i) = track_state;
				track_state = state[i][track_state].past;
				// update the current state a time instant earlier
				list_v_hat(list_proceed, i) = track_state != last_opt_state(i);
				// if state changed, the decoded bit is 1, else 0

				/*last_opt_state(i) = track_state;
				list_v_hat(list_proceed, i) = past_state(track_state, i) != track_state;
				track_state = past_state(track_state, i);	*/
			}
			//cout << "list_v_hat" << list_v_hat;

			// adding the zero position bits
			for (int i = n; i < n_total; ++i) {
				list_v_hat(list_proceed, i) = zero_col_hdr(i - n);
			}
			
			add_zero_col_pos_2_can_opt_set(to_be_found_list_num - 1);
		}
		else {
			// the new list is diverge from zero columns

			//cout << "diverge_list_ind = " << diverge_list_ind << endl;

			// recover diverge_time;
			diverge_time--;				// remember to minus this back to get the correct int_can
			int zero_col_hdr_int_can = -(diverge_time + 10);
			//cout << "zero_col_hdr_int_can = " << zero_col_hdr_int_can << endl;
			
			// adding the zero position bits
			for (int i = n_total - 1; i >= n; --i) {
				list_v_hat(list_proceed, i) = zero_col_hdr_int_can & 1;
				zero_col_hdr_int_can >>= 1;
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
	Matrix<GF2> decode_v_repeatedly(Matrix<T> r_or_hdr, int list_num = 1) {

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

		// 5 times computation complexity to viterbi for cycilic codes


		// to merge the method of hard and soft decoding, if T is ordered structure, 
		// then use soft decoidng, else use hard deocding
		bool is_soft = r_or_hdr.get_is_ordered_structure();

		//cout << "is_soft = " << is_soft << endl;
		//cout << "n = " << n << endl;

		// throw away symbols w.r.t. zero columns, n_total-> n
		r_or_hdr.permute(PM_permutation);
		if (n == n_total);
		else {
			r_or_hdr = r_or_hdr.get_part(0, 0, 0, n - 1);
			// to be fixed, or discarded
		}

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

		// for i==0, the first received symbol (to be discarded)
		//path_metric(0, 0) = dist(0, 0);
		//past_state(0, 0) = 0;					// initial and final state is 0 for a linear block code, set the past state be 0

		// initialize past state
		//path_metric(PM_col_reps(0), 0) = dist(1, 0);
		//past_state(PM_col_reps(0), 0) = 0;		// second path derived from initial 0 state

		//ope_num_after = my_double_auxiliary_storage::operation_number;
		//cout << "ope_num = " << ope_num_after - ope_num_before << endl;

		//cout << "PM_col_reps" << PM_col_reps;

		if (state[0].size() == 0) {
			initialize_state();		// initialize the structure of state
		}

		// for time i, state[i] mapping from state number to state metric and past state pair
		state[0][0] = past_and_metric(0, dist(0, 0));	// if exist, means active

		if (PM_col_reps(0) != 0) {
			// this is only valid if ending_row_ind != 0
			if (state[0].size() != 1) {
				// there are 2 states at time 0
				state[0][PM_col_reps(0)] = past_and_metric(0, dist(1, 0));
			}
		}
		else {
			// contention at state[0]
			state[0][0].metric = dist(0, 0) < dist(0, 1) ? dist(0, 0) : dist(0, 1);
		}

		//Matrix<bool> is_explored(1, pow_2_r, '0');
		//is_explored(0) = true;
		//is_explored(PM_col_reps(0)) = true;
		int no_contention_no_divergence_state_ind = 0;

		for (int i = 1; i < n; ++i) {
			int test_new_state = PM_col_reps(i);
			int former = i - 1;

			//cout << "---------- i= " << i << "----------" << endl;
			//for (auto iter = state[i].begin(); iter != state[i].end(); ++iter) {
			//	// current state, last state, and metric
			//	cout << iter->first << ": " << iter->second << endl;
			//}

			switch (time_seperation(i)) {
			case 0:
				// divergence with no contentioin
				//cout << "i=" << i << ", divergence with no contentioin" << endl;
				for (auto iter = state[former].begin(); iter != state[former].end(); ++iter) {
					int j = iter->first;	// starting state
					my_double last_metric = iter->second.metric;
					state[i][j] = past_and_metric(j, last_metric + dist(0, i));

					/*path_metric(j, i) = path_metric(j, i - 1) + dist(0, i);
					past_state(j, i) = j;*/

					//int new_state = j ^ test_new_state;
					//is_explored(new_state) = true;

					state[i][j ^ test_new_state] = past_and_metric(j, last_metric + dist(1, i));

					/*path_metric(new_state, i) = path_metric(j, i - 1) + dist(1, i);
					past_state(new_state, i) = j;*/
				}
				break;
			case 1:
				// no divergence, no contention, very similar to case 0, but we would not love to combine them
				//cout << "i=" << i << ", no contention" << endl;
				for (auto iter = state[former].begin(); iter != state[former].end(); ++iter) {
					int j = iter->first;	// starting state
					my_double last_metric = iter->second.metric;
					// we assume that state[i]'s ording is not change, i.e., no new insertion and deletion

					// the following may be embarrasing, the condition switching may take time
					switch (no_contention_no_divergence_state[no_contention_no_divergence_state_ind]) {
					case 0:
						state[i][j] = past_and_metric(j, last_metric + dist(0, i));
						break;
					case 1:
						state[i][j ^ test_new_state] = past_and_metric(j, last_metric + dist(1, i));
						//for case 2, do nothing
					}
					no_contention_no_divergence_state_ind++;
				}
				break;
			case 2:
				// contention, the new column is dependent of key of state[former], which forms a linear space in GF2
				//cout << "i=" << i << ", contention" << endl;

				min_contention_col_ind = min_contention_col_ind == -1 ? i : min_contention_col_ind;	// first contention ind
				for (auto iter = state[i].begin(); iter != state[i].end(); ++iter) {
					int j = iter->first;	// starting state
					int new_state = j ^ test_new_state;

					// if j in state[i], j must be in key of state[former], since it is not erased by ending state
					// the key of state[former] always forms a complete linear space of GF2.
					// also j ^ test_new_state in key of state[former], as we analyzed above

					// update correlation distance and past state
					my_double can_1 = state[former][j].metric + dist(0, i);
					// accumulate correlation distance for each path to state j
					my_double can_2 = state[former][new_state].metric + dist(1, i);	// key complexity

					// choose the best path for each merging state
					iter->second = can_1 < can_2 ? past_and_metric(j, can_1) : past_and_metric(new_state, can_2);

					/*bool is_can_1_better = can_1 < can_2;
					path_metric(j, i) = is_can_1_better ? can_1 : can_2;
					past_state(j, i) = is_can_1_better ? j : (j ^ test_new_state);*/
				}
			}

			//cout << "---------- after ----------" << endl;
			//for (auto iter = state[i].begin(); iter != state[i].end(); ++iter) {
			//	// current state, last state, and metric
			//	cout << iter->first << ": " << iter->second << endl;
			//}
		}

		//cout << "state" << endl;
		//for (int i = 0; i < n; ++i) {
		//	for (auto iter = state[i].begin(); iter != state[i].end(); ++iter) {
		//		// current state, last state, and metric
		//		cout << iter->first << ": " << iter->second <<endl;
		//	}
		//	cout << "==========" << endl;
		//}

		/*cout << "path_metric" << path_metric.get_part(0, 0, 3, n - 1);
		cout << "past_state" << past_state.get_part(0, 0, 3, n - 1);*/

		list_v_hat.resize(list_num, n, false);		// to be fixed to n_total
		relative_metric.resize(1, list_num, false);

		// set for listing

		// get the optimum decoded codeword, final state is 0
		int track_state = 0;					// track the current state, for final state, it must be 0
		for (int i = n - 1; i >= 0; --i) {		// path back tracking from time instance n-1 to 0
			last_opt_state(i) = track_state;
			track_state = state[i][track_state].past;
			list_v_hat(0, i) = track_state != last_opt_state(i);

			//list_v_hat(0, i) = past_state(track_state, i) != track_state;	// if state changed, the decoded bit is 1 else it is 0			
			//track_state = past_state(track_state, i);						// update the current state a time instant earlier
		}
		//cout << "list_v_hat" << list_v_hat;

		valid_list_ind.resize(1, 0);		// not releasing the space, just seting row and column
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
			//cout << "relative_metric" << relative_metric;
		}

		//while (!can_opt.empty()) {
		//	cout << can_opt.top();
		//	can_opt.pop();
		//}
		//cout << endl;
		list_v_hat = list_v_hat.get_rows(valid_list_ind);

		//cout << "list_v_hat" << list_v_hat;
#ifdef RUN_MSG
#ifdef count_operation_number
#ifdef use_my_double

		double_ope_num_after = my_double_auxiliary_storage::operation_number;
		cout << "(Viterbi_optimized) double_ope_num = " << double_ope_num_after - double_ope_num_before << endl;
#endif // use_my_double

		GF2_ope_num_after = GF2_auxiliary_storage::operation_number;
		cout << "(Viterbi_optimized) GF2_ope_num = " << GF2_ope_num_after - GF2_ope_num_before << endl;
#endif // count_operation_number
#endif // RUN_MSG



		return list_v_hat.get_part(0, 0, list_num - 1, n - 1);
	}

	void print_rank_distribution() {
		cout << "----rank_distribution----" << endl;
		int s = rank_distribution.size();
		for (int i = 0; i < s; ++i) {
			if (rank_distribution(i) == 0);
			else {
				cout << "rank = " << i + 1 << ", num = " << rank_distribution(i) << endl;
			}
		}
	}

	void PM_preprocessing() {
		PM_permutation = find_opt_PM::column_permute(PM);
	}
};

/**
 * space saving version, using self-developed sorted vector.
 * but slower because of fetching elements induce O(log n) iteration 
 * seems that it is only valid for rate > 1/2, be careful to use it, not fixing, please discard this class
 */
class Viterbi_Sorted_vector {
private:
	int n;
	int r;
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
	int min_contention_col_ind;
	int max_state_num;

	// new variable for list viterbi
	priority_queue<ending_one> fixed_PM_row_end_one_ind;
	vector<Sorted_vector<past_and_metric>> state_linear;
	// use linear container to store past and metric, will use less memory, suitable for repeatedly using same parity matrix as decoder

	Matrix<char> time_seperation;			// the i^th element indicates the state changing at time i, 
	// 0 means divergence with no contention, 1 means no divergence and no contention, 2 means no divergence and contention
	Matrix<char> no_contention_no_divergence_state;		// for no divergence and no contention time instance, 
	// indicates how to update each state, 0 means that update by path 0, 1 menas that updates by path 1, and 2 means that no updates

	Matrix<GF2> unused_PM;
	Matrix<int> valid_list_ind;

	void initialize_size() {
		/* apply space */

		PM_col_reps = Matrix<int>(1, n);
		dist = Matrix<my_double>(2, n);

		/* list decoding */
		last_opt_state = Matrix<int>(1, n);

	}
	void initialize_PM() {
		for (int i = 0; i < n; ++i) {
			PM_col_reps(i) = 0;
			for (int p = 0; p < r; ++p) {
				PM_col_reps(i) <<= 1;
				PM_col_reps(i) += (int)PM(p, i);
			}
		}
		// in may be erasing states when rank right-part of the PM decreases, 
		// for simplicity we only consider zeros ending that states
		fixed_PM_row_end_one_ind = priority_queue<ending_one>();		// claer out
		for (int i = 0; i < r; ++i) {
			int j = n - 1;
			while (j >= 0 && PM(i, j) == 0) {
				--j;
			}
			fixed_PM_row_end_one_ind.push(ending_one(i, j));
		}
		min_contention_col_ind = -1;
	}

public:
	int log2_max_state_num;
	Matrix<GF2> test_zero;

	Viterbi_Sorted_vector(const Matrix<GF2>& _PM) {
		// get the code parameters, for a linear block code of (n,k) and r = n - k, 
		// with parity check matrix PM

		n = _PM.col();
		r = _PM.row();
		initialize_size();

		PM = _PM;
		initialize_PM();

		// this is for repeatly use the same PM, initialize state at constructor
		//initialize_state();
	}
	Viterbi_Sorted_vector(int _n, int _r) {
		n = _n;
		r = _r;
		initialize_size();
	}

	void initialize_state() {
		// initialize the sturcture of state, this is for using the same parity check matrix repeatedly
		priority_queue<ending_one> PM_row_end_one_ind = fixed_PM_row_end_one_ind;
		vector<unordered_map<int, past_and_metric>> state(n);

		int next_ending_one = PM_row_end_one_ind.top().one_pos_ind;
		int ending_row_ind = PM_row_end_one_ind.top().row_ind;
		time_seperation = Matrix<char>(1, n);
		time_seperation.reset(-1);
		no_contention_no_divergence_state = Matrix<char>(1, 0, 'v');

		// for time i, state[i] mapping from state number to state metric and past state pair
		state[0][0] = past_and_metric();	// if exist, means active

		max_state_num = 0;

		// this is only valid if ending_row_ind != 0
		if (next_ending_one != 0) {
			state[0][PM_col_reps(0)] = past_and_metric();
		}
		else {
			while (next_ending_one == 0) {	// clear all rows ending with last non-zero on position 0
				PM_row_end_one_ind.pop();
				if (!PM_row_end_one_ind.empty()) {
					next_ending_one = PM_row_end_one_ind.top().one_pos_ind;
					ending_row_ind = PM_row_end_one_ind.top().row_ind;
				}
				else {
					// this is the ending of decoding
					break;
				}
			}
		}

		max_state_num = max_state_num < (int)state[0].size() ? (int)state[0].size() : max_state_num;

		//Matrix<bool> is_explored(1, pow_2_r, '0');
		//is_explored(0) = true;
		//is_explored(PM_col_reps(0)) = true;

		for (int i = 1; i < n; ++i) {
			int test_new_state = PM_col_reps(i);
			int former = i - 1;

			// copy the former states
			state[i] = state[former];

			if (state[former].find(test_new_state) == state[former].end()) {
				// new state should be opened
				for (auto iter = state[former].begin(); iter != state[former].end(); ++iter) {
					// add the state
					int j = (iter->first) ^ test_new_state;

					// this make sure that state (iter->first) must in key of state[former]
					// if j in key of state[i], then j ^ test_new_state in key of state[former]

					//if (is_explored(j) == false)		// must be false
					state[i][j] = past_and_metric();
				}
			}

			while (i == next_ending_one) {
				// removing the ending state at receiver
				int big_mark = 1 << (r - ending_row_ind - 1);
				auto iter = state[i].begin();
				while (iter != state[i].end()) {
					if (((iter->first) & big_mark) == big_mark) {
						state[i].erase(iter++);
					}
					else {
						iter++;
					}
				}

				PM_row_end_one_ind.pop();
				if (!PM_row_end_one_ind.empty()) {
					next_ending_one = PM_row_end_one_ind.top().one_pos_ind;
					ending_row_ind = PM_row_end_one_ind.top().row_ind;
				}
				else {
					// this is the ending of decoding
					break;
				}
			}

			//cout << "---------- i= " << i << "----------" << endl;
			//for (auto iter = state[i].begin(); iter != state[i].end(); ++iter) {
			//	// current state, last state, and metric
			//	cout << iter->first << ": " << iter->second << endl;
			//}

			if (state[former].find(test_new_state) == state[former].end()) {
				if (state[i].find(test_new_state) != state[i].end()) {
					// divergence with no contentioin
					time_seperation(i) = 0;
					////cout << "i=" << i << ", divergence with no contentioin" << endl;

					//for (auto iter = state[former].begin(); iter != state[former].end(); ++iter) {
					//	int j = iter->first;	// starting state
					//	state[i][j] = past_and_metric();

					//	/*path_metric(j, i) = path_metric(j, i - 1) + dist(0, i);
					//	past_state(j, i) = j;*/

					//	int new_state = j ^ test_new_state;
					//	//is_explored(new_state) = true;

					//	state[i][new_state] = past_and_metric();

					//	/*path_metric(new_state, i) = path_metric(j, i - 1) + dist(1, i);
					//	past_state(new_state, i) = j;*/
					//}
				}
				else {
					// no divergence, no contention
					time_seperation(i) = 1;
					//cout << "i=" << i << ", no contention" << endl;

					for (auto iter = state[former].begin(); iter != state[former].end(); ++iter) {
						// the following may be embarrasing, the condition of if and program switching may take time

						int j = iter->first;	// starting state

						auto update_state = state[i].find(j);
						if (update_state != state[i].end()) {		// path 0, make sure not creating new state
							//update_state->second = past_and_metric();
							no_contention_no_divergence_state.push_back(0);

							//path_metric(j, i) = path_metric(j, i - 1) + dist(0, i);
							//past_state(j, i) = j;
						}
						else {		// at most one path will exist
							int new_state = j ^ test_new_state;
							update_state = state[i].find(new_state);
							if (update_state != state[i].end()) {		// path 1, whic may be not exists
								//update_state->second = past_and_metric();
								no_contention_no_divergence_state.push_back(1);

								/*path_metric(new_state, i) = path_metric(j, i - 1) + dist(1, i);
								past_state(new_state, i) = j;*/
							}
							else {
								no_contention_no_divergence_state.push_back(2);
							}
						}

					}
				}
			}
			else {
				// contention, the new column is dependent of key of state[former], which forms a linear space in GF2
				time_seperation(i) = 2;
				//cout << "i=" << i << ", contention" << endl;

				// to be fixed to -1, or discarded
				min_contention_col_ind = min_contention_col_ind == 0 ? i : min_contention_col_ind;	// first contention ind
				//for (auto iter = state[i].begin(); iter != state[i].end(); ++iter) {
				//	int j = iter->first;	// ending state
				//	int new_state = j ^ test_new_state;

					// if j in state[i], j must be in key of state[former], since it is not erased by ending state
					// the key of state[former] always forms a complete linear space of GF2.
					// also j ^ test_new_state in key of state[former], as we analyzed above

					// update correlation distance and past state
					// accumulate correlation distance for each path to state j

					// choose the best path for each merging state
				//	state[i][j] = past_and_metric();

					/*bool is_can_1_better = can_1 < can_2;
					path_metric(j, i) = is_can_1_better ? can_1 : can_2;
					past_state(j, i) = is_can_1_better ? j : (j ^ test_new_state);*/
					//}
			}

			//cout << "---------- after ----------" << endl;
			//for (auto iter = state[i].begin(); iter != state[i].end(); ++iter) {
			//	// current state, last state, and metric
			//	cout << iter->first << ": " << iter->second << endl;
			//}

			max_state_num = max_state_num < (int)state[i].size() ? (int)state[i].size() : max_state_num;
		}

		state_linear = vector<Sorted_vector<past_and_metric>>(n);
		Matrix<int> valid_index(1, max_state_num);
		for (int i = 0; i < n; ++i) {
			valid_index.resize(1, 0, false);		// starting
			for (auto iter = state[i].begin(); iter != state[i].end(); ++iter) {
				valid_index.push_back(iter->first);
			}
			//cout << "valid_index" << valid_index;
			valid_index.sort();
			//cout << "valid_index" << valid_index;
			int vs = valid_index.size();
			state_linear[i] = Sorted_vector<past_and_metric>(vs);
			for (int j = 0; j < vs; ++j) {
				state_linear[i].assign_ind(j, valid_index(j));
			}

			//cout << "state_linear[" << i << "]" << state_linear[i];
		}
	}
	/**
	 * .change PM not changing the size of code, i.e., n, and k not changed.
	 * not initializing state, for decode_v_once especially.
	 */
	void change_PM(const Matrix<GF2>& _PM) {
		PM = _PM;
		initialize_PM();
	}
	/**
	 * .change PM not changing the size of code, i.e., n, and k not changed.
	 * not initializing state, for decode_v_once especially.
	 */
	void change_unused_PM(const Matrix<GF2>& _unused_PM) {
		unused_PM = _unused_PM;
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
	Matrix<GF2> decode_v_repeatedly_space_saving(const Matrix<T>& r_or_hdr, int list_num = 1) {
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

		// 10 times computation complexity to viterbi for cycilic codes, but space consuming same as viterbi for cycilic codes

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

		// for i==0, the first received symbol (to be discarded)
		//path_metric(0, 0) = dist(0, 0);
		//past_state(0, 0) = 0;					// initial and final state is 0 for a linear block code, set the past state be 0

		// initialize past state
		//path_metric(PM_col_reps(0), 0) = dist(1, 0);
		//past_state(PM_col_reps(0), 0) = 0;		// second path derived from initial 0 state

		//ope_num_after = my_double_auxiliary_storage::operation_number;
		//cout << "ope_num = " << ope_num_after - ope_num_before << endl;

		//cout << "PM_col_reps" << PM_col_reps;

		//if (state_linear[0].size() == 0) {
		//	initialize_state();		// initialize the structure of state
		//}

		// for time i, state[i] mapping from state number to state metric and past state pair
		state_linear[0][0] = past_and_metric(0, dist(0, 0));	// if exist, means active

		if (PM_col_reps(0) != 0) {
			// this is only valid if ending_row_ind != 0
			if (state_linear[0].size() != 1) {
				// there are 2 states at time 0
				state_linear[0][PM_col_reps(0)] = past_and_metric(0, dist(1, 0));
			}
		}
		else {
			// contention at state[0]
			state_linear[0][0].metric = dist(0, 0) < dist(0, 1) ? dist(0, 0) : dist(0, 1);
		}

		//Matrix<bool> is_explored(1, pow_2_r, '0');
		//is_explored(0) = true;
		//is_explored(PM_col_reps(0)) = true;
		int no_contention_no_divergence_state_ind = 0;

		for (int i = 1; i < n; ++i) {
			int test_new_state = PM_col_reps(i);
			int former = i - 1;

			//cout << "---------- i= " << i << "----------" << endl;
			//for (auto iter = state_linear[i].begin(); iter != state_linear[i].end(); ++iter) {
			//	// current state_linear, last state_linear, and metric
			//	cout << iter->first << ": " << iter->second << endl;
			//}

			int sss;
			switch (time_seperation(i)) {
			case 0:
				// divergence with no contentioin
				//cout << "i=" << i << ", divergence with no contentioin" << endl;
				sss = state_linear[former].size();
				for (int q = 0; q < sss; ++q) {
					int j = state_linear[former].key(q);	// starting state_linear
					my_double last_metric = state_linear[former].val(q).metric;
					state_linear[i][j] = past_and_metric(j, last_metric + dist(0, i));

					/*path_metric(j, i) = path_metric(j, i - 1) + dist(0, i);
					past_state(j, i) = j;*/

					//int new_state = j ^ test_new_state;
					//is_explored(new_state) = true;

					state_linear[i][j ^ test_new_state] = past_and_metric(j, last_metric + dist(1, i));

					/*path_metric(new_state, i) = path_metric(j, i - 1) + dist(1, i);
					past_state(new_state, i) = j;*/
				}
				break;
			case 1:
				// no divergence, no contention, very similar to case 0, but we would not love to combine them
				//cout << "i=" << i << ", no contention" << endl;

				sss = state_linear[former].size();
				for (int q = 0; q < sss; ++q) {
					// we assume that state_linear[i]'s ording is not change, i.e., no new insertion and deletion

					// the following may be embarrasing, the condition switching may take time

					int j = state_linear[former].key(q);	// starting state_linear
					my_double last_metric = state_linear[former].val(q).metric;

					switch (no_contention_no_divergence_state[no_contention_no_divergence_state_ind]) {
					case 0:
						state_linear[i][j] = past_and_metric(j, last_metric + dist(0, i));
						break;
					case 1:
						state_linear[i][j ^ test_new_state] = past_and_metric(j, last_metric + dist(1, i));
						//for case 2, do nothing
					}
					no_contention_no_divergence_state_ind++;
				}
				break;
			case 2:
				// contention, the new column is dependent of key of state_linear[former], which forms a linear space in GF2
				//cout << "i=" << i << ", contention" << endl;

				min_contention_col_ind = min_contention_col_ind == -1 ? i : min_contention_col_ind;		// first contention ind
				sss = state_linear[i].size();
				for (int q = 0; q < sss; ++q) {
					int j = state_linear[i].key(q);	// ending state_linear
					int new_state = j ^ test_new_state;

					// if j in state_linear[i], j must be in key of state_linear[former], since it is not erased by ending state_linear
					// the key of state_linear[former] always forms a complete linear space of GF2.
					// also j ^ test_new_state in key of state_linear[former], as we analyzed above

					// update correlation distance and past state_linear
					my_double can_1 = state_linear[former][j].metric + dist(0, i);
					// accumulate correlation distance for each path to state_linear j
					my_double can_2 = state_linear[former][new_state].metric + dist(1, i);	// key complexity

					// choose the best path for each merging state_linear
					state_linear[i].val(q) = can_1 < can_2 ? past_and_metric(j, can_1) : past_and_metric(new_state, can_2);

					/*bool is_can_1_better = can_1 < can_2;
					path_metric(j, i) = is_can_1_better ? can_1 : can_2;
					past_state(j, i) = is_can_1_better ? j : (j ^ test_new_state);*/
				}
			}

			//cout << "---------- after ----------" << endl;
			//for (auto iter = state_linear[i].begin(); iter != state_linear[i].end(); ++iter) {
			//	// current state_linear, last state_linear, and metric
			//	cout << iter->first << ": " << iter->second << endl;
			//}
		}

		//cout << "state_linear" << endl;
		//for (int i = 0; i < n; ++i) {
		//	for (auto iter = state_linear[i].begin(); iter != state_linear[i].end(); ++iter) {
		//		// current state_linear, last state_linear, and metric
		//		cout << iter->first << ": " << iter->second <<endl;
		//	}
		//	cout << "==========" << endl;
		//}

		/*cout << "path_metric" << path_metric.get_part(0, 0, 3, n - 1);
		cout << "past_state" << past_state.get_part(0, 0, 3, n - 1);*/

		list_v_hat.resize(list_num, n, false);
		relative_metric.resize(1, list_num, false);

		// set for listing

		// get the optimum decoded codeword, final state_linear is 0
		int track_state = 0;					// track the current state_linear, for final state_linear, it must be 0
		for (int i = n - 1; i >= 0; --i) {		// path back tracking from time instance n-1 to 0
			last_opt_state(i) = track_state;
			track_state = state_linear[i][track_state].past;
			list_v_hat(0, i) = track_state != last_opt_state(i);
			
			//// if state_linear changed, the decoded bit is 1 else it is 0
			//list_v_hat(0, i) = past_state(track_state, i) != track_state;	
			//track_state = past_state(track_state, i);						// update the current state_linear a time instant earlier
		}
		//cout << "list_v_hat" << list_v_hat;

		valid_list_ind.resize(1, 0);		// not releasing the space, just seting row and column
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
				next_subopt_v_repeatedly_space_saving();

				/* since the parity matrix is determined, we if unused_PM exists, it must be necessary */
				if (!is_unused_PM_exist || (is_unused_PM_exist && unused_PM.check_inner_product(list_v_hat.get_row(i)))) {
					valid_list_ind.push_back(i);
				}
			}
			//cout << "relative_metric" << relative_metric;
		}

		//while (!can_opt.empty()) {
		//	cout << can_opt.top();
		//	can_opt.pop();
		//}
		//cout << endl;

		//cout << "list_v_hat" << list_v_hat;

#ifdef RUN_MSG
#ifdef count_operation_number
#ifdef use_my_double

		double_ope_num_after = my_double_auxiliary_storage::operation_number;
		cout << "(Viterbi_Sorted_vector) double_ope_num = " << double_ope_num_after - double_ope_num_before << endl;
#endif // use_my_double

		GF2_ope_num_after = GF2_auxiliary_storage::operation_number;
		cout << "(Viterbi_Sorted_vector) GF2_ope_num = " << GF2_ope_num_after - GF2_ope_num_before << endl;
#endif // count_operation_number
#endif // RUN_MSG

		return list_v_hat.get_part(0, 0, list_num - 1, n - 1);
	}

	/**
	 * . find the next sub-optimum codeword, with parameter set in private area
	 */
	void next_subopt_v_repeatedly_space_saving() {
		// we have problem here

		//cout << "r=" << r << endl;
		//cout << "min_contention_col_ind=" << min_contention_col_ind << endl;

		// consider the space saving later, i.e., we donot need to store the first r state_linear, since no contention
		for (int i = diverge_time; i >= min_contention_col_ind; --i) {
			int los = last_opt_state(i);
			//if (pointer_pm->list_belonging != -1)
			//	continue;		// merge into proceed state_linear, do not need to continue

			/*if (finished_state(last_opt_state(i), i) != -1)
				break;*/

			int former = i - 1;
			bool is_best_path_0 = list_v_hat(list_proceed, i) == 0;
			// take the worse path, but we should judge whether the worse path exists
			int former_state = is_best_path_0 ? (los ^ PM_col_reps(i)) : los;
			my_double dist_tmp = is_best_path_0 ? dist(1, i) : dist(0, i);

			if (former != -1) {
				int temp_ind = state_linear[former].find(former_state);
				if (temp_ind != -1) {

					can_opt.push(metric_point(former_state, former, state_linear[former].val(temp_ind).metric + dist_tmp \
						- state_linear[i][los].metric + relative_metric(list_proceed), los, list_proceed));

					/*can_opt.push(metric_point(former_state, i - 1, path_metric(former_state, i - 1) + dist_tmp \
						- path_metric(last_opt_state(i), i) + relative_metric(list_proceed), last_opt_state(i)));*/

				}
			}
			else {
				can_opt.push(metric_point(former_state, former, dist_tmp \
					- state_linear[i][los].metric + relative_metric(list_proceed), los, list_proceed));
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
			track_state = state_linear[i][track_state].past;				// update the current state_linear a time instant earlier
			list_v_hat(list_proceed, i) = track_state != last_opt_state(i);	// if state_linear changed, the decoded bit is 1, else 0

			/*last_opt_state(i) = track_state;
			list_v_hat(list_proceed, i) = past_state(track_state, i) != track_state;
			track_state = past_state(track_state, i);	*/
		}
		//cout << "list_v_hat" << list_v_hat;

		can_opt.pop();
		//cout << "second_opt.empty() = " << can_opt.empty() << endl;
	}

};
