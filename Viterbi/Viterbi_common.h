/*****************************************************************//**
 * \file   Viterbi_common.h
 * \brief  common data structure and head file for Viterbi algorithm
 * 
 * \author 26259
 * \date   April 2023
 *********************************************************************/

#pragma once

#include"../GF/GF2.h"
#include"../my_lib/Matrix.h"
#include"../channel/channel.h"
#include<queue>							// will use priority queue for list viterbi decoding
using namespace std;

#ifdef use_my_double
#include"../my_lib/my_double.h"
#include"../my_lib/my_float.h"
#else
#define my_double double
#define my_float float
#endif

static const int invalid_large_metric = 6000000;		// 6 million is enough

struct metric_point {

	int diverge_state;		// the state index
	int diverge_time;		// the merge time instance
	my_double relative_metric;
	int merge_state;		// state at time diverge_time + 1, merge into exist optimum states
	int list_belonging;

	metric_point(int _diverge_state, int _diverge_time, my_double _relative_metric, int _merge_state, int _list_belonging = -1) \
		:diverge_state(_diverge_state), diverge_time(_diverge_time), relative_metric(_relative_metric), merge_state(_merge_state), \
		list_belonging(_list_belonging) {}

	bool operator < (const metric_point& mp) const {
		return relative_metric > mp.relative_metric;		// this make priority queue small value on top
	}
	friend ostream& operator << (ostream& out, const metric_point& mp) {
		out << fixed << setprecision(2);
		out << "diverge_time=" << setw(4) << mp.diverge_time << "    diverge_state=" << setw(4) << mp.diverge_state \
			<< "    merge_state=" << setw(4) << mp.merge_state << "    relative_metric=" << setw(4) << mp.relative_metric \
			<< "    list_belonging=" << setw(4) << mp.list_belonging << endl;

		out.unsetf(ios_base::fixed);
		out << setprecision(6);
		return out;
	}
};

struct ending_one {
	int row_ind;
	int one_pos_ind;
	ending_one(int _row_ind, int _one_pos_ind) :row_ind(_row_ind), one_pos_ind(_one_pos_ind) {}
	bool operator < (const ending_one& eo) const {
		return one_pos_ind > eo.one_pos_ind;		// this make priority queue small value on top
	}
	friend ostream& operator << (ostream& out, const ending_one& eo) {
		out << "(" << eo.row_ind << ", " << eo.one_pos_ind << ")|";
		return out;
	}
};

struct past_and_metric {
	int past;
	my_double metric;
	past_and_metric(int _past = -invalid_large_metric, my_double _metric = invalid_large_metric) \
		: past(_past), metric(_metric){}
	bool operator < (const past_and_metric& pm) const {
		return true;		// not an ordered structure
	}
	friend ostream& operator << (ostream& out, const past_and_metric& eo) {
		out << "(" << eo.past << ", " << eo.metric << ") ";
		return out;
	}
};

/**
 * .metric point designed for Viterbi opt, maybe we can combine it with class metric_point
 */
struct metric_point_opt {

	int diverge_state;		// the state index
	int diverge_time;		// the merge time instance
	int list_belonging;
	my_double relative_metric;

	metric_point_opt() = default;

	metric_point_opt(int _diverge_state, int _diverge_time, int _list_belonging, my_double _relative_metric) \
		:diverge_state(_diverge_state), diverge_time(_diverge_time), list_belonging(_list_belonging), relative_metric(_relative_metric) {}

	bool operator < (const metric_point_opt& mp) const {
		return relative_metric < mp.relative_metric;
	}
	bool operator <= (const metric_point_opt& mp) const {
		return relative_metric <= mp.relative_metric;
	}
	bool operator > (const metric_point_opt& mp) const {
		return relative_metric > mp.relative_metric;
	}
	bool operator >= (const metric_point_opt& mp) const {
		return relative_metric >= mp.relative_metric;
	}

	friend ostream& operator << (ostream& out, const metric_point_opt& mp) {
		out << fixed << setprecision(2);
		out << "diverge_time=" << setw(4) << mp.diverge_time << "    diverge_state=" << setw(4) << mp.diverge_state \
			<< "    relative_metric=" << setw(4) << mp.relative_metric << "    list_belonging=" << setw(4) << mp.list_belonging << endl;

		out.unsetf(ios_base::fixed);
		out << setprecision(6);
		return out;
	}
};

class find_opt_PM {
public:
	static Matrix<int> further_column_permutation_iteration;

	static void solve(const Matrix<GF2>& PM, int iteration, const Matrix<GF2>& GM = Matrix<GF2>()) {

		int r = PM.row();
		int n = PM.col();
		Matrix<GF2> PM_store(r, n);
		Matrix<int> natual(1, n, 'N');
		Matrix<int> rand_permutation(1, n);
		Matrix<int> optimze_permutation(1, n);
		Matrix<int> total_permutation(1, n);
		Matrix<int> state_num(1, n);
		long long min_ope_num = LLONG_MAX;

		int max_column_iteration = 20;
		int column_iteration;
		further_column_permutation_iteration.resize(1, max_column_iteration + 1, false);
		further_column_permutation_iteration.reset(0);

		// the original state number and operation number 
		int PM_form;
		PM_store = PM;
		optimze_permutation = optimize_PM(PM_store, PM_form, column_iteration, max_column_iteration);
		state_num = counting_states(PM_store);
		long long can_ope_num = ope_num_estimation(state_num);
		cout << endl << "----------- i = " << 0 << " (original PM) -------------" << endl;
		min_ope_num = can_ope_num;
		cout << "PM_store" << PM_store;
		cout << "rand_permutation" << natual;
		total_permutation = natual;
		total_permutation.permute(optimze_permutation);
		cout << "total_permutation" << total_permutation;
		cout << "state_num" << state_num;
		cout << "can_ope_num = " << can_ope_num << endl;

		// validation
		if (GM.size() != 0) {
			cout << "(GM * PM.Transpose()).isZero() = " << (GM * PM.Transpose()).isZero() << endl;
		}


		for (int i = 1; i < iteration; ++i) {
			rand_permutation = natual.get_random_element(n);
			PM_store = PM;
			PM_store.permute_col(rand_permutation);
			int PM_form;
			optimze_permutation = optimize_PM(PM_store, PM_form, column_iteration, max_column_iteration);
			state_num = counting_states(PM_store);
			long long can_ope_num = ope_num_estimation(state_num);
			if (can_ope_num > min_ope_num);
			else {
				cout << endl << "----------- i = " << i << " -------------" << endl;
				min_ope_num = can_ope_num;
				cout << "PM_form = " << PM_form << endl;
				cout << "column_iteration = " << column_iteration << endl;
				cout << "PM_store" << PM_store;
				cout << "rand_permutation" << rand_permutation;
				//cout << "optimze_permutation" << optimze_permutation;
				total_permutation = rand_permutation;
				total_permutation.permute(optimze_permutation);
				cout << "total_permutation" << total_permutation;
				cout << "state_num" << state_num;
				cout << "can_ope_num = " << can_ope_num << endl;

				// validation
				if (GM.size() != 0) {
					PM_store.permute_col_back(total_permutation);
					cout << "(GM * PM_store.Transpose()).isZero() = " << (GM * PM_store.Transpose()).isZero() << endl;
				}
			}

			if (i % 2000000 == 0) {
				cout << "//----------- i = " << i << " -------------//" << endl;
			}
		}

		further_column_permutation_iteration.cut_end_0();
	}

	/**
	 * .given a parity check Matrix, transform the Matrix by column permutation and row transformation
	 * to give a desired form of as much '1' at the border of the Matrix, some '0' at the border is allowed
	 *
	 * \return the permutation that lead to the transformed PM
	 */
	static Matrix<int> optimize_PM(Matrix<GF2>& PM, int& PM_form, int& column_iteration, int max_column_iteration = 100) {
		// PM should be randomly permuted first

		int r = PM.row();
		int n = PM.col();
		Matrix<int> final_permutation(1, n, 'N');
		Matrix<int> permutation_record;

		//cout << "PM" << PM;
		PM.row_transformation_to_low_triangle();
		permutation_record = PM.col_permute_to_full_rank_on_right();
		final_permutation.permute(permutation_record);
		//cout << "PM: full rank on right" << PM;		// if the right part of the Matrix is rank defficient, error occur
		PM.row_transformation_right_low_triangle_to_identity();
		//cout << "PM: low triangle" << PM;

		PM.row_transformation_to_up_triangle();

		permutation_record = column_permute(PM);
		final_permutation.permute(permutation_record);

		PM_form = is_final_form(PM);
		int iterations = 0;
		for (; iterations < max_column_iteration && PM_form != 0; ++iterations) {
			if (PM_form == 1) {
				PM.row_transformation_to_up_triangle();
				permutation_record = column_permute(PM);
				final_permutation.permute(permutation_record);
			}
			else if (PM_form == 2) {
				PM.row_transformation_to_low_triangle();
				permutation_record = column_permute(PM);
				final_permutation.permute(permutation_record);
			}
			else {
				// under this case PM_form == 3, we are not sure to solve it

				PM.row_transformation_to_up_triangle();
				permutation_record = column_permute(PM);
				final_permutation.permute(permutation_record);

				PM_form = is_final_form(PM);
				if (PM_form == 0) {
					break;
				}

				PM.row_transformation_to_low_triangle();
				permutation_record = column_permute(PM);
				final_permutation.permute(permutation_record);
			}
			PM_form = is_final_form(PM);
		}
		further_column_permutation_iteration(iterations)++;
		column_iteration = iterations;

		if (iterations == max_column_iteration) {
			// giving up
			//cout << "facing unsolvable Matrix" << endl;
			//cout << "PM" << PM;
			//cout << "PM_form = " << PM_form << endl;
		}

		return final_permutation;
	}

	/**
	 * .column permute of PM, this should after 'row-col-row-row' transformation
	 * to give a desired form of as much '1' at the border of the Matrix, some '0' at the border is allowed
	 *
	 * \retrun the permutation record did to PM
	 */
	static Matrix<int> column_permute(Matrix<GF2>& PM) {

		int r = PM.row();
		int n = PM.col();
		int n_prime = n;

		// continuing permute the columns of Parity check Matrix, it seems meaningless, find a book to read
		Matrix<int> permute_ind(1, n);
		// permute the middle n-2r column
		for (int j = 0; j < n; ++j) {
			// get the starting 1 and ending 1 of the column
			int i = starting_one_ind(PM, j);
			int k = ending_one_ind(PM, j);

			if (i > k) {		// case of zero column
				permute_ind(j) = n * r;
				n_prime--;
			}
			else if (i < r - 1 - k) {
				// permute the column to the front, deciding parameter k and i
				permute_ind(j) = k * r + i;
			}
			else {
				// permute the column to the back, deciding parameter i and k
				permute_ind(j) = (n - r + i) * r + k;
			}
		}

		//cout << "middle_part_permute_ind" << permute_ind;

		Matrix<int> sort_permute = permute_ind.sort_with_ind('<');
		PM.permute_col(sort_permute);
		if (n_prime == n);
		else {
			// shrink PM to erase zero columns
			PM = PM.get_part(0, 0, -1, n_prime - 1);
		}
		return sort_permute;
	}

	/**
	 * . this should be called after checking PM is final form
	 */
	static Matrix<int> counting_states(const Matrix<GF2>& PM) {
		int r = PM.row();
		int n = PM.col();

		// count saving states, note that 1 may not continue at border
		Matrix<int> states_num(1, n);

		int last_ending_1;
		int ending_1 = 0;
		int time_instance = 0;
		for (; time_instance < n && ending_1 != r - 1; ++time_instance) {
			last_ending_1 = ending_1;
			ending_1 = r - 1;		// index of the ending 1 of a column
			for (; ending_1 >= 0; --ending_1) {
				if (PM(ending_1, time_instance) == 0);
				else {
					break;
				}
			}

			// please note that ending 1 never decrease !!
			ending_1 = ending_1 < last_ending_1 ? last_ending_1 : ending_1;

			states_num(time_instance) = ending_1 + 1;
		}


		for (; time_instance < n; ++time_instance) {
			states_num(time_instance) = r;
		}


		int last_starting_1;
		int starting_1 = r - 1;
		states_num(n - 1) = 0;
		time_instance--;
		for (; time_instance >= 0 && starting_1 != 0; --time_instance) {
			last_starting_1 = starting_1;
			starting_1 = 0;			// index of the starting 1 of a column
			for (; starting_1 < r; ++starting_1) {
				if (PM(starting_1, time_instance) == 0);
				else {
					break;
				}
			}
			// please note that starting 1 never increase !!
			starting_1 = starting_1 > last_starting_1 ? last_starting_1 : starting_1;

			states_num(time_instance - 1) -= starting_1;
		}
		//cout << "states_num" << states_num;

		return states_num;
	}

	/**
	 * .a good approaximation to operation number for viterbi algorithm with list 1,
	 * we only return the adding computation for states num, but compare_computation is also estimated, not returned
	 */
	static long long ope_num_estimation(Matrix<int> states_num) {
		int n = states_num.col();

		// compute total computation
		long long total_compare_computation = 0;
		long long total_add_computation = 0;
		total_add_computation = 2;
		for (int j = 1; j < n; ++j) {
			if (states_num(j) == states_num(j - 1) + 1) {
				total_add_computation += (long long)1 << (states_num(j - 1) + 1);		// one trellis for each state
			}
			else if (states_num(j) == states_num(j - 1) || states_num(j) == states_num(j - 1) - 1) {
				total_add_computation += (long long)1 << (states_num(j) + 1);			// 2 trellis for each state
				total_compare_computation += (long long)1 << states_num(j);
			}
			else if (states_num(j) > states_num(j - 1)) {
				total_add_computation += (long long)1 << (states_num(j - 1) + 1);
				//cout << "unexpected state, over expending. j = " << j << endl;
			}
			else {
				total_add_computation += (long long)1 << (states_num(j) + 1);			// 2 trellis for each state
				total_compare_computation += (long long)1 << states_num(j);
				//cout << "unexpected state, over shrinking. j = " << j << endl;
			}
			// for optimized network, there should be no expending and shrinking
		}
		return total_add_computation;

	}

	/**
	 * .condition of final form for PM (size: r*n):
	 *
	 * (condition 1)
	 * column 0 has only 1 at row index 0
	 * for each column in index [1, n-1], the ending 1 is either equals to the former column's, or equals to 1 + the former column's
	 * where the ending 1 never decreases as the column index increases
	 *
	 * (condition 2)
	 * column n-1 has only 1 at row nidex n-1
	 * for each column in index [0, n-2], the starting 1 is either equals to the later column's, or equals to 1 + the later column's
	 * where the starting 1 never increases as the column index decreases
	 *
	 * return 0 -> PM is final form,
	 *		1 -> PM violates condition 1,
	 *		2 -> PM violates condition 2,
	 *		3 -> PM violates both condition 1 and condition 2
	 */
	static int is_final_form(const Matrix<GF2>& PM) {
		int result = 0;
		int r = PM.row();
		int n = PM.col();


		// for the first column
		int ending_1 = r - 1;		// index of the ending 1 of a column
		for (; ending_1 >= 0; --ending_1) {
			if (PM(ending_1, 0) == 0);
			else {
				break;
			}
		}
		if (ending_1 != 0) {
			result = result | 1;
		}
		else {
			for (int time_instance = 1; time_instance < n && ending_1 != r - 1; ++time_instance) {
				int last_ending_1 = ending_1;
				ending_1 = r - 1;		// index of the ending 1 of a column
				for (; ending_1 >= 0; --ending_1) {
					if (PM(ending_1, time_instance) == 0);
					else {
						break;
					}
				}

				if (ending_1 > last_ending_1 + 1) {
					result = result | 1;
					break;
				}

				// please note that ending 1 never decrease !!
				ending_1 = ending_1 < last_ending_1 ? last_ending_1 : ending_1;

			}
		}



		// for the last column
		int starting_1 = 0;
		for (; starting_1 < r; ++starting_1) {
			if (PM(starting_1, n - 1) == 0);
			else {
				break;
			}
		}
		if (starting_1 != r - 1) {
			result = result | 2;
			return result;
		}
		else {
			for (int time_instance = n - 2; time_instance >= 0 && starting_1 != 0; --time_instance) {
				int last_starting_1 = starting_1;
				starting_1 = 0;			// index of the starting 1 of a column
				for (; starting_1 < r; ++starting_1) {
					if (PM(starting_1, time_instance) == 0);
					else {
						break;
					}
				}
				if (starting_1 < last_starting_1 - 1) {
					result = result | 2;
					break;
				}

				// please note that starting 1 never increase !!
				starting_1 = starting_1 > last_starting_1 ? last_starting_1 : starting_1;
			}
		}

		return result;
	}

	// find the starting 1 of column j
	static int starting_one_ind(const Matrix<GF2>& PM, int col_ind) {
		int r = PM.row();
		int i = 0;			// index of the starting 1
		for (; i < r; ++i) {
			if (PM(i, col_ind) == 0);
			else {
				break;		// it must break at some point
			}
		}
		return i;
	}

	// find the ending 1 of column j
	static int ending_one_ind(const Matrix<GF2>& PM, int col_ind) {
		int r = PM.row();
		int k = r - 1;		// index of the ending 1
		for (; k >= 0; --k) {
			if (PM(k, col_ind) == 0);
			else {
				break;		// it must break at some point
			}
		}
		return k;
	}
};

Matrix<int> find_opt_PM::further_column_permutation_iteration;