#pragma once
/*****************************************************************//**
 * \file   Matrix_common.h
 * \brief  static function related to Matrix, stored inside a class
 * 
 * \author 26259
 * \date   July 2023
 *********************************************************************/

#include "Matrix.h"

class Matrix_common {
public:

	/**
	 * .take an example
	 *
	 * input n = 4, k = 2, we have 4 choose 2 equals 6
	 * return a Matrix<int> of size (6*2):
			 0             1
			 0             2
			 0             3
			 1             2
			 1             3
			 2             3
	 * generate the pattern reversely, we get a good OSD TEP ordering
	 */
	static Matrix<int> generating_all_n_choose_k_pattern(int n, int k) {
		int total_row = my::n_choose_k(n, k);
		Matrix<int> ans(total_row, k);
		// give the first row
		for (int j = 0; j < k; ++j) {
			ans(0, j) = j;
		}

		for (int i = 1; i < total_row; ++i) {
			// copy the last row
			for (int j = 0; j < k; ++j) {
				ans(i, j) = ans(i - 1, j);
			}

			// adding one, counting carry bit
			int q = k - 1;
			int decreasing_add_in = 0;
			while (q >= 0) {
				ans(i, q) += 1;
				if (ans(i, q) == n - decreasing_add_in) {
					// adding one to the former bit and this bit to be changed
				}
				else {
					break;
				}
				--q;
				decreasing_add_in++;
			}

			q++;
			// we must have q >= 0
			while (q < k) {
				ans(i, q) = ans(i, q - 1) + 1;
				q++;
			}
		}
		return ans;
	}

	/**
	 * .reverse the result of 'generating_all_n_choose_k_pattern', take an example
	 *
	 * input n = 4, k = 2, we have 4 choose 2 equals 6
	 * return a Matrix<int> of size (6*2):
			 2             3
			 1             3
			 1             2
			 0             3
			 0             2
			 0             1
	 * this is for a good OSD TEP ordering
	 */
	static Matrix<int> generating_all_n_choose_k_pattern_rev(int n, int k) {
		int total_row = my::n_choose_k(n, k);
		Matrix<int> ans(total_row, k);
		// give the first row
		for (int j = 0; j < k; ++j) {
			ans(0, j) = n - k + j;
		}

		for (int i = 1; i < total_row; ++i) {
			// copy the last row
			for (int j = 0; j < k; ++j) {
				ans(i, j) = ans(i - 1, j);
			}

			// minus one, counting carry bit
			int q = k - 1;
			while (q > 0) {
				if (ans(i, q) - 1 == ans(i, q - 1)) {
					// minus one to the former bit and this bit to be changed
				}
				else {
					break;
				}
				--q;
			}

			ans(i, q)--;

			q++;
			// we must have q >= 0
			while (q < k) {
				ans(i, q) = n - k + q;
				++q;
			}

			// the test error pattern generated in the loop, one by one and no storage needed

			// for OSD algorithm, consider write a function that accept the last line and generate the next line, with in place changing
		}
		return ans;
	}

	static void OSD_next_TEP(Matrix<int>& TEP_now, int n) {
		int k = TEP_now.size();

		// minus one, counting carry bit
		int q = k - 1;
		while (q > 0) {
			if (TEP_now(q) - 1 == TEP_now(q - 1)) {
				// minus one to the former bit and this bit to be changed
			}
			else {
				break;
			}
			--q;
		}

		TEP_now(q)--;

		q++;
		// we must have q >= 0
		while (q < k) {
			TEP_now(q) = n - k + q;
			++q;
		}
	}

	/**
	 * . generate random Matrix with k elements in (0,1,...,n-1}
	 *
	 * \param n: number of total elements
	 * \param k: number of choosen elements, k<n
	 * \return random Matrix
	 */
	static Matrix<int> n_randomly_choose_k(int n, int k) {
		Matrix<int> natual(1, n, 'N');
		Matrix<int> ans = natual.get_random_element(k);
		return ans;
	}

	/**
	 * .
	 *
	 * \param n: the number to carry-over
	 * \param len: number's length
	 * \return a matrix contain all n-ary number of length len
	 */
	static Matrix<int> n_ary_incremental_scan(int n, int len) {
		int nn = (int)pow(n, len);
		Matrix<int> ans(nn, len);

		for (int j = 0; j < len; ++j) {
			ans(0, j) = 0;
		}

		for (int i = 1; i < nn; ++i) {

			// copy the last column
			for (int j = 0; j < len; ++j) {
				ans(i, j) = ans(i - 1, j);
			}

			// the last digit + 1
			ans(i, len - 1)++;

			// carry over

			for (int last_co = len - 1; ans(i, last_co) == n; last_co--) {
				ans(i, last_co) = 0;
				ans(i, last_co - 1)++;
			}
		}

		return ans;
	}

	/**
	 * .divide [0,1,...,n-1] randomly into 2 parts, one with size k and the other with size n-k
	 *
	 * \param n: the number indicate the size of a natrual matrix
	 * \param k: one of the division part size, one part with size k, the other with size n-k
	 * \param k_part: the part with size k
	 * \param n_res: the part with size n-k
	 */
	static void n_randomly_divide_into_k_and_res(int n, int k, Matrix<int>& k_part, Matrix<int>& n_res) {
		if (k <= n) {
			if (k <= n - k) {
				n_res.resize(1, n, false);
				for (int i = 0; i < n; ++i) {
					n_res(i) = i;
				}
				k_part.resize(1, k, false);
				k_part.resize(1, 0, false);
				for (int i = 0; i < k; ++i) {
					int choose_ind = my::rand_int_adv(0, n - 1);
					k_part.push_back(n_res(choose_ind));
					n_res.switch_ele(choose_ind, n - 1);
					n_res.pop_back();
					n--;
				}
			}
			else {
				k_part.resize(1, n, false);
				for (int i = 0; i < n; ++i) {
					k_part(i) = i;
				}
				n_res.resize(1, n - k, false);
				n_res.resize(1, 0, false);
				for (int i = 0, imax = n - k; i < imax; ++i) {
					int choose_ind = my::rand_int_adv(0, n - 1);
					n_res.push_back(k_part(choose_ind));
					k_part.switch_ele(choose_ind, n - 1);
					k_part.pop_back();
					n--;
				}
			}
		}
		else {
			cout << "(n_randomly_divide_into_k_and_res) warning: n < k" << endl;
			cout << "n = " << n << ", k = " << k << endl;
		}

	}

	/**
	 * . compute sum of the elements of matrix
	 */
	template<class T>
	static T sum(const Matrix<T>& A) {
		T ans = 0;
		for (int i = 0, imax = A.size(); i < imax; ++i) {
			ans = ans + A(i);
		}
		return ans;
	}

	/**
	 * . compute average of the elements of matrix, return double only
	 */
	template<class T>
	static double ave(const Matrix<T>& A) {
		double ans = (double)sum(A);
		if (A.size() != 0) {
			return ans / A.size();
		}
		else {
			cout << "warning: size of Matrix is 0, Matrix_common::ave failed" << endl;
			return 0;
		}
	}

	/**
	 * . get max of the elements of matrix
	 */
	template<class T>
	static T max(const Matrix<T>& A) {
		return A(A.max_ele_ind());
	}

	/**
	 * .
	 *
	 * \param n: number to divide into 2 parts, i.e., 2  = 1 + 1, 3 = 1 + 2
	 * \return all possible dividion on n
	 */
	Matrix<int> partition_n_in_2_parts(int n) {
		Matrix<int> ans(n / 2, 2);
		for (int i = 0, imax = n / 2; i < imax; ++i) {
			ans(i, 0) = i;
			ans(i, 1) = n - i;
		}
		return ans;			// we don't need this
	}
};

