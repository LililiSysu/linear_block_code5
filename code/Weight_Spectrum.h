/*****************************************************************//**
 * \file   Weight_Spectrum.h
 * \brief  
 * 
 * \author 26259
 * \date   April 2024
 *********************************************************************/

#pragma once

#include"../GF/GF2.h"
#include"../my_lib/Matrix.h"

class Weight_Specturm {
public:

	Matrix<GF2> G_sys;
	int n;
	int k;

	/* store the results */
	Matrix<Matrix<GF2>> codewords;			// the codewords with weight 'i' has index 'i'
	Matrix<unsigned long long> ans;			// the weight specturm corresponding to generator matrix 'G_sys'

	/* store the requirement */
	int w_max_allowed;
	bool codewords_needed;

	/* store the intermediate variables */
	Matrix<GF2> codeword_layer;

	Weight_Specturm(int _n, int _k) {
		n = _n;
		k = _k;

		w_max_allowed = 0;
		codewords_needed = false;
	}

	Weight_Specturm(const Matrix<GF2>& _G_sys) {
		n = _G_sys.col();
		k = _G_sys.row();
		G_sys = _G_sys;
	}

	/**
	 * .generate the weight specturm in 'ans', and codewords with weight <= 'w_max_required' in 'codewords'
	 * 
	 * \param w_max_required: the weight upper bound in 'codewords'
	 */
	void solve(int _w_max_allowed, bool _codewords_needed) {
		w_max_allowed = _w_max_allowed;
		codewords_needed = _codewords_needed;

		codeword_layer.resize(1, n);
		codeword_layer.reset(0);

		codewords.resize(1, w_max_allowed + 1);
		codewords.reset(Matrix<GF2>());				// we must use dynamic push back property, 2 times the storage space

		ans.resize(1, n + 1);
		ans.reset(0);

		solve_recur(-1, 0);

		for (int i = 0; i <= w_max_allowed; ++i) {
			if (codewords(i).size() > 0) {
				codewords(i).resize(codewords(i).size() / n, n);		// let each row be a codeword
			}
		}

		// convert the form of ans into weight - num
		Matrix<unsigned long long> ans_new;
		for (int i = 0; i <= n; ++i) {
			if (ans(i) != 0) {
				ans_new.push_back(i);
				ans_new.push_back(ans(i));
			}
		}
		ans_new.resize(ans_new.size() / 2, 2);
		ans = ans_new;
	}

	void solve_recur(int layer, int layer_weight) {
		if (layer_weight == w_max_allowed + 1) {
			return;
		}

		if (layer == k - 1) {
			// reach the leaf node, generate a valid codeowrd
						
			// check the weight of last layer
			int weight = layer_weight;
			for (int i = k; i < n; ++i) {
				if (codeword_layer(i) == 1) {
					weight++;
				}
			}

			// check if with in the weight constraint
			if (weight <= w_max_allowed) {
				if (codewords_needed == true) {
					for (int i = 0; i < n; ++i) {
						codewords(weight).push_back(codeword_layer(i));
					}
				}
				ans(weight)++;
			}
			return;
		}

		layer++;
		// the next layer is a 0 info bit
		solve_recur(layer, layer_weight);

		// the next layer is a 1 info bit, needs to change the current 'codeword_layer'

		int layer_row_ind = layer * n;

		codeword_layer(layer) = 1;
		for (int i = k; i < n; ++i) {
			codeword_layer(i) += G_sys(layer_row_ind + i);		
		}

		solve_recur(layer, layer_weight + 1);

		// undo the change
		for (int i = k; i < n; ++i) {
			codeword_layer(i) += G_sys(layer_row_ind + i);		// recover the parity-check part
		}
		codeword_layer(layer) = 0;
	}
};


class Weight_Specturm_transform{
public:
	
	/**
	 * .count the weight specturm of code after flipping position, i.e., weight specturm of coset code
	 * 
	 * \param codewords: the codewords classified with weight
	 * \param flip_pos: the positions flipped in the code
	 * \param max_weight_cnt: max considered weight
	 * \return the weight specturm of coset code
	 */
	static Matrix<unsigned long long> solve(const Matrix<Matrix<GF2>>& codewords, const Matrix<int>& flip_pos, int max_weight_cnt) {
		Matrix<unsigned long long> ans(1, max_weight_cnt + 1, '0');

		int flip_pos_size = flip_pos.size();
		int weight_size = codewords.size();
		for (int i = 0; i < weight_size; ++i) {		// for each weight
			int m = codewords(i).row();
			int n = codewords(i).col();
			if (n == 0) {
				continue;
			}
			for (int j = 0; j < m; ++j) {	// for each codeword
				int weight = i;
				int j_row_ind = j * n;
				for (int p = 0; p < flip_pos_size; ++p) {	// for each flip
					int flip = flip_pos(p);
					if (codewords(i)(j_row_ind + flip) == 1) {	// count weight change
						weight--;
					}
					else {
						weight++;
					}
				}
				if (weight <= max_weight_cnt){
					ans(weight)++;
				}
			}
		}
		//cout << "ans" << ans;

		// convert the form of ans into weight - num
		Matrix<unsigned long long> ans_new;
		for (int i = 0; i <= max_weight_cnt; ++i) {
			if (ans(i) != 0) {
				ans_new.push_back(i);
				ans_new.push_back(ans(i));
			}
		}
		ans_new.resize(ans_new.size() / 2, 2);
		ans = ans_new;

		return ans;
	}
};
