#pragma once
/*****************************************************************//**
 * \file   BCH_dual_trace.h
 * \brief  utilize the trace function to trun RS dual code into a BCH dual code
 * 
 * \author 26259
 * \date   June 2023
 *********************************************************************/

#include"../GF/GF2e.h"
#include"../my_lib/Matrix.h"
#include<unordered_set>
using namespace std;

struct binary_tree{
	binary_tree* left;
	binary_tree* right;
	binary_tree() :left(NULL), right(NULL) {}
	static void destroy_recur(binary_tree* head) {
		if (head->left != NULL) {
			destroy_recur(head->left);
		}
		if (head->right != NULL) {
			destroy_recur(head->right);
		}
		// after that, head has non child and can be destroyed without memory leak

		//head->left = NULL;
		//head->right = NULL;		// this is actually unnecessary
		delete head;
	}
};

/**
 * .use the systematic generator of super code to produce generator of its subfield subcode
 */
template<int m>
class tree_store_method{
public:

	// set by constructor
	int code_length_minus_1;		// this function is for fixed code length only, please decide it in advance
	int info_length;
	int max_code_num;
	Matrix<GF2e<m>> symbol_coordinate;
	
	// set by 'init'
	binary_tree* head;
	int valid_row_num;
	Matrix<GF2> codewords;
	Matrix<int> encode_index_and_xor_index;
	Matrix<GF2e<m>> super_generator_partial;
	Matrix<GF2> codeword_can;

	tree_store_method(int _code_length, int _info_length , const Matrix<GF2e<m>>& _symbol_coordinate) \
		: code_length_minus_1(_code_length - 1), info_length(_info_length), symbol_coordinate(_symbol_coordinate) {
		
		max_code_num = info_length * symbol_coordinate.size();
		head = NULL;
		valid_row_num = 0;
		codewords.resize(max_code_num, code_length_minus_1 + 1, false);
		encode_index_and_xor_index.resize(max_code_num, code_length_minus_1 + 1, false);
			// xor index is no more than code_length_minus_1

		codeword_can.resize(1, code_length_minus_1 + 1, false);
	}
	~tree_store_method() {
		binary_tree::destroy_recur(head);
	}

	void init(const Matrix<GF2e<m>>& _super_generator_partial) {
		if (head != NULL)
			binary_tree::destroy_recur(head);
		head = new binary_tree();
		valid_row_num = 0;

		super_generator_partial = _super_generator_partial;
	}

	/**
	 * .insert a vector to the tree, store the inserted vector in a leaf
	 * .\return bool variable of whether a new leaf is inserted
	 */
	bool insert(const Matrix<GF2>& v) {
		binary_tree* node = head;
		binary_tree** node_child;
		for (int i = 0; i < code_length_minus_1; ++i) {
			node_child = (v(i) == 0) ? &(node->left) : &(node->right);
			if (*node_child != NULL);
			else {
				// adding new node to tree
				*node_child = new binary_tree();
			}
			node = *node_child;			// the node move forward		
		}

		// for the leaves
		node_child = (v(code_length_minus_1) == 0) ? &(node->left) : &(node->right);
		if (*node_child == NULL) {		// we observed that it is not always true

			// adding new node to tree
			*node_child = new binary_tree();

			// adding vector to codewords
			codewords.set_row(valid_row_num, v);
			valid_row_num++;
			return true;
		}
		else {
			return false;
		}
	}
	/**
	 * . the interface may change
	 */
	void get_generator() {
		GF2e<m> tmp;
		for (int i = 0; i < info_length; ++i) {
			for (int j = 0, jmax = symbol_coordinate.size(); j < jmax; ++j) {
				tmp = symbol_coordinate(j);			// note that this is not set by alpha

				// generate a candidate codeword to add on the tree
				for (int p = 0; p <= code_length_minus_1; ++p) {
					codeword_can(p) = (tmp * super_generator_partial(i, p)).trace();
				}

				// to be upgraded
				insert(codeword_can);
			}
		}
		// a naive method, to be upgraded
		codewords.resize(valid_row_num, code_length_minus_1 + 1, true);		// confusing, but unimportant
	}
};

// please write the code step by step, this is a bin class, never use it
//
/**
 * .use the systematic generator of super code to produce generator of its subfield subcode
 */
template<int m>
class subfield_subcode_sys_partial_generator {
public: 

	// set by constructor
	int n;								// n
	int considered_len;					// a little greater than k'-k, since first n-k bits may be dependent 
	int k_prime_dual;					// n-k'
	int k_dual;							// n-k
	int k_prime_minus_k;				// k'-k
	Matrix<GF2e<m>> symbol_coordinate;	// symbols that fill the whole space of a GF2e<m>, size m-1 because of '1' is excluded
	int symbol_coordinate_size;

	// set by 'init'
	Matrix_flex_col<int> xor_index;		// if the first int is -1, indicating this column in 'generator_info_part' is not occupied
	Matrix<GF2> codeword_can;			// length of considered length
	Matrix<GF2e<m>> super_generator;	// generator of super code, RS dual (n,n-k')

	// set by 'insert'
	Matrix_flex_col<GF2> generator_info_part;	// store the info part of generator, by Matrix flex col
	Matrix<int> xor_index_tmp;					// store each line of xor_index temporally

	// set by 'get_generator'
	Matrix<GF2> generator_ans;

	subfield_subcode_sys_partial_generator(int _n, int _k_prime_dual, int _k_dual, int _external_considered_length, const Matrix<GF2e<m>>& _symbol_coordinate) : n(_n), k_prime_dual(_k_prime_dual), k_dual(_k_dual), considered_len(_k_dual - _k_prime_dual + _external_considered_length), symbol_coordinate(_symbol_coordinate) {

		k_prime_minus_k = k_dual - k_prime_dual;
		symbol_coordinate_size = symbol_coordinate.size();
			// in usual symbol_coordinate.size() equals to m-1
		xor_index.resize(Matrix<int>(1, 1, considered_len, 'd'));
			// xor index is no more than considered_len

		// alocate enough size what so ever
		codeword_can.resize(1, considered_len, false);
		generator_info_part.resize(Matrix<int>(considered_len - 1, -1, 1, 'd'));
		//cout << "generator_info_part" << generator_info_part;

		xor_index_tmp.resize(1, considered_len - 1, false);
		generator_ans.resize(k_prime_minus_k, n, false);
	}
	void init(const Matrix<GF2e<m>>& _super_generator) {
		super_generator = _super_generator;
		xor_index.reset(-1);
		generator_info_part.reset(-1);
	}

	/**
	 * .insert a vector to the tree, store the inserted vector in a leaf
	 * 
	 * \return bool variable of whether a new leaf is inserted
	 */
	bool insert(Matrix<GF2>& v, int encode_index) {
		xor_index_tmp.resize(1, 0, false);
		for (int i = 0; i < considered_len; ++i) {
			if (v(i) == 0);
			else if (xor_index(i, 0) == -1) {
				// this sub space is not occupied, store the vector into this sub space
				for (int j = i + 1; j < considered_len; ++j) {
					generator_info_part(i, j - i - 1) = v(j);
				}
				xor_index(i, 0) = encode_index;
				for (int j = 0, jmax = xor_index_tmp.size(); j < jmax; ++j) {
					xor_index(i, j + 1) = xor_index_tmp(j);
				}
				return true;
			}
			else {
				// xor the ocupied vector in 'generator_info_part', and continue
				v(i) = 0;
				for (int j = i + 1; j < considered_len; ++j) {
					v(j) += generator_info_part(i, j - i - 1);
				}
				xor_index_tmp.push_back(i);
			}
		}
		return false;
	}

	/**
	 * . the interface may change, too complex and not right, discarded this function
	 */
	void get_generator() {
		GF2e<m> tmp; 
		unordered_set<int> unoccupied_ind;
		for (int i = 0; i < k_prime_minus_k; ++i) {
			unoccupied_ind.insert(i);
		}
		int occupied_vectors = 0;
		int dependent_vectors = 0;
		bool flag_continue = true;
		for (int i = 0; i < k_prime_dual && flag_continue; ++i) {
			int sym_pos_ind = i * symbol_coordinate_size;
			for (int j = 0; j < symbol_coordinate_size && flag_continue; ++j) {
				tmp = symbol_coordinate(j);			// note that this is not set by alpha

				// generate a candidate codeword to add on the tree
				xor_index_tmp.resize(1, 0, false);
				for (int p = 0; p < considered_len; ++p) {
					codeword_can(p) = (tmp * super_generator(i, p + k_prime_dual)).trace();
				}

				int insert_i = -10;
				for (int i = 0; i < considered_len; ++i) {
					if (codeword_can(i) == 0);
					else if (xor_index(i, 0) == -1) {
						insert_i = i;
						break;
					}
					else {
						// xor the ocupied vector in 'generator_info_part', and continue
						codeword_can(i) = 0;
						for (int j = i + 1; j < considered_len; ++j) {
							codeword_can(j) += generator_info_part(i, j - i - 1);
						}
						xor_index_tmp.push_back(i);
					}
				}

				if (insert_i != -10) {
					// this sub space is not occupied, store the vector into this sub space
					for (int j = insert_i + 1; j < considered_len; ++j) {
						generator_info_part(insert_i, j - insert_i - 1) = codeword_can(j);
					}
					xor_index(insert_i, 0) = sym_pos_ind + j;
					for (int j = 0, jmax = xor_index_tmp.size(); j < jmax; ++j) {
						xor_index(insert_i, j + 1) = xor_index_tmp(j);
					}
					occupied_vectors++;
					unoccupied_ind.erase(insert_i);
					if (occupied_vectors >= k_prime_minus_k && unoccupied_ind.empty()) {
						// if the first k_prime_minus_k position of xor_index is occupied, break,
						flag_continue = false;
					}
				}
				else if (occupied_vectors >= k_prime_minus_k) {
					dependent_vectors++;
					if (dependent_vectors == 3) {
						flag_continue = false;
					}
				}
			}
		}

		// retrive the generator through xor_index

		// process the first column of xor_index
		int row_ind = 0;
		int row_cnt = 0;
		
		while (row_cnt < k_prime_minus_k && row_ind < considered_len) {
			int xor_index_tmp = xor_index(row_ind, 0);
			if (xor_index_tmp != -1) {
				int sym_ind = xor_index_tmp / symbol_coordinate_size;
				int sym_coo = xor_index_tmp % symbol_coordinate_size;
				tmp = symbol_coordinate(sym_coo);

				for (int p = 0; p < n; ++p) {
					generator_ans(row_ind, p) = (tmp * super_generator(sym_ind, p)).trace();
				}
				row_cnt++;
			}
			row_ind++;
		}
		cout << "xor_index" << xor_index;
		
		if (row_cnt == k_prime_minus_k) {
			row_cnt--; 
			row_ind--;
			while (row_cnt > 0) {
				if (xor_index(row_ind, 0) != -1) {
					int col_ind = 1;
					int xor_index_tmp = xor_index(row_ind, col_ind);
					while (xor_index_tmp != -1) {
						for (int p = 0; p < n; ++p) {
							generator_ans(row_ind, p) += generator_ans(xor_index_tmp, p);
						}
						col_ind++;
						xor_index_tmp = xor_index(row_ind, col_ind);
					}
					row_cnt--;
				}
				row_ind--;
			}
		}
		else{
			cout << "error of generating k_dual generator, for rows not enough" << endl;
			cout << "row_cnt = " << row_cnt << "\tk_dual = " << k_dual << endl;
		}

	}
};

