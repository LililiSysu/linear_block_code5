/*****************************************************************//**
 * \file   Complex.h
 * \brief  class Complex, must include "my.h", a seperate file can help manage code
 * 
 * \todo   the compare of complex number should be fixed to only compare the real part
 *		   remove the 'isFA_form' of the complex, or write a new simple class
 *		   this is recommendation for simplicity and speed of the program
 * 
 * \author lilili
 * \date   October 2022
 *********************************************************************/
#pragma once

#include"my.h"		// must include
#include<math.h>
#include<iostream>

// Note: don't write a destructor function unless necessary

#define Complex_class
// Complex class
class Complex
{
	// the angle range from (-my::pi, my::pi]
public:
	double x;
	double y;
	bool isFA_form; //0: Real-Imaginary. 1: Fabs-Angle.

	//default: Real-Imaginary type.
	Complex() = default;
	Complex(double _x, double _y = 0, bool _isFA_form = false) \
		:x(_x), y(_y), isFA_form(_isFA_form) {}

	static const Complex J;
	static double ang_abs_in_pi(double ang) {
		if (ang > -my::pi && ang <= my::pi);
		else if (ang > my::pi) {
			do {
				ang -= 2 * my::pi;
			} while (ang > my::pi);
		}
		else{
			do {
				ang += 2 * my::pi;
			} while (ang <= -my::pi);
		}
		return ang;
	};
	double get_real() const {
		if (!isFA_form)
			return x;
		else
			return x * cos(y);
	}
	double get_imag() const {
		if (!isFA_form)
			return y;
		else
			return x * sin(y);
	}
	double get_abs() const {
		if (!isFA_form)
			return sqrt(x * x + y * y);
		else
			return my::abs(x);
	}
	double get_angle() const {
		if (!isFA_form)
			if (x > 0)
				return atan(y / x);
			else if (x < 0 && y >= 0)
				return atan(y / x) + my::pi;
			else if (x < 0 && y < 0)
				return atan(y / x) - my::pi;
			else if (x == 0 && y != 0)
				return (y > 0) * my::pi / 2;
			else
				return 0;	//set the angle of (0 + 0j) be 0
		else
			return y;
	}
	void set_RI() {
		if (isFA_form) {
			double x_store = x;
			x = x_store * cos(y);
			y = x_store * sin(y);
			isFA_form = false;
		}
	}
	void set_FA() {
		if (!isFA_form) {
			double y_store = y;
			y = get_angle();
			x = sqrt(x * x + y_store * y_store);
		}
	}
	// if imag over real abs less than my::zero_approximation, approach it as real
	bool is_Real() {
		return my::abs(get_imag() / get_real()) < my::zero_approximation;		
	}

	explicit operator double() const {
		if (!isFA_form)
			return x;
		else
			return x * cos(y);
	}

	void operator += (const Complex& c2) {
		if (!isFA_form) {
			x += c2.get_real();
			y += c2.get_imag();
		}
		else{
			double tmp = x;
			x = tmp * cos(y) + c2.get_real();
			y = tmp * sin(y) + c2.get_imag();
			isFA_form = false;
			set_FA();
		}
	}
	void operator -= (const Complex& c2) {
		if (!isFA_form) {
			x -= c2.get_real();
			y -= c2.get_imag();
		}
		else {
			double tmp = x;
			x = tmp * cos(y) - c2.get_real();
			y = tmp * sin(y) - c2.get_imag();
			isFA_form = false;
			set_FA();
		}
	}
	void operator *= (const Complex& c2) {
		if (!isFA_form) {
			double c2_real = c2.get_real();
			double c2_imag = c2.get_imag();
			double tmp = x;
			x = tmp * c2_real - y * c2_imag;
			y = tmp * c2_imag + y * c2_real;
		}
		else {
			x *= c2.get_abs();
			y = ang_abs_in_pi(y + c2.get_angle());
		}
	}
	void operator /= (const Complex& c2) {
		if (!isFA_form) {
			double c2_real = c2.get_real();
			double c2_imag = c2.get_imag();
			double temp = c2_real * c2_real + c2_imag * c2_imag;
			double tmp = x;
			x = (tmp * c2_real + y * c2_imag) / temp;
			y = (y * c2_real - tmp * c2_imag) / temp;
		}
		else {
			x /= c2.get_abs();
			y = ang_abs_in_pi(y - c2.get_angle());
		}
	}

	//the following answer form is same with the [c1/c] form.
	friend Complex operator + (const Complex& c1, const Complex& c2) {
		if (!c1.isFA_form)
			return Complex(c1.x + c2.get_real(), c1.y + c2.get_imag());
		else
		{
			Complex ans(c1.x * cos(c1.y) + c2.get_real(), c1.x * sin(c1.y) + c2.get_imag());
			ans.set_FA();
			return ans;
		}
	}
	friend Complex operator - (const Complex& c1, const Complex& c2) {
		if (!c1.isFA_form)
			return Complex(c1.x - c2.get_real(), c1.y - c2.get_imag());
		else {
			Complex ans(c1.x * cos(c1.y) - c2.get_real(), c1.x * sin(c1.y) - c2.get_imag());
			ans.set_FA();
			return ans;
		}
	}
	friend Complex operator * (const Complex& c1, const Complex& c2) {
		if (!c1.isFA_form) {
			double c2_real = c2.get_real();
			double c2_imag = c2.get_imag();
			return Complex(c1.x * c2_real - c1.y * c2_imag,	c1.x * c2_imag + c1.y * c2_real);
		}
		else
			return Complex(c1.x * c2.get_abs(), ang_abs_in_pi(c1.y + c2.get_angle()), true);
	}
	friend Complex operator / (const Complex& c1, const Complex& c2) {
		if (!c1.isFA_form) {
			double c2_real = c2.get_real();
			double c2_imag = c2.get_imag();
			double temp = c2_real * c2_real + c2_imag * c2_imag;
			return Complex((c1.x * c2_real + c1.y * c2_imag) / \
				temp, (c1.y * c2_real - c1.x * c2_imag) / temp);
		}
		else
			return Complex(c1.x / c2.get_abs(), ang_abs_in_pi(c1.y - c2.get_angle()), true);
	}

	friend Complex operator - (const Complex& c) {
		if (!c.isFA_form)
			return Complex(-c.x, -c.y);
		else
			return Complex(c.x, ang_abs_in_pi(c.y + my::pi), true);
	}
	friend bool operator == (const Complex& c1, const Complex& c2) {
		return ((c1.isFA_form == c2.isFA_form) ? (c1.x == c2.x && c1.y == c2.y) : \
			(c1.get_real() == c2.get_real() && c1.get_imag() == c2.get_imag()));
	}
	friend bool operator != (const Complex& c1, const Complex& c2) {
		return !(c1 == c2);
	}

	// ordered with real value, keep same as matlab
	friend bool operator > (const Complex& c1, const Complex& c2) {
		return c1.get_real() > c2.get_real();
	}
	friend bool operator < (const Complex& c1, const Complex& c2) {
		return c1.get_real() < c2.get_real();
	}
	friend bool operator >= (const Complex& c1, const Complex& c2) {
		return c1.get_real() >= c2.get_real();
	}
	friend bool operator <= (const Complex& c1, const Complex& c2) {
		return c1.get_real() <= c2.get_real();
	}

	friend std::ostream& operator << (std::ostream& out, const Complex& c) {
		// Note: approach my::zero_approximation as 0, WARNING: this may make misunderstanding, in displaying small numbers (discarded)
		if (!c.isFA_form) {
#ifdef Complex_output_supress_approx_0
			bool x_non_approx_zero = my::abs(c.x) > my::zero_approximation;
			bool y_non_approx_zero = my::abs(c.y) > my::zero_approximation;
			if (x_non_approx_zero && y_non_approx_zero) {
#endif 
				if (c.y >= 0)
					out << c.x << "+" << c.y << 'j';
				else
					out << c.x << "" << c.y << 'j';
#ifdef Complex_output_supress_approx_0
			}
			else if (my::abs(c.y) > my::zero_approximation) {
				out << c.y << 'j';
			}
			else if (my::abs(c.x) > my::zero_approximation) {
				out << c.x;
			}
			else {
				out << 0;
			}
#endif
		}
		else
			out << c.x << "e^(" << c.y << "j)";
		return out;
	}
	friend std::istream& operator >> (std::istream& in, Complex& c) {
		in >> c.x >> c.y;
		c.isFA_form = false;	// default case is real-imaginary input
		return in;
	}

	// friend comply with <cmath> functions
	friend Complex sqrt(const Complex& c) {
		if (!c.isFA_form) {
			Complex ans(sqrt(c.get_abs()), c.get_angle() / 2, true);
			ans.set_RI();
			return ans;
		}
		else
			return Complex(sqrt(c.x), c.y / 2, true);
	}
	friend Complex exp(const Complex& c) {
		if (!c.isFA_form)
			return Complex(exp(c.x) * cos(c.y), exp(c.x) * sin(c.y));
		else
			return Complex(exp(c.x * cos(c.y)), ang_abs_in_pi(c.x * sin(c.y)), true);
	}
	friend Complex log(const Complex& c) {
		if (!c.isFA_form) {
			Complex ans(log(sqrt(c.x * c.x + c.y * c.y)), c.get_angle());
			return ans;
		}
		else {
			Complex ans(log(c.x), c.y);
			ans.set_FA();
			return ans;
		}
	}
	friend Complex pow(const Complex& c1, const Complex& c2) {
		Complex temp(c2);
		temp.set_RI();
		Complex ans = c1;
		ans.set_FA();
		ans = Complex(pow(ans.x, temp.x), ang_abs_in_pi(ans.y * temp.x), true);
		if (!c1.isFA_form)
			ans.set_RI();
		return ans * exp(temp.y * J * log(c1));
	}
	friend Complex log10(const Complex& c) {
		return 1 / log(10) * log(c);
	}
	friend Complex round(const Complex& c) {
		Complex ans(round(c.get_real()), round(c.get_imag()));
		if (c.isFA_form)
			ans.set_FA();
		return ans;
	}

	friend Complex sin(const Complex& c) {
		return (exp(c * J) - exp(-c * J)) / (2 * J);
	}
	friend Complex cos(const Complex& c) {
		return (exp(c * J) + exp(-c * J)) / 2;
	}
	friend Complex tan(const Complex& c) {
		return (exp(c * J) - exp(-c * J)) / (exp(c * J) + exp(-c * J)) / J;
	}
	friend Complex cot(const Complex& c) {
		return (exp(c * J) + exp(-c * J)) * J / (exp(c * J) - exp(-c * J));
	}
	friend Complex sec(const Complex& c) {
		return 2 / (exp(c * J) + exp(-c * J));
	}
	friend Complex asin(const Complex& c) {
		return log(c * J + sqrt(1 - c * c)) * (-J);
	}
	friend Complex acos(const Complex& c) {
		return log(c + sqrt(1 - c * c) * J) * (-J);
	}
	friend Complex atan(const Complex& c) {
		return log((-c + J) / (c + J)) / (2 * J);
	}
	friend Complex sinh(const Complex& c) {
		return (exp(c) - exp(-c)) / 2;
	}
	friend Complex cosh(const Complex& c) {
		return (exp(c) + exp(-c)) / 2;
	}
	friend Complex tanh(const Complex& c) {
		return (exp(c) + exp(-c)) / (exp(c) - exp(-c));
	}
	friend Complex asinh(const Complex& c) {
		return log(c + sqrt(c * c + 1));
	}
	friend Complex acosh(const Complex& c) {
		return log(c + sqrt(c * c - 1));
	}
	friend Complex atanh(const Complex& c) {
		return log((1 + c) / (1 - c)) / 2;
	}
};
const Complex Complex::J = Complex(0, 1);

// implementation of my definition function
Complex my::abs(const Complex& c) {
	if (!c.isFA_form)
		return sqrt(c.x * c.x + c.y * c.y);
	else
		return my::abs(c.x);
}
Complex my::conj(const Complex& c) {
	return Complex(c.x, ((c.isFA_form && c.y == my::pi) ? c.y : -c.y), c.isFA_form);
}