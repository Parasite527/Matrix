#pragma once
#include <iostream>
#include <vector>
#include <string>
#include <cmath>

class BigInteger {
public:
    explicit BigInteger(const bool b) {
        _size = 1;
        sign = false;
        number.push_back(b ? 1 : 0);
    }

    BigInteger(const int d) {
        _size = 1;
        if(d < 0) {
            sign = true;
            number.push_back(static_cast<unsigned> (-d % base));
            if(static_cast<size_t>(-d) >= base) {
                number.push_back(static_cast<unsigned> (-d / base));
                _size++;
            }
        }
        else {
            sign = false;
            number.push_back(static_cast<unsigned> (d % base));
            if(static_cast<size_t>(d) >= base) {
                number.push_back(static_cast<unsigned> (d / base));
                _size++;
            }
        }
    }

    explicit BigInteger(const long long d) {
        long long tmp = d;
        _size = 0;
        if(tmp < 0) {
            sign = true;
            tmp *= -1;
        }
        else {
            sign = false;
        }
        while(tmp > 0){
            number.push_back(static_cast<unsigned>(tmp % static_cast<long long>(base)));
            tmp /= base;
            _size++;
        }
    }

    explicit BigInteger(const unsigned d) {
        sign = false;
        _size = 1;
        number.push_back(d % base);
        if(static_cast<size_t>(d) >= base) {
            number.push_back(d / base);
            _size++;
        }

    }

    explicit BigInteger(const unsigned long long d) {
        unsigned tmp = d;
        _size = 0;
        sign = false;
        while(tmp > 0) {
            number.push_back(static_cast<unsigned>(tmp % static_cast<unsigned long long>(base)));
            tmp /= base;
            _size++;
        }
    }

    BigInteger():number(), sign(false), _size(0) {}

    [[nodiscard]] std::string toString() const {
        std::string s, s1;
        if(sign)
            s += '-';
        s += std::to_string(number[_size - 1]);
        for(size_t i = 1; i < _size; i++) {
            s1 = std::to_string(number[_size - 1 - i]);
            for(size_t j = s1.size(); j < 9; j++) {
                s += '0';
            }
            s += s1;
        }

        return s;
    }

    void negative() {
        sign = !sign;
    }

    BigInteger& operator=(const BigInteger& a) = default;

    BigInteger& operator=(const int d) {
        number.clear();
        _size = 1;
        if(d < 0) {
            sign = true;
            number.push_back(static_cast<unsigned> (-d % base));
            if(static_cast<size_t>(-d) >= base) {
                number.push_back(static_cast<unsigned> (-d / base));
                _size++;
            }
        }
        else {
            sign = false;
            number.push_back(static_cast<unsigned> (d % base));
            if(static_cast<size_t>(d) >= base) {
                number.push_back(static_cast<unsigned> (d / base));
                _size++;
            }
        }

        return *this;
    }

    BigInteger& operator=(const unsigned d) {
        number.clear();
        sign = false;
        _size = 1;
        number.push_back(d % base);
        if(static_cast<size_t>(d) >= base) {
            number.push_back(d / base);
            _size++;
        }

        return *this;
    }

    BigInteger& operator=(const unsigned long long d1) {
        number.clear();
        _size = 0;
        unsigned long long d = d1;
        sign = false;
        while(d > 0){
            number.push_back(static_cast<unsigned>(d % static_cast<unsigned long long>(base)));
            d /= base;
            _size++;
        }

        return *this;
    }

    BigInteger& operator=(const long long a) {
        number.clear();
        long long d = a;
        _size = 0;
        if(d < 0) {
            sign = true;
            d *= -1;
        }
        else {
            sign = false;
        }
        while(d > 0){
            number.push_back(static_cast<unsigned>(d % static_cast<long long>(base)));
            d /= base;
            _size++;
        }

        return *this;
    }

    BigInteger& operator++() {
        if(sign) {
            size_t i = 0;
            for( ; number[i] == 0; i++){
                number[i] = base - 1;
            }
            number[i]--;
            delete_zero();
        } else {
            size_t i = 0;

            number.push_back(0);
            _size++;

            for( ; number[i] == base - 1; i++){
                number[i] = 0;

            }
            number[i]++;
            delete_zero();
        }

        return *this;
    }

    BigInteger& operator--() {
        if(_size == 1 && number[0] == 0){
            sign = true;
            number[0] = 1;
        } else if(!sign) {
            size_t i = 0;
            for( ; number[i] == 0; i++){
                number[i] = base - 1;
            }
            number[i]--;
            delete_zero();
        } else {
            size_t i = 0;

            number.push_back(0);
            _size++;

            for( ; number[i] == base - 1; i++){
                number[i] = 0;

            }
            number[i]++;
            delete_zero();
        }

        return *this;
    }

    const BigInteger operator++(int) {
        BigInteger tmp(*this);
        if(sign) {
            size_t i = 0;
            for( ; number[i] == 0; i++){
                number[i] = base - 1;
            }
            number[i]--;
            delete_zero();
        } else {
            size_t i = 0;

            number.push_back(0);
            _size++;

            for( ; number[i] == base - 1; i++){
                number[i] = 0;

            }
            number[i]++;
            delete_zero();
        }

        return tmp;
    }

    const BigInteger operator--(int) {
        BigInteger tmp(*this);
        if(_size == 1 && number[0] == 0){
            sign = true;
            number[0] = 1;
        } else if(!sign) {
            size_t i = 0;
            for( ; number[i] == 0; i++){
                number[i] = base - 1;
            }
            number[i]--;
            delete_zero();
        } else {
            size_t i = 0;

            number.push_back(0);
            _size++;

            for( ; number[i] == base - 1; i++){
                number[i] = 0;

            }
            number[i]++;
            delete_zero();
        }

        return tmp;
    }

    BigInteger& operator+=(const BigInteger& a) {
        if(sign == a.sign) {
            unsigned carry = 0;
            for (size_t i = 0; i < std::max(_size, a._size); i++){
                if(i == _size){
                    number.push_back(0);
                    _size++;
                }
                number[i] += carry + (i < a._size ? a.number[i] : 0);
                carry = number[i] / base;
                number[i] %= base;
            }

            if(carry == 1){
                _size++;
                number.push_back(1);
            }

            return *this;
        } else if((sign ? -*this : *this) > (a.sign ? -a : a)){

            unsigned carry = 0;
            for (size_t i = 0; i < std::max(_size, a._size); i++) {
                number[i] +=  base;
                number[i] -= (carry + (i < a._size ? a.number[i] : 0));
                carry = 1 - number[i] / base;
                number[i] %= base;
            }

            delete_zero();

        } else {
            BigInteger tmp;
            tmp = *this;
            *this = a;
            sign = a.sign;
            unsigned carry = 0;
            for (size_t i = 0; i < std::max(_size, tmp._size); i++) {
                number[i] +=  base;
                number[i] -= (carry + (i < tmp._size ? tmp.number[i] : 0));
                carry = 1 - number[i] / base;
                number[i] %= base;
            }

            delete_zero();
        }

        return *this;
    }

    BigInteger& operator-=(const BigInteger& a) {
        *this += -a;
        return *this;
    }

    BigInteger& operator*=(const BigInteger& a) {

        sign = (a.sign && !sign) || (!a.sign && sign);

        std::vector<unsigned long long> res(a._size + _size);
        for (size_t i = 0; i < _size; ++i) {
            for (size_t j = 0; j < a._size; ++j) {
                res[i + j] += static_cast<unsigned long long>(number[i]) * static_cast<unsigned long long>(a.number[j]);
                res[i + j + 1] += res[i + j] / static_cast<unsigned long long>(base);
                res[i + j] %= static_cast<unsigned long long>(base);
            }
        }

        _size = a._size + _size;

        number.resize(_size);
        for(size_t i = 0; i < _size; ++i) {
            number[i] = static_cast<unsigned>(res[i]);
        }

        delete_zero();

        return *this;
    }

    BigInteger& operator/=(const BigInteger& a) {
        if(_size < a._size) { //if this.size < a.size => return 0
            sign = false;
            _size = 1;
            number.clear();
            number.push_back(0);

            return *this;
        } else if (_size == a._size) {
            bool q_sign = (a.sign && !sign) || (!a.sign && sign);
            sign = false;
            BigInteger quotient, divider = a;
            divider.sign = false;
            unsigned l = 0, r = base;

            while (r - l > 1) {
                unsigned m = (l + r) / 2;
                if (*this < divider * static_cast<BigInteger>(m))
                    r = m;
                else
                    l = m;
            }

            quotient.number.push_back(l);
            quotient._size++;
            sign = q_sign;
            number.clear();
            _size = quotient._size;

            for(size_t j = 0; j < quotient._size; j++) {
                number.push_back(quotient.number[quotient._size - 1 - j]);
            }

            delete_zero();

            return *this;
        }

        bool q_sign = (a.sign && !sign) || (!a.sign && sign);
        sign = false;
        BigInteger quotient, divider = a, carry;
        divider.sign = false;
        size_t i = divider._size;

        for (size_t j = 0; j < divider._size; j++) {
            carry.number.push_back(number[_size - divider.size() + j]);
            carry._size++;
        }

        while (i < _size) {
            //copy digit u need to
            for (size_t k = 0; (carry < divider) && i < _size; k++) {
                carry <<= 1;
                carry += static_cast<BigInteger>(number[_size - i - 1]);

                if (k > 0) {
                    quotient.number.push_back(0);
                    quotient._size++;
                }

                i++;

                if(carry._size == 1 && carry.number[0] == 0) {
                    carry.number.clear();
                    carry._size = 0;
                }
            }

            unsigned l = 0, r = base;
            while (r - l > 1) {
                unsigned m = (l + r) / 2;
                if (carry < (divider * static_cast<BigInteger>(m)))
                    r = m;
                else
                    l = m;
            }

            quotient.number.push_back(l);
            quotient._size++;

            carry -= static_cast<BigInteger>(l) * divider;
            if (carry._size == 1 && carry.number[0] == 0) {
                carry.number.clear();
                carry._size = 0;
            }

        }

        sign = q_sign;
        number.clear();
        _size = quotient._size;
        for (size_t j = 0; j < quotient._size; j++) {
            number.push_back(quotient.number[quotient._size - 1 - j]);
        }

        delete_zero();

        return *this;
    }

    BigInteger& operator%=(const BigInteger& a) {
        return *this -= *this / a * a;;
    }

    friend bool operator!=(const BigInteger& a, const BigInteger& b);

    friend bool operator==(const BigInteger& a, const BigInteger& b);

    friend bool operator<(const BigInteger& a, const BigInteger& b);

    friend bool operator<=(const BigInteger& a, const BigInteger& b);

    friend bool operator>(const BigInteger& a, const BigInteger& b);

    friend bool operator>=(const BigInteger& a, const BigInteger& b);

    friend BigInteger operator-(const BigInteger& a);

    friend BigInteger operator-(const BigInteger& a, const BigInteger& b);

    friend BigInteger operator+(const BigInteger& a, const BigInteger& b);

    friend BigInteger operator*(const BigInteger& a, const BigInteger& b);

    friend BigInteger operator/(const BigInteger& a, const BigInteger& b);

    explicit operator bool() {
        return !(_size == 1 && number[0] == 0);
    }

    friend std::ostream& operator<<(std::ostream &out, const BigInteger& a);

    friend std::istream& operator>>(std::istream &in, BigInteger& a);

    friend BigInteger operator"" _bi(const char*);

    friend class Rational;

    [[nodiscard]] size_t size() const {
        return _size;
    }

    [[nodiscard]] bool is_negative() const {
        return sign;
    }

    BigInteger& operator<<=(const unsigned n) {
        BigInteger production;

        for(size_t i = 0; i < n; i++) {
            production.number.push_back(0);
            production._size++;
        }

        for(size_t i = 0; i < _size; i++) {
            production.number.push_back(number[i]);
            production._size++;
        }
        production.sign = sign;

        *this = production;
        return *this;
    } // fast *= base

    BigInteger& operator>>=(const unsigned n) {
        if(n >= _size) {
            _size = 1;
            sign = false;
            number.clear();
            number.push_back(0);

            return *this;
        }

        for(size_t i = 0; i < _size - n; i++) {
            number[i] = number[i + n];
        }
        for(size_t i = 0; i < n; i++) {
            number.pop_back();
        }

        _size -= n;

        return *this;
    } // fast /= base

private:
    std::vector<unsigned> number;
    bool sign;
    size_t _size;
    static const size_t base = 1000000000;

    void delete_zero() {
        while (_size > 1 && number[_size - 1] == 0){
            number.pop_back();
            _size--;
        }

        if(number[0] == 0 && _size == 1)
            sign = false;
    }

    static int char_to_int(const char c) {
        return static_cast<int>(c) - 48;
    }

    BigInteger str_to_Bigint(std::string s) {
        if(s[0] == '-') {
            sign = true;
            s.erase(0, 1);
        }

        for(size_t i = 0; i < s.size() / 9; i++) {
            unsigned digit = 0;
            for(int j = 0; j < 9; j++) {
                digit = 10 * digit + char_to_int(s[s.size() - 9 * (i + 1) + j]);
            }
            number.push_back(digit);
        }

        if(s.size() % 9 != 0) {
            unsigned digit = 0;
            for (size_t i = 0; i < s.size() % 9; i++) {
                digit = 10 * digit + char_to_int(s[i]);
            }

            number.push_back(digit);
        }

        _size = (s.size() - 1) / 9 + 1;
        return *this;
    }

};

BigInteger operator"" _bi(const char* s) {
    BigInteger a;
    std::string tmp(s);
    return a.str_to_Bigint(tmp);
}

bool operator==(const BigInteger& a, const BigInteger& b) {
    if(a.sign != b.sign || a._size != b._size)
        return false;
    for(size_t i = 0; i < a._size; i++) {
        if (a.number[a._size - i - 1] != b.number[a._size - i - 1])
            return false;
    }
    return true;
}

bool operator!=(const BigInteger& a, const BigInteger& b) {
    return !(a == b);
}

bool operator<(const BigInteger& a, const BigInteger& b) {
    if(a.sign != b.sign)
        return a.sign;
    if(a._size != b._size)
        return (a._size < b._size && !a.sign) || (a._size > b._size && a.sign);

    for(size_t i = 0; i < a._size; i++){
        if(a.number[a._size - i - 1] < b.number[a._size - i - 1])
            return true;
        if(a.number[a._size - i - 1] > b.number[a._size - i - 1])
            return false;
    }

    return false;
}

bool operator<=(const BigInteger& a, const BigInteger& b) {
    return !(b < a);
}

bool operator>(const BigInteger& a, const BigInteger& b) {
    return b < a;
}

bool operator>=(const BigInteger& a, const BigInteger& b) {
    return !(a < b);
}

BigInteger operator+(const BigInteger& a, const BigInteger& b) {
    BigInteger sum(a);
    sum += b;
    return sum;
}

BigInteger operator-(const BigInteger& a, const BigInteger& b) {
    BigInteger sum(a);
    sum -= b;
    return sum;
}

BigInteger operator*(const BigInteger& a, const BigInteger& b) {
    BigInteger composition(a);
    composition *= b;
    return composition;
}

BigInteger operator/(const BigInteger& a, const BigInteger& b) {
    BigInteger quotient(a);
    quotient /= b;
    return quotient;
}

BigInteger operator%(const BigInteger& a, const BigInteger& b) {
    BigInteger remainder(a);
    remainder %= b;
    return remainder;
}

BigInteger operator-(const BigInteger& a) {
    BigInteger tmp;
    tmp = a;
    tmp.sign ^= 1;
    if(tmp._size == 1 && tmp.number[0] == 0)
        tmp.sign = false;
    return tmp;
}

BigInteger operator<<(const BigInteger& a,const unsigned n) {
    BigInteger production(a);
    production <<= n;
    return production;
} /// * base fast

BigInteger operator>>(const BigInteger& a,const unsigned n) {
    BigInteger production(a);
    production >>= n;
    return production;
} /// / base fast

std::ostream& operator<<(std::ostream& out, const BigInteger& a) {
    out << a.toString();
    return out;
}

std::istream& operator>>(std::istream& in, BigInteger& a) {
    a.number.clear();
    a._size = 0;
    a.sign = false;

    std::string tmp;
    in >> tmp;

    a.str_to_Bigint(tmp);
    return in;
}

class Rational {
public:
    Rational(): numerator(), denominator(1) {};

    Rational(const Rational& a) = default;

    Rational(const int a): numerator(a), denominator(1) {}

    explicit Rational(const unsigned a): numerator(a), denominator(1) {}

    explicit Rational(const long long a): numerator(a), denominator(1) {}

    explicit Rational(const unsigned long long a): numerator(a), denominator(1) {}

    Rational(const BigInteger& a): numerator(a), denominator(1) {}

    [[nodiscard]]std::string toString() const {
        this->reduction();

        if(denominator == 1_bi)
            return numerator.toString();

        return numerator.toString() + '/' + denominator.toString();
    } //returns string-fraction

    [[nodiscard]]std::string asDecimal(size_t precision = 0) const {//wrong
        this->reduction();

        std::string s, s1;
        BigInteger tmp(numerator);
        if(tmp.sign) {
            s += '-';
            tmp.negative();
        }

        tmp /= denominator;

        s += tmp.toString();

        size_t trunc_size = s.size();
        tmp = numerator;
        if(tmp.sign)
            tmp.negative();

        tmp %= denominator;
        if(tmp == 0)
            return s;

        s += '.';
        tmp <<= (precision - 1 + 9) / 9;

        tmp /= denominator;

        std::string last_digit = std::to_string(tmp.number[tmp.size() - 1]);
        for(size_t i = last_digit.size() + 9 * (tmp.size() - 1); i < (precision - 1 + 9) / 9 * 9; i++)
            s += '0';
        s += tmp.toString();

        while (s.size() - 1 > trunc_size + precision) {
            s.pop_back();
        }

        if(s[s.size() - 1] == '.')
            s.pop_back();

        if(s.size() == 2 && s[0] == '-' && s[1] == '0')
            s = "0";

        return s;
    }

    explicit operator double() { ///// sure ok
        std::string s;

        this->reduction();

        if(!numerator.sign)
            s = asDecimal(324);
        else{
            numerator.sign = false;
            s = asDecimal(324);
            numerator.sign = true;
        }

        double res = 0;
        size_t dot_pos = s.find('.');

        s.erase(dot_pos, 1);

        for (size_t i = 0; i < s.size() - 1; ++i) {
            res += static_cast<double>(s[i] - 48) * pow(10, static_cast<int>(dot_pos - i - 1));
        }

        if (numerator.sign)
            res *= -1;

        return res;
    }

    BigInteger gcd() const {
        BigInteger num = numerator, denom = denominator;
        if(numerator.sign)
            num.negative();
        while(true) {
            num %= denom;
            if(num == 0)
                return denom;

            denom %= num;
            if(denom == 0) {
                return num;
            }
        }
    }

    void reduction() const {
        if(denominator.is_negative()){
            denominator.negative();
            numerator.negative();
        }
        BigInteger gcd = this->gcd();
        denominator /= gcd;
        numerator /= gcd;
    }

    Rational& operator+=(const Rational& a) {
        numerator = numerator * a.denominator + a.numerator * denominator;
        denominator *= a.denominator;
        this->reduction();
        return *this;
    }

    Rational& operator-=(const Rational& a){
        numerator = numerator * a.denominator - a.numerator * denominator;
        denominator *= a.denominator;
        this->reduction();
        return *this;
    }

    Rational& operator*=(const Rational& a){
        numerator *= a.numerator;
        denominator *= a.denominator;
        //this->reduction();
        return *this;
    }

    Rational& operator/=(Rational a){
        numerator *= a.denominator;
        denominator *= a.numerator;
        this->reduction();
        return *this;
    }

    Rational& operator=(const int a) {
        Rational tmp(a);
        std::swap(numerator, tmp.numerator);
        std::swap(denominator, tmp.denominator);
        return *this;
    }

    Rational& operator=(const unsigned a) {
        Rational tmp(a);
        std::swap(numerator, tmp.numerator);
        std::swap(denominator, tmp.denominator);
        return *this;
    }

    Rational& operator=(const long long a) {
        Rational tmp(a);
        std::swap(numerator, tmp.numerator);
        std::swap(denominator, tmp.denominator);
        return *this;
    }

    Rational& operator=(const unsigned long long a) {
        Rational tmp(a);
        std::swap(numerator, tmp.numerator);
        std::swap(denominator, tmp.denominator);
        return *this;
    }

    Rational& operator=(const Rational& a) = default;

    Rational& operator=(const BigInteger& a) {
        Rational tmp(a);
        std::swap(numerator, tmp.numerator);
        std::swap(denominator, tmp.denominator);
        return *this;
    }

    friend bool operator==(const Rational& a, const Rational& b);

    friend bool operator!=(const Rational& a, const Rational& b);

    friend bool operator<=(const Rational& a, const Rational& b);

    friend bool operator>=(const Rational& a, const Rational& b);

    friend bool operator<(const Rational& a, const Rational& b);

    friend bool operator>(const Rational& a, const Rational& b);

    friend bool operator==(const Rational& a, const Rational& b);

    friend Rational operator-(const Rational& a);

    friend Rational operator"" _rd(const char* s);

    friend Rational operator"" _frac(const char* s);

    friend std::istream& operator>>(std::istream& in, Rational& a);

    friend std::ostream& operator<<(std::ostream& out, const Rational& a);
private:
    mutable BigInteger numerator;
    mutable BigInteger denominator;

    Rational strDecimal_to_Rational(std::string& s) {
        unsigned point_pos = s.find('.');
        if(point_pos == static_cast<unsigned>(-1))
            point_pos = s.find(',');

        if(point_pos == static_cast<unsigned>(-1)) {
            numerator.str_to_Bigint(s);
            return *this;
        }

        if (s[0] == '-') {
            numerator.sign = true;
            s.erase(0, 1);
        }

        s.erase(point_pos, 1);
        unsigned denominator_size = s.size() - point_pos;

        numerator.str_to_Bigint(s);

        denominator <<= denominator_size / 9;
        for (size_t i = 0; i < denominator_size % 9; i++) {
            denominator *= 10;
        }

        this->reduction();

        return *this;
    }

    Rational strFrac_to_Rational(std::string& num, std::string& den) {
        denominator.number.clear();
        denominator._size = 0;
        denominator.sign = false;

        numerator.str_to_Bigint(num);
        denominator.str_to_Bigint(den);

        this->reduction();

        return *this;
    }

    void clear() {
        numerator.number.clear();
        numerator._size = 0;
        numerator.sign = false;
        denominator = 1_bi;
    }
};

Rational operator+(const Rational& a, const Rational& b) {
    Rational tmp(a);
    tmp += b;
    return tmp;
}

Rational operator-(const Rational& a, const Rational& b) {
    Rational tmp(a);
    tmp -= b;
    return tmp;
}

Rational operator*(const Rational& a, const Rational& b) {
    Rational tmp(a);
    tmp *= b;
    return tmp;
}

Rational operator/(const Rational& a, const Rational& b) {
    Rational tmp(a);
    tmp /= b;
    return tmp;
}

Rational operator-(const Rational& a) {
    Rational tmp(a);
    tmp.numerator.negative();
    return tmp;
}

bool operator==(const Rational& a, const Rational& b) {
    a.reduction();
    b.reduction();
    return a.numerator == b.numerator && a.denominator == b.denominator;
}

bool operator!=(const Rational& a, const Rational& b) {
    return !(a == b);
}

bool operator<=(const Rational& a, const Rational& b) {
    return !(b < a);
}

bool operator>=(const Rational& a, const Rational& b) {
    return !(a < b);
}

bool operator<(const Rational& a, const Rational& b) {
    if(a.numerator > 0_bi && b.numerator < 0_bi)
        return false;
    if(a.numerator < 0_bi && b.numerator > 0_bi)
        return true;

    return a.numerator * b.denominator < b.numerator * a.denominator;
}

bool operator>(const Rational& a, const Rational& b) {
    return b < a;
}

Rational operator"" _rd(const char* s) {
    std::string tmp(s);
    Rational a;
    a.strDecimal_to_Rational(tmp);
    return a;
} // returns Rational from decimal with point

Rational operator"" _frac(const char* s) {
    std::string tmp(s), den;
    Rational a;
    unsigned frac_pos = tmp.find('/');

    if(frac_pos == static_cast<unsigned>(-1)) {
        a.strDecimal_to_Rational(tmp);
        return a;
    }

    size_t size = tmp.size();

    for(size_t i = frac_pos + 1; i < tmp.size(); i++){
        den.push_back(tmp[i]);
    }

    for(size_t i = frac_pos; i < size; i++){
        tmp.pop_back();
    }

    a.strFrac_to_Rational(tmp, den);

    return a;
} // returns Rational from fraction

std::ostream& operator<<(std::ostream& out, const Rational& a){
    out << a.toString();
    return out;
}

std::istream& operator>>(std::istream& in, Rational& a) {
    std::string tmp;
    in >> tmp;
    a.clear();

    if(tmp.find('/') != static_cast<unsigned>(-1)) {
        std::string den;
        unsigned frac_pos = tmp.find('/');

        if(frac_pos == static_cast<unsigned>(-1)) {
            a.strDecimal_to_Rational(tmp);
            return in;
        }

        for(size_t i = frac_pos + 1; i < tmp.size(); i++){
            den.push_back(tmp[i]);
        }

        size_t size = tmp.size();

        for(size_t i = frac_pos; i < size; i++){
            tmp.pop_back();
        }

        a.strFrac_to_Rational(tmp, den);

        return in;
    }

    a.strDecimal_to_Rational(tmp);
    a.reduction();

    return in;
}

template <size_t N>
class Residue {
public:
    Residue(): number(0) {};

    Residue(int a): number(static_cast<size_t>((a % static_cast<int>(N) + static_cast<int>(N))) % N) {}

    Residue(size_t a): number(a % N) {}

    Residue& operator=(const int a) {
        number = static_cast<size_t>((a % static_cast<int>(N) + static_cast<int>(N))) % N;
        return *this;
    }

    Residue& operator=(const Residue<N> a) {
        number = a.number;
        return *this;
    }

    ~Residue() = default;

    Residue& operator+=(const Residue& a) {
        number = (number + a.number + N) % N;
        return *this;
    }

    Residue& operator-=(const Residue& a) {
        number = (number - a.number + N) % N;
        return *this;
    }

    Residue& operator*=(const Residue& a) {
        number = (number * a.number) % N;
        return *this;
    }

    Residue& operator^=(size_t indicator) {
        Residue<N> a(1);

        while(indicator > 1) {
            if(indicator & 1) {
                a *= *this;
                indicator--;
            } else {
                *this *= *this;
                indicator >>= 1;
            }
        }

        return *this *= a;
    }

    Residue& operator/=(const Residue& a) {
        if(Prime<N>::is_prime) {
            *this *= (a ^ (N - 2));
        }
        return *this;
    }

    explicit operator int() const {
        return static_cast<int>(number);
    }

    explicit operator size_t() const {
        return number;
    }

    bool operator==(const Residue<N>& a) const {
        return a.number == number;
    }

    bool operator!=(const Residue<N>& a) const {
        return a.number != number;
    }
private:
    size_t number;

    template <size_t A, size_t B, bool exit>
    struct Is_prime {
        static const bool is_prime = (A == 2) || (A == 3) || ((A % B != 0) && Is_prime<A, B + 1, B * B >= A>::is_prime);
    };

    template <size_t A>
    struct Is_prime<A, A, true> {
        static const bool is_prime = Is_prime<A, A, false>::is_prime;
    };

    template <size_t A, size_t B>
    struct Is_prime<A, B, true> {
        static const bool is_prime = true;
    };

    template <size_t A>
    struct Prime {
        static const bool is_prime = Is_prime<static_cast<size_t>(Is_prime<A, 2, false>::is_prime), 0u, true>::is_prime;
    };
};

template <size_t N>
Residue<N> operator+(const Residue<N>& a, const Residue<N>& b) {
    Residue<N> sum(a);
    sum += b;
    return sum;
}

template <size_t N>
Residue<N> operator-(const Residue<N>& a, const Residue<N>& b) {
    Residue<N> difference(a);
   difference -= b;
    return difference;
}

template <size_t N>
Residue<N> operator*(const Residue<N>& a, const Residue<N>& b) {
    Residue<N> product(a);
    product *= b;
    return product;
}

template <size_t N>
Residue<N> operator/(const Residue<N>& a, const Residue<N>& b) {
    Residue<N> division(a);
    division /= b;
    return division;
}

template <size_t N>
Residue<N> operator^(const Residue<N>& a, size_t b) {
    Residue<N> res(a);
    res ^= b;
    return res;
}

template <size_t N>
std::istream& operator>>(std::istream &in, Residue<N>& a) {
    int tmp;
    std::cin >> tmp;
    a = tmp;
    return in;
}

template <size_t N>
std::ostream& operator<<(std::ostream& out, const Residue<N>& a) {
    out << (size_t)a;
    return out;
}


template <size_t N, size_t M, typename Field = Rational>
class Matrix{
public:
    Matrix(): matrix(N, std::vector<Field>(M, static_cast<Field>(0))) {
        if(N == M) {
            for(size_t i = 0; i < N; i++){
                matrix[i][i] = static_cast<Field>(1);
            }
        }
    }

    Matrix(Field n): matrix(N, std::vector<Field>(M, static_cast<Field>(0))) {
        if(N == M) {
            for(size_t i = 0; i < N; i++){
                matrix[i][i] = n;
            }
        }
    }

    Matrix(const std::vector<std::vector<Field>> matrix): matrix(matrix) {};

    Matrix(std::initializer_list<std::vector<Field>> list): matrix(list) {};

    Matrix(const Matrix& m) = default;

    Matrix& operator=(const Matrix& m) = default;

    Matrix& operator=(const std::vector<std::vector<Field>>& m) {
        matrix = m;
        return *this;
    }

    Matrix<N, M, Field>& operator+=(const Matrix<N, M, Field>& m) {
        for(size_t i = 0; i < N; i++) {
            for(size_t j = 0; j < M; j++) {
                matrix[i][j] += m[i][j];
            }
        }
        return *this;
    }

    Matrix<N, M, Field>& operator-=(const Matrix<N, M, Field>& m) {
        for(size_t i = 0; i < N; i++) {
            for(size_t j = 0; j < M; j++) {
                matrix[i][j] -= m[i][j];
            }
        }
        return *this;
    }

    Matrix operator*=(Field a) {
        for(size_t i = 0; i < N; i++) {
            for(size_t j = 0; j < M; j++) {
                matrix[i][j] *= a;
            }
        }
        return *this;
    }

    template <size_t K>
    Matrix<N, K, Field> operator*(const Matrix<M, K, Field>& m) const {
        Matrix<N, K, Field> product(static_cast<Field>(0));
        for(size_t i = 0; i < N; i++) {
            for(size_t k = 0; k < K; k++) {
                for(size_t j = 0; j < M; j++) {
                    product[i][k] += matrix[i][j] * m[j][k];
                }
            }
        }
        return product;
    }

    Matrix<N, N, Field>& operator*=(const Matrix<N, N, Field>& m) {
        Matrix<N, N, Field> product(static_cast<Field>(0));
        for(size_t i = 0; i < N; i++) {

            for(size_t k = 0; k < N; k++) {
                for(size_t j = 0; j < N; j++) {
                    product[i][k] += matrix[i][j] * m[j][k];
                }
            }
        }
        *this = product;

        return *this;
    }

    Matrix<N, N, Field>& operator^=(size_t indicator) {
        Matrix<N, N,Field> E;

        while(indicator > 1) {
            if(indicator & 1) {
                E *= *this;
                indicator--;
            } else {
                *this *= *this;
                indicator >>= 1;
            }
        }

        return *this *= E;
    }

    Matrix<N, M, Field> operator-(){
        return *this * static_cast<Matrix<N, M, Field>>(-1);
    }

    std::vector<Field> getRow(size_t k) const {
        return matrix[k];
    }

    std::vector<Field> getColumn(size_t k) const {
        std::vector<Field> tmp(N);
        for(size_t i = 0; i < N; i++){
            tmp[i] = matrix[i][k];
        }
        return tmp;
    }

    Field det() const {
        Field deter = static_cast<Field>(1);
        Matrix<N, N, Field> m(matrix);
        for(size_t i = 0; i < N - 1; i++) {
            if(m.matrix[i][i] == static_cast<Field>(0)) {
                for(size_t j = i + 1; j < N; j++) {//
                    if (m.matrix[j][i] != static_cast<Field>(0)){
                        std::swap(m.matrix[i], m.matrix[j]);
                        deter *= static_cast<Field>(-1);
                        break;
                    }
                }
            }
            if(m.matrix[i][i] == static_cast<Field>(0)) {
                return static_cast<Field>(0);
            }

            for(size_t q = i + 1; q < N; q++) {//
                Field k = m.matrix[q][i] / m.matrix[i][i];
                if(m.matrix[q][i] != static_cast<Field>(0)) {
                    for (size_t j = i; j < N; j++) {
                        m.matrix[q][j] -= m.matrix[i][j] * k;
                    }
                }
            }
        }

        for(size_t i = 0; i < N; i++) {
            deter *= m.matrix[i][i];
        }

        return deter;
    }

    [[nodiscard]] size_t rank() const {
        Matrix m(matrix);
        for(size_t i = 0; i < N - 1; i++) {
            if(m.matrix[i][i] == static_cast<Field>(0)) {
                for(size_t j = i + 1; j < N; j++) {
                    if (m.matrix[j][i] != static_cast<Field>(0)){
                        std::swap(m.matrix[i], m.matrix[j]);
                        break;
                    }
                }
            }
            if(m.matrix[i][i] == static_cast<Field>(0)) {
                continue;
            }


            for(size_t q = i + 1; q < N; q++) {
                Field k = m.matrix[q][i] / m.matrix[i][i];
                if(m.matrix[q][i] != static_cast<Field>(0)) {
                    for (size_t j = i; j < M; j++) {
                        m.matrix[q][j] -= m.matrix[i][j] * k;
                    }
                }
            }
        }

        size_t rank = 0;
        for(size_t i = 0; i < N; i++) {
            for(size_t j = 0; j < M; j++) {
                if(m.matrix[i][j] != static_cast<Field>(0)) {
                    rank++;
                    break;
                }
            }
        }
        return rank;
    }

    Field trace() const {
        if(equal<N, M>::is_equal) {
            Field trace;
            for (size_t i = 0; i < N; i++) {
                trace += matrix[i][i];
            }
            return trace;
        }
    }

    Matrix<M, N, Field> transposed() const {
        Matrix<M, N, Field> trans;
        for(size_t i = 0; i < M; i++) {
            trans[i] = getColumn(i);
        }
        return trans;
    }

    Matrix<N, N, Field> inverted() const {
        Matrix<N, N, Field> invert, copy(*this);

        for(size_t i = 0; i < N - 1; i++) {
            if(copy.matrix[i][i] == static_cast<Field>(0)) {
                for(size_t j = i + 1; j < N; j++) {
                    if (copy.matrix[j][i] != static_cast<Field>(0)){
                        std::swap(copy.matrix[i], copy.matrix[j]);
                        std::swap(invert.matrix[i], invert.matrix[j]);
                        break;
                    }
                }
            }
            if(copy.matrix[i][i] == static_cast<Field>(0)) {
                continue;
            }

            for(size_t q = i + 1; q < N; q++) {
                Field k = copy.matrix[q][i] / copy.matrix[i][i];
                if(copy.matrix[q][i] != static_cast<Field>(0)) {
                    for (size_t j = 0; j < N; j++) {
                        copy.matrix[q][j] -= copy.matrix[i][j] * k;
                        invert.matrix[q][j] -= invert.matrix[i][j] * k;
                    }
                }
            }
        }


/////////////////////////////////////////////////////////////////////////////

        for(int i = static_cast<int>(N) - 1; i > 0; i--) {
            if(copy.matrix[i][i] == static_cast<Field>(0)) {
                for(int j = i - 1; j >= 0; j--) {
                    if (copy.matrix[j][i] != static_cast<Field>(0)){
                        std::swap(copy.matrix[i], copy.matrix[j]);
                        std::swap(invert.matrix[i], invert.matrix[j]);
                        break;
                    }
                }
            }

            if(copy.matrix[i][i] == static_cast<Field>(0)) {
                continue;
            }

            for(int q = i - 1; q >= 0; q--) {
                Field k = copy.matrix[q][i] / copy.matrix[i][i];
                if(copy.matrix[q][i] != static_cast<Field>(0)) {
                    for (int j = N - 1; j >= 0; j--) {//////////////////////////
                        copy.matrix[q][j] -= copy.matrix[i][j] * k;/////////////
                        invert.matrix[q][j] -= invert.matrix[i][j] * k;

                    }
                }

            }
        }


        for(size_t i = 0; i < N; i++) {
            Field k = copy.matrix[i][i];
            for (size_t j = 0; j < N; j++) {
                invert.matrix[i][j] /= k;
            }
        }

        return invert;
    }

    Matrix<N, N, Field>& invert() {
        *this = inverted();
        return *this;
    }

    std::vector<Field>& operator[](size_t Row) {
        return matrix[Row];
    }

   const std::vector<Field>& operator[](size_t Row) const {
        return matrix[Row];
    }

    template<size_t U, size_t V, typename F>
    friend std::ostream& operator<<(std::ostream& out, const Matrix<U, V, F>& m);

    template<size_t U, size_t V, typename F>
    friend std::istream& operator>>(std::istream& in, Matrix<U, V, F>& m);

    template<size_t U, size_t V, typename F>
    friend bool operator==(const Matrix<N, M, Field>& m1, const Matrix<N, M, Field>& m2);

    template<size_t U, size_t V, typename F>
    friend bool operator!=(const Matrix<N, M, Field>& m1, const Matrix<N, M, Field>& m2);
private:
    std::vector<std::vector<Field>> matrix;

    template<size_t V, size_t U>
    struct equal {
        static const bool is_equal = equal<V, U>::is_equal;
    };

    template<size_t V>
    struct equal<V, V> {
        static const bool is_equal = true;
    };
};

template <size_t N, size_t M, typename Field = Rational>
Matrix<N, M, Field> operator+(const Matrix<N, M, Field>& m1, const Matrix<N, M, Field>& m2) {
    Matrix sum(m1);
    sum += m2;
    return sum;
}

template <size_t N, size_t M, typename Field = Rational>
Matrix<N, M, Field> operator-(const Matrix<N, M, Field>& m1, const Matrix<N, M, Field>& m2) {
    Matrix difference(m1);
    difference -= m2;
    return difference;
}

template <size_t N, size_t M, typename Field = Rational>
Matrix<N, M, Field> operator*(const Matrix<N, M, Field>& m1, Field a) {
    Matrix product(m1);
    product *= a;
    return product;
}

template <size_t N, size_t M, typename Field = Rational>
Matrix<N, M, Field> operator*(Field a, const Matrix<N, M, Field>& m1) {
    Matrix product(m1);
    product *= a;
    return product;
}

template<size_t N, typename Field = Rational>
Matrix<N, N, Field> operator^(const Matrix<N, N, Field>& m, size_t indicator) {
    Matrix res(m);
    res ^= indicator;
    return res;
}

template <size_t N, size_t M, typename Field = Rational>
bool operator==(const Matrix<N, M, Field>& m1, const Matrix<N, M, Field>& m2){
    for(size_t i = 0; i < N; i++){
        for(size_t j = 0; j < M; j++){
            if(m1[i][j] != m2[i][j]){
                return false;
            }
        }
    }
    return true;
}

template <size_t N, size_t M, typename Field = Rational>
bool operator!=(const Matrix<N, M, Field>& m1, const Matrix<N, M, Field>& m2){
    return !(m1 == m2);
}

template <size_t N, size_t M, typename Field = Rational>
std::ostream& operator<<(std::ostream& out, const Matrix<N, M, Field>& m) {
    for(size_t i = 0; i < N; i++){
        for(size_t j = 0; j < M; j++){
            out << m.matrix[i][j] << ' ';
        }
        out << '\n';
    }
    return out;
}

template<size_t N, size_t M, typename Field = Rational>
std::istream& operator>>(std::istream& in, Matrix<N, M, Field>& m) {
    for(size_t i = 0; i < N; i++) {
        for (size_t j = 0; j < M; j++) {
            in >> m.matrix[i][j];
        }
    }
    return in;
}

template <size_t N, typename Field = Rational>
using SquareMatrix = Matrix<N, N, Field>;