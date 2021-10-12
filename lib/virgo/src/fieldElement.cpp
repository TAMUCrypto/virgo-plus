#include "fieldElement.hpp"
#include <cstdio>
#include <chrono>
#include <immintrin.h>

namespace virgo {
    const unsigned long long fieldElement::mod = 2305843009213693951LL;
    const unsigned long long fieldElementPacked::mod = 2305843009213693951LL;
    __m256i fieldElementPacked::packed_mod, fieldElementPacked::packed_mod_minus_one;
    bool fieldElement::initialized = false;
    int fieldElement::multCounter, fieldElement::addCounter;
    bool fieldElement::isCounting;
    bool fieldElement::isSumchecking;

    fieldElement::fieldElement() {
        real = img = 0;
    }

    fieldElement::fieldElement(const fieldElement &b) {
        real = b.real;
        img = b.img;
    }

    fieldElement::fieldElement(long long x) {
        real = x >= 0 ? x : mod + x;
        img = 0;
    }

    fieldElement::fieldElement(long long int x, long long int y) {
        real = x >= 0 ? x : mod + x;
        img = y >= 0 ? y : mod + y;
    }

    fieldElement fieldElement::operator+(const fieldElement &other) const {
        if (isCounting)
        {
            ++addCounter;
        }
        fieldElement ret;
        ret.img = other.img + img;
        ret.real = other.real + real;
        if (mod <= ret.img)
            ret.img = ret.img - mod;
        if (mod <= ret.real)
            ret.real = ret.real - mod;
        return ret;
    }

    fieldElement fieldElement::operator*(const fieldElement &other) const {
        if (isCounting)
        {
            ++multCounter;
        }
        fieldElement ret;
        auto all_prod = mymult(img + real, other.img + other.real); //at most 6 * mod
        //unsigned long long ac, bd;
        //mymult_2vec(real, b.real, img, b.img, ac, bd);
        auto ac = mymult(real, other.real), bd = mymult(img, other.img); //at most 1.x * mod
        auto nac = ac;
        if (bd >= mod)
            bd -= mod;
        if (nac >= mod)
            nac -= mod;
        nac ^= mod; //negate
        bd ^= mod; //negate

        auto t_img = all_prod + nac + bd; //at most 8 * mod
        t_img = myMod(t_img);
        if (t_img >= mod)
            t_img -= mod;
        ret.img = t_img;
        auto t_real = ac + bd;

        while (t_real >= mod)
            t_real -= mod;
        ret.real = t_real;
        return ret;
    }

    fieldElement fieldElement::operator-(const fieldElement &other) const {
        if (isCounting)
        {
            ++addCounter;
        }
        fieldElement ret;
        auto tmp_r = other.real ^ mod; //tmp_r == -b.real is true in this prime field
        auto tmp_i = other.img ^ mod; //same as above
        ret.real = real + tmp_r;
        ret.img = img + tmp_i;
        if (ret.real >= mod)
            ret.real -= mod;
        if (ret.img >= mod)
            ret.img -= mod;

        return ret;
    }

    fieldElement fieldElement::operator-() const {
        if (isCounting)
        {
            ++addCounter;
        }
        return zero() - *this;
    }

    void fieldElement::init() {
        initialized = true;
        srand(3396);
        isCounting = false;
        fieldElementPacked::init();
    }

    void fieldElementPacked::init() {
        packed_mod = _mm256_set_epi64x(mod, mod, mod, mod);
        packed_mod_minus_one = _mm256_set_epi64x(mod - 1, mod - 1, mod - 1, mod - 1);
    }


    fieldElement fieldElement::random() {
        fieldElement ret;
        ret.real = randomNumber() % mod;
        ret.img = randomNumber() % mod;
        return ret;
    }

    bool fieldElement::operator!=(const fieldElement &other) const {
        return real != other.real || img != other.img;
    }

    bool fieldElement::operator==(const fieldElement &other) const {
        return !(*this != other);
    }

    fieldElement &fieldElement::operator=(const fieldElement &other) {
        real = other.real;
        img = other.img;
        return *this;
    }

    fieldElement &fieldElement::operator+=(const fieldElement &other) {
        *this = *this + other;
        return *this;
    }

    fieldElement &fieldElement::operator-=(const fieldElement &other) {
        *this = *this - other;
        return *this;
    }

    fieldElement &fieldElement::operator*=(const fieldElement &other) {
        *this = *this * other;
        return *this;
    }

    fieldElement::operator bool() const {
        return real || img;
    }

    bool fieldElement::isNegative() const {
        return ((real >> 60) & 1) && img == 0;
    }

    unsigned char fieldElement::getBitWidth() const {
        assert(img == 0);
        auto dat = real;
        unsigned char res = 0;
        for (int i = 32; i && dat; i >>= 1) {
            if (dat >> i) {
                res += i;
                dat >>= i;
            }
        }
        return res + 1;
    }

    unsigned char fieldElement::getBit(unsigned int i) const {
        assert(img == 0);
        return (real >> i) & 1;
    }

    __int128_t fieldElement::toint128() const {
        assert(img == 0);
        if ((real >> 60) & 1) return i64(real) - mod;
        return real;
    }

    bool fieldElement::operator<(const fieldElement &other) const {
        assert(img == 0);
        return real < other.real;
    }

    bool fieldElement::isZero() {
        return !real && !img;
    }

    fieldElement fieldElement::abs() const {
        assert(img == 0);
        fieldElement res = -*this;
        return res.real < this -> real ? res : real;
    }

    fieldElement fieldElement::sqr() const {
        return (*this) * (*this);
    }

    fieldElement fieldElement::inv() const {
        __uint128_t p = 2305843009213693951LL;
        return fastPow(*this, p * p - 2);
    }

    void fieldElement::setAbs() {
        *this = this -> abs();
    }

    void fieldElement::setSqr() {
        *this = this -> sqr();
    }

    void fieldElement::setInv() {
        *this = this -> inv();
    }

    void fieldElement::print(FILE *fileno) const {
        fprintf(fileno, "(%llu %llu)\n", real, img);
    }

    fieldElement fieldElement::maxWithZero(const fieldElement &a, const fieldElement &b) {
        assert(a.img == 0 && b.img == 0);
        if (a.isNegative() && b.isNegative()) return fieldElement::zero();
        return fieldElement(a.isNegative() ? b : b.isNegative() ? a : std::max(a.real, b.real));
    }

    fieldElement fieldElement::maxUnsigned(const fieldElement &a, const fieldElement &b) {
        return a < b ? b : a;
    }

    fieldElement fieldElement::getRootOfUnity(int log_order) {
        fieldElement rou;
        //general root of unity, have log_order 2^61
        rou.img = 1033321771269002680L;
        rou.real = 2147483648L;

        assert(log_order <= 61);

        for (int i = 0; i < __max_order - log_order; ++i)
            rou = rou * rou;

        return rou;
    }

    fieldElement fieldElement::zero() {
        return fieldElement(0);
    }

    fieldElement fieldElement::one() {
        return fieldElement(1);
    }

    vector<fieldElement> fieldElement::generateRandomness(u64 size) {
        int k = size;
        vector<fieldElement> ret(k);

        for (int i = 0; i < k; ++i)
            ret[i] = fieldElement::random();
        return ret;
    }

    fieldElement fieldElement::innerProd(vector<fieldElement>::iterator a, vector<fieldElement>::iterator b, u64 n) {
        fieldElement res = fieldElement::zero();
        for (int i = 0; i < n; ++i)
            res += a[i] * b[i];
        return res;
    }

    double fieldElement::self_speed_test_add(int repeat) {
        fieldElement a, b;
        a = random();
        b = random();
        fieldElement c;
        auto t0 = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < repeat; ++i)
        {
            c = a + b;
            b = c + a;
            a = c + b;
        }
        printf("%llu %llu\n", c.img, c.real);
        auto t1 = std::chrono::high_resolution_clock::now();
        return std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0).count();
    }

    double fieldElement::self_speed_test_mult(int repeat) {
        fieldElement a, b;
        a = random();
        b = random();
        fieldElement c;
        auto t0 = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < repeat; ++i)
        {
            c = a * b;
            b = c * a;
            a = c * b;
        }
        printf("%llu %llu\n", c.img, c.real);
        auto t1 = std::chrono::high_resolution_clock::now();
        return std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0).count();
    }

    char *fieldElement::toString() const {
        char *s = new char[50];
        sprintf(s, "(%llu %llu)", this -> real, this -> img);
//    if (isNegative()) sprintf(s, "-%llu", (-(*this)).real);
//    else sprintf(s, "%llu", real);
        return s;
    }

    ostream &operator << (ostream &out, const fieldElement &c) {
        out << '(' << c.real << ' ' << c.img << ')';
        return out;
    }

    fieldElement fieldElement::fastPow(fieldElement x, __uint128_t p) {
        fieldElement ret = fieldElement(1), tmp = x;
        while (p)
        {
            if (p & 1)
            {
                ret = ret * tmp;
            }
            tmp = tmp * tmp;
            p >>= 1;
        }
        return ret;
    }

    unsigned long long fieldElement::myMod(unsigned long long int x) {
        return (x >> 61) + (x & mod);
    }

    unsigned long long fieldElement::mymult(const unsigned long long int x, const unsigned long long int y) {
        //return a value between [0, 2PRIME) = x * y mod PRIME
        /*
        unsigned long long lo, hi;
        lo = _mulx_u64(x, y, &hi);
        return ((hi << 3) | (lo >> 61)) + (lo & PRIME);
        */
        unsigned long long hi;
        asm(
        "mov %[x_read], %%rdx;\n"
        "mulx %[y_read], %%r9, %%r10;"
        "shld $0x3, %%r9, %%r10;\n"
        "and %[mod_read], %%r9;\n"
        "add %%r10, %%r9;\n"
        "mov %%r9, %[hi_write]"
        : [hi_write] "=r"(hi)
        : [x_read] "r"(x), [y_read]"r"(y), [mod_read]"r"(mod)
        : "rdx", "r9", "r10"
        );
        return hi;
    }

    unsigned long long fieldElement::randomNumber() {
        unsigned long long ret = ::random() % 10;
        for (int i = 1; i < 20; ++i)
            ret = (ret * 10ull + (unsigned long long)(::random() % 10)) % mod;
        return ret;
    }


    fieldElementPacked::fieldElementPacked() {
        real = _mm256_set_epi64x(0, 0, 0, 0);
        img = _mm256_set_epi64x(0, 0, 0, 0);
    }

    fieldElementPacked::fieldElementPacked(const fieldElement &x0, const fieldElement &x1, const fieldElement &x2,
                                           const fieldElement &x3) {
        real = _mm256_set_epi64x(x3.real, x2.real, x1.real, x0.real);
        img = _mm256_set_epi64x(x3.img, x2.img, x1.img, x0.img);
    }

    fieldElementPacked fieldElementPacked::operator+(const fieldElementPacked &b) const {
        fieldElementPacked ret;
        ret.img = b.img + img;
        ret.real = b.real + real;
        __m256i msk0, msk1;
        msk0 = _mm256_cmpgt_epi64(ret.img, packed_mod_minus_one);
        msk1 = _mm256_cmpgt_epi64(ret.real, packed_mod_minus_one);
        ret.img = ret.img - _mm256_and_si256(msk0, packed_mod);
        ret.real = ret.real - _mm256_and_si256(msk1, packed_mod);
        return ret;
    }

    fieldElementPacked fieldElementPacked::operator*(const fieldElementPacked &b) const {
        fieldElementPacked ret;
        __m256i all_prod = packed_mymult(img + real, b.img + b.real); //at most 6 * mod
        __m256i ac = packed_mymult(real, b.real), bd = packed_mymult(img, b.img); //at most 1.x * mod
        __m256i nac = ac;
        __m256i msk;
        msk = _mm256_cmpgt_epi64(bd, packed_mod_minus_one);
        bd = _mm256_sub_epi64(bd, _mm256_and_si256(packed_mod, msk));

        msk = _mm256_cmpgt_epi64(nac, packed_mod_minus_one);
        nac = _mm256_sub_epi64(nac, _mm256_and_si256(packed_mod, msk));

        nac = _mm256_xor_si256(nac, packed_mod);
        bd = _mm256_xor_si256(bd, packed_mod);

        __m256i t_img = _mm256_add_epi64(_mm256_add_epi64(all_prod, nac), bd);
        t_img = packed_myMod(t_img);

        msk = _mm256_cmpgt_epi64(t_img, packed_mod_minus_one);
        t_img = _mm256_sub_epi64(t_img, _mm256_and_si256(packed_mod, msk));

        ret.img = t_img;
        __m256i t_real = _mm256_add_epi64(ac, bd);
        while (1)
        {
            msk = _mm256_cmpgt_epi64(t_real, packed_mod_minus_one);
            int res = _mm256_testz_si256(msk, msk);
            if (res)
                break;
            t_real = _mm256_sub_epi64(t_real, _mm256_and_si256(packed_mod, msk));
        }

        ret.real = t_real;
        return ret;
    }

    fieldElementPacked fieldElementPacked::operator-(const fieldElementPacked &b) const {
        fieldElementPacked ret;
        __m256i tmp_r = b.real ^ packed_mod; //tmp_r == -b.real is true in this prime field
        __m256i tmp_i = b.img ^ packed_mod; //same as above
        ret.real = real + tmp_r;
        ret.img = img + tmp_i;
        __m256i msk0, msk1;
        msk0 = _mm256_cmpgt_epi64(ret.real, packed_mod_minus_one);
        msk1 = _mm256_cmpgt_epi64(ret.img, packed_mod_minus_one);

        ret.real = ret.real - _mm256_and_si256(msk0, packed_mod);
        ret.img = ret.img - _mm256_and_si256(msk1, packed_mod);

        return ret;
    }

    void fieldElementPacked::getFieldElement(fieldElement *dst) const {
        static unsigned long long real_arr[packed_size], img_arr[packed_size];
        _mm256_store_si256((__m256i *)real_arr, real);
        _mm256_store_si256((__m256i *)img_arr, img);
        for (int i = 0; i < 4; ++i)
        {
            dst[i].real = real_arr[i];
            dst[i].img = img_arr[i];
        }
    }

    __mmask8 fieldElementPacked::operator == (const fieldElementPacked &b) const {
        __m256i res_real = real ^ b.real;
        __m256i res_img = img ^ b.img;
        return _mm256_testz_si256(res_real, res_real) && _mm256_testz_si256(res_img, res_img);
    }

    __mmask8 fieldElementPacked::operator != (const fieldElementPacked &b) const {
        return !(*this == b);
    }

    __m256i fieldElementPacked::packed_mymult(const __m256i x, const __m256i y) {
        __m256i ac, ad, bc, bd;
        __m256i x_shift, y_shift;
        x_shift = _mm256_srli_epi64(x, 32);
        y_shift = _mm256_srli_epi64(y, 32);
        bd = _mm256_mul_epu32(x, y);
        ac = _mm256_mul_epu32(x_shift, y_shift);
        ad = _mm256_mul_epu32(x_shift, y);
        bc = _mm256_mul_epu32(x, y_shift);

        __m256i ad_bc = _mm256_add_epi64(ad, bc);
        __m256i bd_srl32 = _mm256_srli_epi64(bd, 32);
        __m256i ad_bc_srl32 = _mm256_srli_epi64(_mm256_add_epi64(ad_bc, bd_srl32), 32);
        __m256i ad_bc_sll32 = _mm256_slli_epi64(ad_bc, 32);
        __m256i hi = _mm256_add_epi64(ac, ad_bc_srl32);

        __m256i lo = _mm256_add_epi64(bd, ad_bc_sll32);


        //return ((hi << 3) | (lo >> 61)) + (lo & PRIME);
        return _mm256_add_epi64(_mm256_or_si256(_mm256_slli_epi64(hi, 3), _mm256_srli_epi64(lo, 61)), _mm256_and_si256(lo, packed_mod));
    }

    __m256i fieldElementPacked::packed_myMod(const __m256i x) {
        //return (x >> 61) + (x & mod);
        __m256i srl64 = _mm256_srli_epi64(x, 61);
        __m256i and64 = _mm256_and_si256(x, packed_mod);
        return _mm256_add_epi64(srl64, and64);
    }
}