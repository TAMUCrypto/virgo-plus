//
// Created by 69029 on 8/24/2021.
//

#ifndef VIRGO_PLUS_TRANSCRIPTCACHE_HPP
#define VIRGO_PLUS_TRANSCRIPTCACHE_HPP
#define MEMPOOL_CAPACITY 100000
#include <cstring>

extern "C"{
#include "../lib/libXKCP.a.headers/SimpleFIPS202.h"
}

class transcriptCache {
public:
    transcriptCache() { mempool_sz = 0; }

    template<class T>
    void store(const T *in, size_t in_sz) {
        assert(mempool_sz + in_sz <= MEMPOOL_CAPACITY);
        memcpy(mempool + mempool_sz, reinterpret_cast<const unsigned char *>(in), in_sz);
        mempool_sz += in_sz;
    }

    template<class T>
    void store(const T &in) {
        size_t in_sz = sizeof (in);
        assert(mempool_sz + in_sz <= MEMPOOL_CAPACITY);
        memcpy(mempool + mempool_sz, reinterpret_cast<const unsigned char *>(&in), in_sz);
        mempool_sz += in_sz;
    }

    void load(unsigned char *out) {
        assert(mempool_sz);
        SHA3_256(out, mempool, mempool_sz);
        mempool_sz = 32;
        memcpy(mempool, out, 32);
    }

    virgo::fieldElement random() {
        static unsigned char out[32];
        load(out);
        u64 real = *reinterpret_cast<u64*>(out) % virgo::fieldElement::mod;
        u64 imag = *reinterpret_cast<u64*>(out + 8) % virgo::fieldElement::mod;
        return virgo::fieldElement((i64) real, (i64) imag);
    }
private:
    unsigned char mempool[MEMPOOL_CAPACITY]{};
    size_t mempool_sz;
};


#endif //VIRGO_PLUS_TRANSCRIPTCACHE_HPP
