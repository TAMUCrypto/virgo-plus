#pragma once

#include <vector>
#include <unordered_map>
#include <utility>
#include <virgo/polyCommit.hpp>
#include <unordered_set>
#include "inputCircuit.hpp"
#include "config_pc.hpp"

class gate {
public:
	gateType ty;
	int l;
    u64 u, v, lv;
    F c;
    bool is_assert;
	gate() { }
	gate(gateType t, int ll, u64 uu, u64 vv, const F &cc, bool is_assert_zero):
	    ty(t), l(ll), u(uu), v(vv), lv(0), c(cc), is_assert(is_assert_zero) {
	}
};


class layer {
public:
	vector<gate> gates;
	int bitLength;
	u64 size;

	vector<vector<u64>> dadId;  // map from subset id to real id
	vector<int> dadBitLength;   // subset bit length
	vector<u64> dadSize;        // subset size
    u64 maxDadSize;             // max subset size
	int maxDadBitLength;        // max subset bit length
};

class layeredCircuit {
public:
	vector<layer> circuit;
	int size;

	static layeredCircuit readFromStream(char *);
	static layeredCircuit randomize(int layer, int eachLayer);
	void subsetInit();
};

