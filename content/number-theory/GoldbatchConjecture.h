/**
 * Author: Borworntat D.
 * Date: 2023-12-02
 * License: CC0
 * Source: https://www.geeksforgeeks.org/program-for-goldbachs-conjecture-two-primes-with-given-sum/ 
 * Description: Find two prime numbers which sum equals $s$ 
 * Time: $O(N\log{N})$
 * Status: Tested
 */

#include "FastEratosthenes.h"

pair<int, int> goldbatchConjecture(int s, vector<int> p={}) {
	if(s <= 2 || s % 2 != 0) {
		return make_pair(-1, -1);
	}
	if(p.size() == 0) {
		p = eratosthenes();
	}
	for(auto x: p) {
		if(x > s / 2) {
			break;
		}
		int d = s - x;
		if(binary_search(p.begin(), p.end(), d)) {
			return make_pair(min(x, d), max(x, d));
		}
	}
	return make_pair(-1, -1);
}
