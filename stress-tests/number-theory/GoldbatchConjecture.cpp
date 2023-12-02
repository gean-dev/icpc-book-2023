#include "../utilities/template.h"

#include "../../content/number-theory/FastEratosthenes.h"
#include "../../content/number-theory/GoldbatchConjecture.h"

int main(int argc, char* argv[]) {
	pair<int, int> fail = make_pair(-1, -1);
	int limit = 1000;
	vector<int> p = eratosthenes();
	for(int i=4; i<=limit; i+=2) {
		pair<int, int> r = goldbatchConjecture(i, p);
		assert(r != fail && r.first + r.second == i && isPrime[r.first] && isPrime[r.second]);
	}
	cout << "Tests passed!" << endl;
}
