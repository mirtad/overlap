#define _LCP_CPP_

#include <cstdlib>

#include "lcp.hpp"

void lcp(int *LCP, int *SA, int *s, int n) {
	int* POS = (int*) malloc(n * sizeof(int));

	for (int i = 0; i < n; i++) POS[SA[i]] = i;

	int lcp = LCP[0] = 0;

	for (int i = 0; i < n; i++) {
		int pos = POS[i];

		if (pos > 0) {
			int j = SA[pos-1];

			while (s[i+lcp] == s[j+lcp]) lcp++;

			LCP[pos] = lcp;

			if (lcp > 0) lcp--;
		}
	}

	free(POS);
}
