#include <math.h>

// optimizes the order of the fft
// p is the maximum prime allowed in n
void fft_optimize(int *n, int p, int par)
{
	if(*n <= 2) return;

  for(;; (*n)++){
		int i, n2;

		if((par > 0) && (*n % 2 != par)) continue;

		n2 = *n;
		for(i = 2; i<=n2; i++){
			if(n2 % i == 0){
				if(i > p) break;
				n2 = n2 / i; 
				i--; 
			}
		}
		if(i > n2) return;
	}
}
