/* quadlods.cpp
 * Quadratic low-discrepancy sequence
 * 
 * This generates sequences of vectors, up to 6542-dimensional, as follows:
 * Each prime is assigned to a dimension: 2 to the 0th, 3 to the 1st, etc.
 * If the prime is congruent to 1 mod 4, take (sqrt(p)+1)/2, else sqrt(p).
 * Add the quadratic number to an accumulator mod 1.
 * Copy the accumulator and exclusive-or it with 0x...9669699669969669.
 * Assemble a vector of all these accumulators xored with the bit pattern.
 * 
 * The quadratic numbers are approximated by rational numbers whose denominator
 * is at least some specified limit. The exclusive-oring is done in such a way
 * that the result will not exceed the denominator.
 */

#include <iostream>
#include "quadlods.h"

using namespace std;

vector<unsigned short> primes;
quadlods quads;

void initprimes()
{
  int i,j;
  bool prime;
  primes.clear();
  for (i=2;i<65535;i++)
  {
    for (j=0,prime=true;j<primes.size() && prime && (unsigned int)primes[j]*primes[j]<=i;j++)
      if (i%primes[j]==0)
	prime=false;
    if (prime)
      primes.push_back(i);
  }
}

void quadlods::init(int dimensions,double resolution)
{
  if (primes.size()==0)
    initprimes();
}

vector<mpq_class> quadlods::gen()
{
}

int main(int argc,char **argv)
{
  int i;
  quads.init(5,1e10);
  for (i=0;i<100;i++)
  {
    cout<<primes[i]<<' ';
    if (i%10==9)
      cout<<endl;
  }
  return 0;
}
