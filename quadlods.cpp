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
/* Sets num[i]/denom[i] to a rational approximation of an integer in Q(sqrt(primes[i])).
 * If p mod 4 is 1, q=(sqrt(p)+1)/2, else q=sqrt(p).
 * 2: p²-2=0
 * 3: p²-3=0
 * 5: p²-p-1=0 p-1=1/p
 * 7: p²-7=0
 * 11: p²-11=0
 * 13: p²-p-3=0 p-1=3/p
 * 17: p²-p-4=0 p-1=4/p
 * Then sets num=num%denom.
 */
{
  int i,p;
  mpz_class n,d,ntmp,dtmp;
  if (primes.size()==0)
    initprimes();
  for (i=denom.size();i<dimensions;i++)
  {
    p=primes[i];
    for (n=d=1;d<resolution;)
    {
      if ((p-1)&3)
      {
	ntmp=d*p+n;
	dtmp=n+d;
      }
      else
      {
	ntmp=(p/4)*d+n;
	dtmp=n;
      }
      n=ntmp;
      d=dtmp;
    }
    n%=d;
    denom.push_back(d);
    num.push_back(n);
  }
  num.resize(dimensions);
  denom.resize(dimensions);
  acc.resize(dimensions);
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
  for (i=0;i<5;i++)
    cout<<quads.num[i]<<'/'<<quads.denom[i]<<' '<<quads.acc[i]<<endl;
  return 0;
}
