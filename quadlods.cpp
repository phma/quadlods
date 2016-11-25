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
/* Copyright 2014,2016 Pierre Abbat.
 * This file is part of the Quadlods library.
 * 
 * The Quadlods library is free software: you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * Quadlods is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License and Lesser General Public License
 * for more details.
 *
 * You should have received a copy of the GNU General Public License
 * and Lesser General Public License along with Quadlods. If not, see
 * <http://www.gnu.org/licenses/>.
 */

#include <iostream>
#include <cfloat>
#include <cmath>
#include "quadlods.h"

using namespace std;

vector<unsigned short> primes;
mpz_class thue(0x69969669),third(0x55555555);
int morse(32),b2adic(32);

quadirr::quadirr()
{
  a=c=0;
  b=d=p=1;
}

quadirr::quadirr(int A,int B,int C,int D,int P)
{
  a=A;
  b=B;
  c=C;
  d=D;
  p=P;
}

double quadirr::realval()
{
  return a/(double)b+c*sqrt(p)/d;
}

mpz_class thuemorse(int n)
{
  while (morse<=n)
  {
    thue+=(mpz_class)(((thue&((mpz_class)1<<(morse>>5)))>0)?(unsigned)0x96696996:0x69969669)<<morse;
    morse+=32;
  }
  return thue&(((mpz_class)1<<n)-1);
}

mpz_class minusthird(int n)
{
  while (b2adic<=n)
  {
    third+=(mpz_class)0x55555555<<b2adic;
    b2adic+=32;
  }
  return third&(((mpz_class)1<<n)-1);
}

mpz_class graydecode(mpz_class n)
{
  int i;
  i=mpz_sizeinbase(n.get_mpz_t(),2);
  while (i&(i-1))
    i+=i&-i;
  while (i)
  {
    n^=n>>i;
    i/=2;
  }
  return n;
}

mpz_class jumble(mpz_class acc,mpz_class denom,int jumbletype)
{
  int i;
  mpz_class bitdiff,hibits,lobits,ret;
  bitdiff=denom&~acc;
  i=mpz_sizeinbase(bitdiff.get_mpz_t(),2)-1;
  switch (jumbletype)
  {
    case QL_JUMBLE_THIRD:
      ret=acc^minusthird(i);
      break;
    case QL_JUMBLE_THUEMORSE:
      ret=acc^thuemorse(i);
      break;
    case QL_JUMBLE_GRAY:
      lobits=acc&(((mpz_class)1<<i)-1);
      hibits=acc-lobits;
      ret=hibits+graydecode(lobits);
      break;
    default:
      ret=acc;
  }
  return ret;
}

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

void compquad(int p,double resolution,mpz_class &nmid,mpz_class &dmid)
{
  mpz_class nhi,dhi,nlo,dlo,comp;
  for (nhi=dlo=1,nlo=dhi=dmid=0;dmid<resolution;)
  {
    dmid=dhi+dlo;
    nmid=nhi+nlo;
    if ((p-1)&3)
      comp=nmid*nmid-p*dmid*dmid;
    else
      comp=nmid*(nmid-dmid)-(p/4)*dmid*dmid;
    if (comp>0)
    {
      dhi=dmid;
      nhi=nmid;
    }
    else
    {
      dlo=dmid;
      nlo=nmid;
    }
  }
  nmid%=dmid;
}

int nthprime(int n)
{
  if (n<0 || n>=primes.size())
    return 0;
  else
    return primes[n];
}

double nthquad(int n)
{
  mpz_class nmid,dmid;
  if (n<0 || n>=primes.size())
    return NAN;
  else
  {
    compquad(primes[n],27/DBL_EPSILON,nmid,dmid);
    return mpq_class(nmid,dmid).get_d();
  }
}

void quadlods::init(int dimensions,double resolution,int j)
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
 * (99/70)²-99/70=99²/70²-99*70/70²
 */
{
  int i,p;
  mpz_class nmid,dmid;
  if (primes.size()==0)
    initprimes();
  if (dimensions>primes.size())
    dimensions=primes.size();
  if (dimensions<0)
    dimensions=0;
  for (i=denom.size();i<dimensions;i++)
  {
    p=primes[i];
    compquad(p,resolution,nmid,dmid);
    denom.push_back(dmid);
    num.push_back(nmid);
  }
  if (j>=0)
    jumbletype=j;
  if (jumbletype<0 || jumbletype>QL_JUMBLE_GRAY)
    jumbletype=QL_JUMBLE_GRAY;
  num.resize(dimensions);
  denom.resize(dimensions);
  acc.resize(dimensions);
}

void quadlods::init(vector<int> dprimes,double resolution,int j)
/* Used for testing, setting up a generator with a particularly
 * bad set of primes.
 */
{
  int i,p;
  mpz_class nmid,dmid;
  if (primes.size()==0)
    initprimes(); // TODO check that dprimes are actually prime and distinct
  for (i=denom.size();i<dprimes.size();i++)
  {
    p=dprimes[i];
    compquad(p,resolution,nmid,dmid);
    denom.push_back(dmid);
    num.push_back(nmid);
  }
  if (j>=0)
    jumbletype=j;
  if (jumbletype<0 || jumbletype>QL_JUMBLE_GRAY)
    jumbletype=QL_JUMBLE_GRAY;
  num.resize(dprimes.size());
  denom.resize(dprimes.size());
  acc.resize(dprimes.size());
}

vector<mpq_class> quadlods::readout()
{
  int i;
  vector<mpq_class> ret;
  for (i=0;i<num.size();i++)
  {
    ret.push_back(mpq_class((jumble(acc[i],denom[i],jumbletype)<<1)|1,denom[i]<<1));
    ret[i].canonicalize();
  }
  return ret;
}

vector<double> quadlods::dreadout()
{
  int i;
  vector<double> ret;
  for (i=0;i<num.size();i++)
  {
    ret.push_back(mpq_class((jumble(acc[i],denom[i],jumbletype)<<1)|1,denom[i]<<1).get_d());
  }
  return ret;
}

void quadlods::setmiddle()
/* Set the point to the middle of the square/cube/etc. This is used in testing
 * to find the successively smaller distances between q[n] and q[n+h].
 */
{
  int i;
  for (i=0;i<num.size();i++)
    acc[i]=denom[i]>>1;
}

void quadlods::advance(mpz_class n)
{
  int i;
  for (i=0;i<num.size();i++)
    if (n<0)
      acc[i]=(acc[i]-n*(denom[i]-num[i]))%denom[i];
    else
      acc[i]=(acc[i]+n*num[i])%denom[i];
}

unsigned int quadlods::seedsize()
{
  unsigned i,maxlen,len;
  for (i=maxlen=0;i<denom.size();i++)
  {
    len=(mpz_sizeinbase(denom[i].get_mpz_t(),2)+7)/8;
    if (len>maxlen)
      maxlen=len;
  }
  return maxlen*denom.size();
}

void quadlods::seed(char *s,unsigned int n)
{
  unsigned i,sz;
  sz=denom.size();
  for (i=0;i<n;i++)
    acc[i%sz]=((acc[i%sz]<<8)+(s[i]&0xff))%denom[i%sz];
}

vector<mpq_class> quadlods::gen()
{
  advance(1);
  return readout();
}

vector<double> quadlods::dgen()
{
  advance(1);
  return dreadout();
}
