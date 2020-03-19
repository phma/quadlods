/******************************************************/
/*                                                    */
/* quadlods.cpp - quadratic low-discrepancy sequence  */
/*                                                    */
/******************************************************/
/* This generates sequences of vectors, up to 6542-dimensional, as follows:
 * Each prime is assigned to a dimension: 5 to the 0th, 3 to the 1st, etc.
 * If the prime is congruent to 1 mod 4, take (sqrt(p)+1)/2, else sqrt(p).
 * Add the quadratic number to an accumulator mod 1.
 * Copy the accumulator and exclusive-or it with 0x...9669699669969669.
 * Assemble a vector of all these accumulators xored with the bit pattern.
 * 
 * The quadratic numbers are approximated by rational numbers whose denominator
 * is at least some specified limit. The exclusive-oring is done in such a way
 * that the result will not exceed the denominator.
 */
/* Copyright 2014,2016-2020 Pierre Abbat.
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
#include <fstream>
#include <cfloat>
#include <cmath>
#include <string>
#include <array>
#include "quadlods.h"
#include "config.h"

#define M_1PHI 0.6180339887498948482046

using namespace std;
using namespace quadlods;

namespace quadlods
{
  map<unsigned,unsigned> relprimes;
  vector<unsigned short> primes;
  map<unsigned,vector<unsigned short> > reverseScrambleTable;
  vector<PrimeContinuedFraction> primesCfSorted;
  mpz_class thue(0x69969669),third(0x55555555);
  int morse(32),b2adic(32);
  int primePowerTable[][2]=
  {
    {16,65536},{10,59049},{8,65536},{6,15625},{6,46656},{5,16807},
    {5,32768},{5,59049},{4,10000},{4,14641},{4,20736},{4,28561}
  };
  array<int,2> primePower(unsigned short p);
  int reverseScramble(int limb,int p,int scrambletype);
  mpz_class thuemorse(int n);
  mpz_class minusthird(int n);
  mpz_class graydecode(mpz_class n);
  mpz_class scramble(mpz_class acc,mpz_class denom,int scrambletype);
  short readshort(std::istream &file);
  void initprimes();
  void compquad(int p,double resolution,mpz_class &nmid,mpz_class &dmid);
  void compquad(ContinuedFraction cf,double resolution,mpz_class &nmid,mpz_class &dmid);
}

unsigned quadlods::gcd(unsigned a,unsigned b)
{
  while (a&&b)
  {
    if (a>b)
    {
      b^=a;
      a^=b;
      b^=a;
    }
    b%=a;
  }
  return a+b;
}

unsigned quadlods::relprime(unsigned n)
// Returns the integer closest to n/φ of those relatively prime to n.
{
  unsigned ret,twice;
  double phin;
  ret=relprimes[n];
  if (!ret)
  {
    phin=n*M_1PHI;
    ret=rint(phin);
    twice=2*ret-(ret>phin);
    while (gcd(ret,n)!=1)
      ret=twice-ret+(ret<=phin);
    relprimes[n]=ret;
  }
  return ret;
}

unsigned quadlods::scrambledig(unsigned dig,unsigned p)
/* Raises dig to the power relprime(p-1) mod p. This is used for scrambling
 * the Halton generator. -1, 0, and 1 are unaffected.
 */
{
  unsigned acc=1,pow=relprime(p-1);
  int i;
  for (i=15;i>=0;i--)
  {
    acc=(acc*acc)%p;
    if ((pow>>i)&1)
      acc=(acc*dig)%p;
  }
  return acc;
}

array<int,2> quadlods::primePower(unsigned short p)
/* Returns the number of digits in base p that fit in one 16-bit limb
 * and the number of different limbs. p should be prime. If p is 0, 1,
 * 14, 15, or 16, returns garbage.
 */
{
  array<int,2> ret;
  if (p<17)
  {
    ret[0]=primePowerTable[p-2][0];
    ret[1]=primePowerTable[p-2][1];
  }
  else
  {
    ret[0]=1;
    ret[1]=p;
    if (p<=256)
    {
      ret[0]++;
      ret[1]*=p;
    }
    if (p<=40)
    {
      ret[0]++;
      ret[1]*=p;
    }
  }
  return ret;
}

void quadlods::fillReverseScrambleTable(int p,int scrambletype)
/* If p=5 and scrambletype=QL_SCRAMBLE_NONE, fills reverseScrambleTable[5].
 * Entry 320041 (base 5) is 140023 (base 5).
 * If p=5 and scrambletype=QL_SCRAMBLE_POWER, fills reverseScrambleTable[0x40005].
 * Entry 320041 (base 5) is 140032 (base 5).
 * Without scrambling, there is no need to fill the table for primes >256.
 */
{
  array<int,2> pp=primePower(p);
  int i,j,dec,acc;
  vector<unsigned short> scrambleTable;
  int inx=(scrambletype<<16)+p;
  if ((scrambletype==QL_SCRAMBLE_POWER || p<256) && reverseScrambleTable[inx].size()==0)
  {
    for (i=0;i<p;i++)
      if (scrambletype==QL_SCRAMBLE_POWER)
	scrambleTable.push_back(scrambledig(i,p));
      else
	scrambleTable.push_back(i);
    for (i=0;i<pp[1];i++)
    {
      acc=0;
      dec=i;
      for (j=0;j<pp[0];j++)
      {
	acc=p*acc+scrambleTable[dec%p];
	dec/=p;
      }
      reverseScrambleTable[inx].push_back(acc);
    }
  }
}

int quadlods::reverseScramble(int limb,int p,int scrambletype)
{
  int inx=(scrambletype<<16)+p;
  fillReverseScrambleTable(p,scrambletype);
  if (reverseScrambleTable[inx].size())
    return reverseScrambleTable[inx][limb];
  else
    return limb;
}

bool quadlods::incHacc(std::vector<unsigned short> &hacc,int pp,int increment,int pos,bool sign)
/* Increments a Halton accumulator, whose prime power is pp, by increment,
 * which should be less in absolute value than pp, starting at the posth limb,
 * and returns the sign of the result (true is negative).
 */
{
  int i,limb,carry;
  for (i=pos;increment;i++)
  {
    while (hacc.size()<=i)
      hacc.push_back(sign?(pp-1):0);
    limb=hacc[i]+increment;
    increment=0;
    while (limb>=pp)
    {
      increment++;
      limb-=pp;
    }
    while (limb<0)
    {
      increment--;
      limb+=pp;
    }
    if (i==hacc.size()-1)
    {
      if (increment==1 && sign)
      {
	increment=0;
	sign=false;
      }
      if (increment==-1 && !sign)
      {
	increment=0;
	sign=true;
      }
    }
    hacc[i]=limb;
  }
  return sign;
}

bool quadlods::incHacc(std::vector<unsigned short> &hacc,int pp,mpz_class increment,bool sign)
{
  int ex=0,limb;
  mpz_class ppex=1;
  while (fabs(increment.get_d()/ppex.get_d())>2*pp)
  {
    ex++;
    ppex*=pp;
  }
  while (ex>=0)
  {
    limb=lrint(increment.get_d()/ppex.get_d());
    sign=incHacc(hacc,pp,limb,ex,sign);
    increment-=limb*ppex;
    ex--;
    ppex/=pp;
  }
  return sign;
}

mpz_class quadlods::haccValue(vector<unsigned short> &hacc,int pp,bool sign)
{
  mpz_class ret=-sign;
  int i;
  for (i=hacc.size()-1;i>=0;i--)
    ret=ret*pp+hacc[i];
  return ret;
}

mpz_class quadlods::thuemorse(int n)
{
  while (morse<=n)
  {
    thue+=(mpz_class)(((thue&((mpz_class)1<<(morse>>5)))>0)?(unsigned)0x96696996:0x69969669)<<morse;
    morse+=32;
  }
  return thue&(((mpz_class)1<<n)-1);
}

mpz_class quadlods::minusthird(int n)
{
  while (b2adic<=n)
  {
    third+=(mpz_class)0x55555555<<b2adic;
    b2adic+=32;
  }
  return third&(((mpz_class)1<<n)-1);
}

mpz_class quadlods::graydecode(mpz_class n)
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

mpz_class quadlods::scramble(mpz_class acc,mpz_class denom,int scrambletype)
{
  int i;
  mpz_class bitdiff,hibits,lobits,ret;
  bitdiff=denom&~acc;
  i=mpz_sizeinbase(bitdiff.get_mpz_t(),2)-1;
  switch (scrambletype)
  {
    case QL_SCRAMBLE_THIRD:
      ret=acc^minusthird(i);
      break;
    case QL_SCRAMBLE_THUEMORSE:
      ret=acc^thuemorse(i);
      break;
    case QL_SCRAMBLE_GRAY:
      lobits=acc&(((mpz_class)1<<i)-1);
      hibits=acc-lobits;
      ret=hibits+graydecode(lobits);
      break;
    default:
      ret=acc;
  }
  return ret;
}

short quadlods::readshort(std::istream &file)
{
  char buf[2];
  file.read(buf,2);
  return *(short *)buf;
}

void quadlods::initprimes()
{
  int i,j,n;
  int primeCheck=0;
  bool prime;
  PrimeContinuedFraction pcf;
  ifstream primeFile(string(SHARE_DIR)+"/primes.dat",ios::binary);
  primes.clear();
  for (i=2;i<65535;i++)
  {
    for (j=0,prime=true;j<primes.size() && prime && (unsigned int)primes[j]*primes[j]<=i;j++)
      if (i%primes[j]==0)
	prime=false;
    if (prime)
    {
      primeCheck+=i;
      primes.push_back(i);
    }
  }
  primesCfSorted.clear();
  for (i=0;i<QL_MAX_DIMS;i++)
  {
    pcf.prime=(unsigned short)readshort(primeFile);
    primeCheck-=pcf.prime;
    n=readshort(primeFile);
    pcf.cf.period=readshort(primeFile);
    pcf.cf.terms.clear();
    for (j=0;j<n;j++)
      pcf.cf.terms.push_back(readshort(primeFile));
    primesCfSorted.push_back(pcf);
  }
  if (primeCheck)
  {
    cerr<<"Prime file is corrupt or failed to load"<<endl;
    primesCfSorted.clear();
  }
}

void quadlods::compquad(int p,double resolution,mpz_class &nmid,mpz_class &dmid)
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

void quadlods::compquad(ContinuedFraction cf,double resolution,mpz_class &nmid,mpz_class &dmid)
{
  mpz_class nhi,dhi,nlo,dlo;
  int i=0,j=0;
  bool comp=false;
  for (nhi=dlo=1,nlo=dhi=dmid=0;dmid<resolution;)
  {
    dmid=dhi+dlo;
    nmid=nhi+nlo;
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
    if (++j>=cf.terms[i])
    {
      comp=!comp;
      j=0;
      if (++i>=cf.terms.size())
	i-=cf.period;
    }
  }
  nmid%=dmid;
}

int quadlods::nthprime(int n)
{
  if (primes.size()==0)
    initprimes();
  if (n<0 || n>=primes.size())
    return 0;
  else if (primesCfSorted.size())
    return primesCfSorted[n].prime;
  else
    return primes[n];
}

double quadlods::nthquad(int n)
{
  mpz_class nmid,dmid;
  if (n<0 || n>=primes.size())
    return NAN;
  else
  {
    if (primesCfSorted.size())
      compquad(primesCfSorted[n].cf,27/DBL_EPSILON,nmid,dmid);
    else
    {
      compquad(nthprime(n),27/DBL_EPSILON,nmid,dmid);
    }
    return mpq_class(nmid,dmid).get_d();
  }
}

double ContinuedFraction::averageTerm() const
{
  int i,sixteens=0;
  double product=1;
  if (period<=0 || period>terms.size())
    return NAN;
  else
  {
    for (i=0;i<period;i++)
    {
      product*=terms[terms.size()-i-1];
      if (product>16)
      {
	product/=16;
	sixteens++;
      }
    }
    return pow(product,1./period)*pow(16,(double)sixteens/period);
  }
}

int ContinuedFraction::maximumTerm() const
{
  int i,max=0;
  for (i=0;i<period && i<terms.size();i++)
    if (max<terms[terms.size()-i-1])
      max=terms[terms.size()-i-1];
  return max;
}

bool operator<(const PrimeContinuedFraction a,const PrimeContinuedFraction b)
{
  double aavg=a.cf.averageTerm(),bavg=b.cf.averageTerm();
  int amax=a.cf.maximumTerm(),bmax=b.cf.maximumTerm();
  if (amax!=bmax)
    return amax<bmax;
  else if (fabs(aavg-bavg)>1e-10)
    return aavg<bavg;
  else
    return a.prime<b.prime;
}

Quadlods::Quadlods()
{
  mode=QL_MODE_RICHTMYER;
  scrambletype=QL_SCRAMBLE_NONE;
  sign=false;
}

void Quadlods::init(int dimensions,double resolution,int j)
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
  int i,p,newmode;
  mpz_class nmid,dmid;
  if (dimensions>QL_MAX_DIMS)
    dimensions=QL_MAX_DIMS;
  if (dimensions<-QL_MAX_DIMS)
    dimensions=-QL_MAX_DIMS;
  newmode=resolution?QL_MODE_RICHTMYER:QL_MODE_HALTON;
  if (mode!=newmode)
    primeinx.clear();
  if (newmode==QL_MODE_HALTON)
  {
    num.clear();
    denom.clear();
    acc.clear();
    for (i=hacc.size();i<dimensions;i++)
    {
      primeinx.push_back(i);
      p=nthprime(i);
      hacc.resize(primeinx.size());
      if (i)
	incHacc(hacc[i],primePower(p)[1],
		haccValue(hacc[i-1],primePower(nthprime(primeinx[i-1]))[1],sign),false);
    }
    for (i=hacc.size();i>dimensions;i--)
    {
      primeinx.push_back(QL_MAX_DIMS+i-1);
      p=nthprime(QL_MAX_DIMS+i-1);
      hacc.resize(primeinx.size());
      if (i)
	incHacc(hacc[i],primePower(p)[1],
		haccValue(hacc[i-1],primePower(nthprime(primeinx[i-1]))[1],sign),false);
    }
  }
  else
  {
    hacc.clear();
    for (i=denom.size();i<dimensions;i++)
    {
      primeinx.push_back(i);
      if (primesCfSorted.size())
	compquad(primesCfSorted[i].cf,resolution,nmid,dmid);
      else
      {
	p=nthprime(i);
	compquad(p,resolution,nmid,dmid);
      }
      denom.push_back(dmid);
      num.push_back(nmid);
    }
    for (i=-denom.size();i>dimensions;i--)
    {
      primeinx.push_back(QL_MAX_DIMS+i-1);
      if (primesCfSorted.size())
	compquad(primesCfSorted[QL_MAX_DIMS+i-1].cf,resolution,nmid,dmid);
      else
      {
	p=nthprime(QL_MAX_DIMS+i-1);
	compquad(p,resolution,nmid,dmid);
      }
      denom.push_back(dmid);
      num.push_back(nmid);
    }
  }
  if (j>=0)
    scrambletype=j;
  if (scrambletype<0 || scrambletype>QL_SCRAMBLE_GRAY)
    scrambletype=QL_SCRAMBLE_GRAY;
  mode=newmode;
  dimensions=abs(dimensions);
  if (mode==QL_MODE_RICHTMYER)
  {
    num.resize(dimensions);
    denom.resize(dimensions);
    acc.resize(dimensions);
  }
}

void Quadlods::init(vector<int> dprimes,double resolution,int j)
/* Used for testing, setting up a generator with a particularly
 * bad set of primes.
 */
{
  int i,k,p,newmode;
  mpz_class nmid,dmid;
  // TODO check that dprimes are actually distinct
  newmode=resolution?QL_MODE_RICHTMYER:QL_MODE_HALTON;
  if (mode!=newmode)
    primeinx.clear();
  for (i=primeinx.size();i<dprimes.size();i++)
    for (k=0;k<QL_MAX_DIMS;k++)
      if (dprimes[i]==nthprime(k))
      {
        primeinx.push_back(k);
        k=8191;
      }
  for (i=newmode?hacc.size():denom.size();i<primeinx.size();i++)
  {
    if (newmode==QL_MODE_RICHTMYER)
    {
      if (primesCfSorted.size())
	compquad(primesCfSorted[primeinx[i]].cf,resolution,nmid,dmid);
      else
      {
	p=nthprime(primeinx[i]);
	compquad(p,resolution,nmid,dmid);
      }
      denom.push_back(dmid);
      num.push_back(nmid);
    }
    if (newmode==QL_MODE_HALTON)
    {
      if (primesCfSorted.size())
	p=primesCfSorted[primeinx[i]].prime;
      else
	p=nthprime(primeinx[i]);
      hacc.push_back(vector<unsigned short>());
    }
  }
  if (j>=0)
    scrambletype=j;
  if (scrambletype<0 || scrambletype>QL_SCRAMBLE_GRAY)
    scrambletype=QL_SCRAMBLE_GRAY;
  mode=newmode;
  if (mode==QL_MODE_RICHTMYER)
  {
    num.resize(primeinx.size());
    denom.resize(primeinx.size());
    acc.resize(primeinx.size());
  }
}

mpz_class Quadlods::gethacc(int n)
// All Halton accumulators should have the same value, so n shouldn't matter.
{
  mpz_class ret=sign?-1:0;
  int i,limbbase;
  if (hacc.size())
  {
    n%=hacc.size();
    if (n<0)
      n+=hacc.size();
    limbbase=primePower(nthprime(primeinx[n]))[1];
    for (i=hacc[n].size()-1;i>=0;i--)
      ret=ret*limbbase*hacc[n][i];
  }
  return ret;
}

vector<mpq_class> Quadlods::readout()
{
  int i;
  vector<mpq_class> ret;
  for (i=0;i<num.size();i++)
  {
    ret.push_back(mpq_class((scramble(acc[i],denom[i],scrambletype)<<1)|1,denom[i]<<1));
    ret[i].canonicalize();
  }
  return ret;
}

vector<double> Quadlods::dreadout()
{
  int i,pp;
  vector<double> ret;
  for (i=0;mode==QL_MODE_RICHTMYER && i<num.size();i++)
  {
    ret.push_back(mpq_class((scramble(acc[i],denom[i],scrambletype)<<1)|1,denom[i]<<1).get_d());
  }
  for (i=0;mode==QL_MODE_HALTON && i<hacc.size();i++)
  {
    pp=primePower(nthprime(primeinx[i]))[1];
    ret.push_back(mpq_class(haccValue(hacc[i],pp,sign),mpz_class(1)).get_d());
  }
  return ret;
}

void Quadlods::setmiddle()
/* Set the point to the middle of the square/cube/etc. This is used in testing
 * to find the successively smaller distances between q[n] and q[n+h].
 */
{
  int i;
  for (i=0;i<num.size();i++)
    acc[i]=denom[i]>>1;
}

void Quadlods::advance(mpz_class n)
{
  int i,pp;
  bool newsign=sign;
  for (i=0;i<num.size();i++)
    if (n<0)
      acc[i]=(acc[i]-n*(denom[i]-num[i]))%denom[i];
    else
      acc[i]=(acc[i]+n*num[i])%denom[i];
  for (i=0;i<hacc.size();i++)
  {
    pp=primePower(nthprime(primeinx[i]))[1];
    newsign=incHacc(hacc[i],pp,n,sign);
  }
  sign=newsign;
}

unsigned int Quadlods::seedsize()
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

void Quadlods::seed(char *s,unsigned int n)
{
  unsigned i,sz;
  sz=denom.size();
  for (i=0;i<n;i++)
    acc[i%sz]=((acc[i%sz]<<8)+(s[i]&0xff))%denom[i%sz];
}

vector<mpq_class> Quadlods::gen()
{
  advance(1);
  return readout();
}

vector<double> Quadlods::dgen()
{
  advance(1);
  return dreadout();
}

Quadlods select(Quadlods& b,vector<int> dimensions)
{
  Quadlods ret;
  int i,j;
  ret.scrambletype=b.scrambletype;
  for (i=0;i<dimensions.size();i++)
    if (dimensions[i]>=0 && dimensions[i]<b.size())
    {
      for (j=0;j<ret.size() && ret.primeinx[j]!=b.primeinx[dimensions[i]];j++);
      if (j==ret.size())
      {
        ret.num.     push_back(b.num     [dimensions[i]]);
        ret.denom.   push_back(b.denom   [dimensions[i]]);
        ret.acc.     push_back(b.acc     [dimensions[i]]);
        ret.primeinx.push_back(b.primeinx[dimensions[i]]);
      }
    }
  return ret;
}
