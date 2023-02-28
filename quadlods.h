/******************************************************/
/*                                                    */
/* quadlods.h - quadratic low-discrepancy sequence    */
/*                                                    */
/******************************************************/
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
#ifndef QUADLODS_H
#define QUADLODS_H
#include <vector>
#include <map>
#include <gmpxx.h>

#define QL_MODE_RICHTMYER 0
#define QL_MODE_HALTON 1
#define QL_SCRAMBLE_NONE 0
#define QL_SCRAMBLE_THIRD 1
#define QL_SCRAMBLE_THUEMORSE 2
#define QL_SCRAMBLE_GRAY 3
#define QL_SCRAMBLE_POWER 4
#define QL_SCRAMBLE_FAURE 5
#define QL_SCRAMBLE_TIPWITCH 6
#define QL_SCRAMBLE_DEFAULT 255
/* The scrambletype controls how to scramble the bits of the accumulator when
 * reading the generator. If acc is 0xc0de and denom is 0x10000, it returns:
 * QL_SCRAMBLE_NONE:      1100000011011110
 * QL_SCRAMBLE_THIRD:     1001010110001011 (xor with 5555, 1/3=0.5555... hex)
 * QL_SCRAMBLE_THUEMORSE: 0101011010110111 (xor with the Thue-Morse word)
 * QL_SCRAMBLE_GRAY:      1000000010010100 (inverse Gray code).
 * Of these, Gray code is the best, as it leaves no tendency to slope
 * one way or the other. It is therefore the default. However, Gray code
 * takes about 5/3 times as much time as the others.
 *
 * QL_SCRAMBLE_POWER is for Halton. It scrambles each digit by raising it to
 * a power relatively prime to the totient, leaving -1, 0, and 1 untouched.
 * QL_SCRAMBLE_FAURE is the Faure permutation scrambling for Halton. The
 * permutation is defined as follows: if n is even, take the permutation of n/2
 * twice, double the numbers, and add 1 to the second half; if n is odd, stick
 * a point in the middle, spreading the four quadrants apart.
 * QL_SCRAMBLE_TIPWITCH uses a sequence of permutations computed by combining
 * base reversal, flipping upside down, and interdigitating permutations.
 * It is named for the tipitiwitchet, whose leaves have spikes on the edges,
 * which interdigitate when it closes. I tried base-reversing "tippitiwitchet"
 * and "aldrovanda" (a related plant), but the results were unpronounceable.
 * Tipwitch is the default scrambling for Halton.
 */
#define QL_MAX_DIMS 6542

namespace quadlods
{
  extern std::vector<unsigned short> primes;
  extern std::map<unsigned,std::vector<unsigned short> > reverseScrambleTable;
  unsigned gcd(unsigned a,unsigned b);
  int nthprime(int n);
  double nthquad(int n,bool mod1=false);
  unsigned relprime(unsigned n);
  unsigned scrambledig(unsigned dig,unsigned p);
  unsigned faureperm(unsigned dig,unsigned p);
  void fillReverseScrambleTable(int p,int scrambletype);
  bool incHacc(std::vector<unsigned short> &hacc,int pp,int increment,int pos,bool sign);
  bool incHacc(std::vector<unsigned short> &hacc,int pp,mpz_class increment,bool sign);
  mpz_class haccValue(std::vector<unsigned short> &hacc,int pp,bool sign);
}

class ContinuedFraction
// Represents a periodic continued fraction, i.e. a quadratic number.
{
public:
  std::vector<int> terms;
  int period; // 0 means it terminates
  double averageTerm() const;
  int maximumTerm() const;
};

class PrimeContinuedFraction
{
public:
  int prime;
  ContinuedFraction cf;
  friend bool operator<(const PrimeContinuedFraction &a,const PrimeContinuedFraction &b);
};

class Quadlods
{
protected:
  std::vector<mpz_class> num,denom,acc;
  std::vector<std::vector<unsigned short> > hacc;
  std::vector<short> primeinx;
  int scrambletype;
  int mode;
  bool sign;
public:
  Quadlods();
  void init(int dimensions,double resolution,int j=QL_SCRAMBLE_DEFAULT);
  /* If dimensions>6542, it is silently truncated to 6542.
   * If dimensions<0, primes are taken from the end of the list.
   * Since Halton does not use the resolution, setting the resolution to 0
   * sets the mode to Halton.
   */
  void init(std::vector<int> dprimes,double resolution,int j=QL_SCRAMBLE_DEFAULT);
  int size()
  {
    return mode?hacc.size():acc.size();
  }
  int getMode()
  {
    return mode;
  }
  mpz_class getnum(int n)
  {
    return num[n];
  }
  mpz_class getdenom(int n)
  {
    return denom[n];
  }
  mpz_class getacc(int n)
  {
    return acc[n];
  }
  mpz_class gethacc(int n=0);
  int getprimeinx(int n)
  {
    if (n<0 || n>primeinx.size())
      return -1;
    else
      return primeinx[n];
  }
  int getprime(int n)
  {
    return quadlods::nthprime(getprimeinx(n));
  }
  std::vector<mpq_class> readout();
  std::vector<mpq_class> readoutUnscrambled();
  void setmiddle();
  void setscramble(int j);
  int getscramble()
  {
    return scrambletype;
  }
  void advance(mpz_class n);
  unsigned int seedsize();
  void seed(char *s,unsigned int n);
  std::vector<mpq_class> gen();
  std::vector<double> dreadout();
  std::vector<double> dreadoutUnscrambled();
  std::vector<double> dgen();
  friend Quadlods select(Quadlods& b,std::vector<int> dimensions);
};
#endif
