/* Copyright 2014,2016,2017 Pierre Abbat.
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
#include <gmpxx.h>

#define QL_JUMBLE_NONE 0
#define QL_JUMBLE_THIRD 1
#define QL_JUMBLE_THUEMORSE 2
#define QL_JUMBLE_GRAY 3
/* The jumbletype controls how to jumble the bits of the accumulator when
 * reading the generator. If acc is 0xc0de and denom is 0x10000, it returns:
 * QL_JUMBLE_NONE:      1100000011011110
 * QL_JUMBLE_THIRD:     1001010110001011 (xor with 5555, 1/3=0.5555... hex)
 * QL_JUMBLE_THUEMORSE: 0101011010110111 (xor with the Thue-Morse word)
 * QL_JUMBLE_GRAY:      1000000010010100 (inverse Gray code).
 * Of these, Gray code is the best, as it leaves no tendency to slope
 * one way or the other. It is therefore the default. However, Gray code
 * takes about 5/3 times as much time as the others.
 */
#define QL_MAX_DIMS 6542

class quadirr
/* Represents the quadratic irrational a/b+c*sqrt(p)/d.
 * Used for computing the continued fraction representation of sqrt(p) or
 * (sqrt(p)+1)/2 and, thus, order the numbers by their discrepancy constant.
 */
{
private:
  int a,b,c,d,p;
public:
  quadirr();
  quadirr(int A,int B,int C,int D,int P);
  double realval();
  bool operator=(const quadirr &r) const;
  quadirr& operator-=(int n);
  void recip();
};

int nthprime(int n);
double nthquad(int n);

class quadlods
{
protected:
  std::vector<mpz_class> num,denom,acc;
  std::vector<short> primeinx;
  int jumbletype;
public:
  void init(int dimensions,double resolution,int j=QL_JUMBLE_GRAY);
  // If dimensions>6542, it is silently truncated to 6542.
  void init(std::vector<int> dprimes,double resolution,int j=QL_JUMBLE_GRAY);
  int size()
  {
    return acc.size();
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
  int getprimeinx(int n)
  {
    if (n<0 || n>primeinx.size())
      return -1;
    else
      return primeinx[n];
  }
  int getprime(int n)
  {
    return nthprime(getprimeinx(n));
  }
  std::vector<mpq_class> readout();
  void setmiddle();
  void setjumble(int j)
  {
    jumbletype=j;
  }
  void advance(mpz_class n);
  unsigned int seedsize();
  void seed(char *s,unsigned int n);
  std::vector<mpq_class> gen();
  std::vector<double> dreadout();
  std::vector<double> dgen();
  friend quadlods select(quadlods& b,std::vector<int> dimensions);
};
#endif
