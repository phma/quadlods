/******************************************************/
/*                                                    */
/* contfrac.h - continued fraction expansions         */
/*                                                    */
/******************************************************/
/* Copyright 2018 Pierre Abbat.
 * This file is part of the Quadlods program.
 * 
 * The Quadlods program is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
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
#include "quadlods.h"

#define OVERFLOW 1
#define ZERODIV 2
#define IMAGINARY 3

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
  int geta()
  {
    return a;
  };
  int getb()
  {
    return b;
  };
  int getc()
  {
    return c;
  };
  int getd()
  {
    return d;
  };
  int getp()
  {
    return p;
  };
  double realval();
  std::string stringval();
  bool operator==(const quadirr &r) const;
  bool is0();
  quadirr& operator-=(int n);
  void recip();
};

struct QuadMax
{
  quadirr qi;
  int max;
};

quadirr nthquadQi(int n);
ContinuedFraction contFrac(quadirr q);
QuadMax equivClass(quadirr q);
