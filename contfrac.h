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
  bool operator==(const quadirr &r) const;
  quadirr& operator-=(int n);
  void recip();
};

quadirr nthquadQi(int n);
ContinuedFraction contFrac(quadirr q);
