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
#include <iostream>
#include <cfloat>
#include <cmath>
#include "contfrac.h"

using namespace std;

quadirr nthquadQi(int n)
{
  int p=nthprime(n);
  if ((p-1)&3)
    return quadirr(0,1,1,1,p);
  else
    return quadirr(1,2,1,2,p);
}

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

bool quadirr::operator==(const quadirr &r) const
{
  if (c==0 && r.c==0)
    return a*r.b==b*r.a;
  else
    return a*r.b==b*r.a && p==r.p && c*r.d==d*r.c;
}

quadirr& quadirr::operator-=(int n)
{
  a-=b*n;
  return *this;
}

unsigned gcd(unsigned a,unsigned b)
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

/* To compute the reciprocal of a/b+c*sqrt(p)/d:
 * 1/(a/b+c*sqrt(p)/d)
 * (a/b-c*sqrt(p)/d)/(a²/b²-c²*p/d²)
 * 
 * 1/2+1*sqrt(5)/2-1=-1/2+1*sqrt(5)/2
 * (-1/2-1*sqrt(5)/2)/(1/4-1*5/4)=1/2+1*sqrt(5)/2
 * 
 * (1,2,1,2,65029)-128=(-255,2,1,2,65029)
 * 1/(-255,2,1,2,65029)=(-255,2,-1,2,65029)/(65025/4-65029/4)=(255,2,1,2,65029)
 * (255,2,1,2,65029)-255=(-255,2,1,2,65029)
 * (128;255,255,...)
 * 
 * (0,1,1,1,65027)-255=(-255,1,1,1,65027)
 * 1/(-255,1,1,1,65027)=(-255,1,-1,1,65027)/(65025-65027)=(255,2,1,2,65027)
 * (255,2,1,2,65027)-255=(-255,2,1,2,65027)
 * 1/(-255,2,1,2,65027)=(-255,2,-1,2,65027)/(65025/4-65027/4)=(255,1,1,1,65027)
 * (255,1,1,1,65027)-510=(-255,1,1,1,65027)
 * (255;255,510,255,510,...)
 */

void quadirr::recip()
{
  mpq_class denom,e(a,b),f(c,d);
  int g,h,gc;
  e.canonicalize();
  f.canonicalize();
  e*=e;
  f=f*f*p;
  denom=e-f;
  g=denom.get_num().get_si();
  h=denom.get_den().get_si();
  a*=h;
  b*=g;
  c*=-h;
  d*=g;
  gc=gcd(abs(a),abs(b));
  a/=gc;
  b/=gc;
  gc=gcd(abs(c),abs(d));
  c/=gc;
  d/=gc;
  if (b<0)
  {
    a=-a;
    b=-b;
  }
  if (d<0)
  {
    c=-c;
    d=-d;
  }
}

ContinuedFraction contFrac(quadirr q)
{
  ContinuedFraction ret;
  vector<quadirr> partials;
  int i;
  bool done=false;
  while (!done)
  {
    if (partials.size())
    {
      partials.push_back(partials.back());
      partials.back()-=ret.terms.back();
      partials.back().recip();
    }
    else
      partials.push_back(q);
    ret.terms.push_back(floor(partials.back().realval()));
    for (i=1;!done && i<partials.size();i++)
      if (partials[partials.size()-1]==partials[partials.size()-1-i])
      {
	done=true;
	ret.period=i;
      }
  }
  return ret;
}
