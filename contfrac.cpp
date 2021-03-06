/******************************************************/
/*                                                    */
/* contfrac.cpp - continued fraction expansions       */
/*                                                    */
/******************************************************/
/* Copyright 2018,2019 Pierre Abbat.
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
using namespace quadlods;

quadirr nthquadQi(int n)
/* Returns the quadratic irrational of the nth prime.
 * If the prime is congruent to 1 mod 4, it returns (sqrt(p)+1)/2, whose maximum
 * continued fraction term is half or less that of sqrt(p).
 * If the prime is congruent to 3 mod 4, it returns sqrt(p). The continued
 * fraction terms of sqrt(p) and (sqrt(p)+1)/2 are the same, except for the
 * first term and a rotation.
 * If the prime is 2, sqrt(p) has the smaller maximum term.
 */
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
  if (p<0)
    throw IMAGINARY;
  return a/(double)b+c*sqrt(p)/d;
}

string quadirr::stringval()
{
  string ratpart,quadpart;
  if (b==1)
    ratpart=to_string(a);
  else
    ratpart=to_string(a)+'/'+to_string(b);
  quadpart="√"+to_string(p);
  if (d!=1)
    quadpart+='/'+to_string(d);
  if (abs(c)!=1)
    quadpart=to_string(abs(c))+quadpart;
  quadpart=((c<0)?'-':'+')+quadpart;
  return ratpart+quadpart;
}

bool quadirr::operator==(const quadirr &r) const
{
  if (c==0 && r.c==0)
    return a*r.b==b*r.a;
  else
    return a*r.b==b*r.a && p==r.p && c*r.d==d*r.c;
}

bool quadirr::is0()
{
  int i;
  for (i=2;i*i<=p;i++)
    if (p%(i*i)==0)
    {
      c*=i;
      p/=i*i;
      i--;
    }
  if (p==1)
    return b!=0 && d!=0 && a*d+b*c==0;
  else
    return (a==0 && b!=0) && ((c==0 || p==0) && d!=0);
}

quadirr& quadirr::operator-=(int n)
{
  a-=b*n;
  return *this;
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
  long long A,B,C,D;
  if (b==0 || d==0)
    throw ZERODIV;
  e.canonicalize();
  f.canonicalize();
  e*=e;
  f=f*f*p;
  denom=e-f;
  g=denom.get_num().get_si();
  h=denom.get_den().get_si();
  if (g!=denom.get_num() || h!=denom.get_den())
    throw OVERFLOW;
  A=(long long)a*h;
  B=(long long)b*g;
  C=(long long)c*-h;
  D=(long long)d*g;
  a=A;
  b=B;
  c=C;
  d=D;
  if (a!=A || b!=B || c!=C || d!=D)
    throw OVERFLOW;
  gc=gcd(abs(a),abs(b));
  if (gc==0)
    throw ZERODIV;
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
      if (partials.back().is0())
      {
	done=true;
	ret.period=0;
      }
      else
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
  ret.terms.pop_back();
  return ret;
}

QuadMax equivClass(quadirr q)
/* There is an equivalence relation on numbers defined by ξ≡η iff ξ=(aη+b)/(cη+d)
 * for integers a,b,c,d where ad-bc=±1. If the numbers have continued fraction
 * expansions, this is the same as saying that the expansions have identical
 * tails. All rational numbers are equivalent, and the representative of this
 * equivalence class is 0. For periodic continued fractions (real quadratic
 * irrationals), the representative can be chosen as the least or greatest
 * number in the periodic part; this function returns the least. For numbers
 * in general, defining a representative requires the axiom of choice.
 */
{
  QuadMax ret;
  vector<quadirr> partials;
  int i,period;
  bool done=false;
  double realFrac;
  while (!done)
  {
    if (partials.size())
    {
      partials.push_back(partials.back());
      partials.back()-=floor(partials.back().realval());
      realFrac=partials.back().realval();
      partials.back().recip();
      if (fabs(partials.back().realval()*realFrac-1)>0.1)
	cerr<<"Reciprocal failed\n"; // This should never happen; it should throw.
    }
    else
      partials.push_back(q);
    for (i=0;!done && i<partials.size()-1;i=2*i+1)
      if (partials[partials.size()-1]==partials[i])
      {
	done=true;
	period=partials.size()-1-i;
      }
  }
  ret.qi=partials.back();
  ret.max=floor(partials.back().realval());
  for (i=1;i<period;i++)
  {
    if (partials[partials.size()-1-i].realval()<ret.qi.realval())
      ret.qi=partials[partials.size()-1-i];
    if (partials[partials.size()-1-i].realval()>ret.max)
      ret.max=floor(partials[partials.size()-1-i].realval());
  }
  return ret;
}
