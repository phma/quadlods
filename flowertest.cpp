/******************************************************/
/*                                                    */
/* flowertest.cpp - draw flower diagrams of sequence  */
/*                                                    */
/******************************************************/
/* Copyright 2018-2020 Pierre Abbat.
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
#include <iomanip>
#include <ctime>
#include <cmath>
#include "flowertest.h"
#include "histogram.h"

using namespace std;
using namespace quadlods;

void flowertest(Quadlods &quad,int iters,PostScript &ps)
/* Draw a flower diagram of the sequence. The flower diagram of an unscrambled
 * sequence with step φ (from prime 5) is the pattern of flowers in an
 * asteraceous flower head.
 */
{
  int i,j,inx,allinx;
  time_t now,then;
  Quadlods sel1;
  vector<int> pinx1;
  vector<double> point;
  double ang,r;
  ps.setpaper(a4land,0);
  ps.prolog();
  allinx=iters*quad.size();
  for (j=0;j<quad.size();j++)
  {
    pinx1.clear();
    pinx1.push_back(j);
    sel1=select(quad,pinx1);
    ps.startpage();
    ps.setscale(-sqrt(iters),-sqrt(iters),sqrt(iters),sqrt(iters));
    for (i=0;i<iters;i++)
    {
      point=sel1.dgen();
      r=sqrt(i+0.5);
      ang=2*M_PI*point[0];
      ps.dot(r*cos(ang),r*sin(ang));
      if (((i-iters)&255)==255)
      {
	now=time(nullptr);
	inx=j*iters+i;
	if (now!=then)
	{
	  cout<<rint((double)inx/allinx*100)<<"% \r";
	  cout.flush();
	  then=now;
	}
      }
    }
    ps.endpage();
  }
}

void quadplot(PostScript &ps)
/* Draw all the quadratic irrationals mod 1 in the same manner as
 * the flower diagram.
 *
 * At a certain point, the curves of the graph bend. This happens when
 * primes greater than 65536 would be interspersed with primes less than
 * 65536, but the greater primes are not computed by the program and hence
 * not counted.
 *
 * If m is odd, the primes with maximum term m are between m² and (m+2)² and
 * are congruent to 1 mod 4. If m is even, the primes with maximum term m are
 * between (m/2)² and (m/2+1)² and are congruent to 2 or 3 mod 4.
 */
{
  int i,j;
  double point;
  double ang,r,x,y;
  double a,comb=1.625;
  histogram hist(-0.1,0.1);
  ps.setpaper(a4land,0);
  ps.prolog();
  for (j=0;j<1;j++)
  {
    ps.startpage();
    ps.setscale(-sqrt(QL_MAX_DIMS),-sqrt(QL_MAX_DIMS),sqrt(QL_MAX_DIMS),sqrt(QL_MAX_DIMS));
    ps.setcolor(1,0,1);
    a=-floor(sqrt(QL_MAX_DIMS)/comb);
    for (i=-a;i<=a;i++)
      ps.line2p(xy(-sqrt(QL_MAX_DIMS),i*comb),xy(sqrt(QL_MAX_DIMS),i*comb));
    ps.setcolor(0,0,0);
    for (i=0;i<QL_MAX_DIMS;i++)
    {
      point=nthquad(i);
      r=sqrt(i+0.5);
      ang=2*M_PI*point;
      x=r*cos(ang);
      y=r*sin(ang);
      if (i==4228) // max term increases from 256 to 258; 257 is missing
	ps.setcolor(1,0,1);
      ps.dot(x,y);
      if (x>2*abs(y))
      {
	hist<<y;
	//if (y<0 && y>-2.6*comb)
	  //cout<<i<<' '<<y/comb<<endl;
      }
    }
    ps.endpage();
  }
  //ps.startpage();
  //hist.plot(ps,HISTO_LINEAR);
  //ps.endpage();
}
