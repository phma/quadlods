/******************************************************/
/*                                                    */
/* flowertest.cpp - draw flower diagrams of sequence  */
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
#include <iomanip>
#include <ctime>
#include <cmath>
#include "flowertest.h"

using namespace std;

void flowertest(Quadlods &quad,int iters,PostScript &ps)
/* Draw a flower diagram of the sequence. The flower diagram of an unjumbled
 * sequence with step φ (from prime 5) is the pattern of flowers in an
 * asteraceous flower head.
 */
{
  int i,j,k,inx,allinx;
  char buf[24];
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
