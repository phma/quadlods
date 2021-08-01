/******************************************************/
/*                                                    */
/* fourier.cpp - plot Fourier transforms of sequence  */
/*                                                    */
/******************************************************/
/* Copyright 2021 Pierre Abbat.
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
#include "discrepancy.h"
#include "ldecimal.h"

using namespace std;
using namespace quadlods;

void fouriertest(Quadlods &quad,int iters,PostScript &ps)
/* Plot the Fourier transforms of the sequence.
 */
{
  int i,j,inx,allinx;
  time_t now,then;
  Quadlods sel1;
  vector<int> pinx1;
  vector<double> point,fpoint;
  vector<vector<double> > points;
  double ang,r,disc;
  setFlowerDisc(true);
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
    ps.write(0.8*sqrt(iters),0.8*sqrt(iters),to_string(sel1.getprime(0)));
    points.clear();
    for (i=0;i<iters;i++)
    {
      point=sel1.dgen();
      r=sqrt(i+0.5);
      ang=2*M_PI*point[0];
      ps.dot(r*cos(ang),r*sin(ang));
      fpoint.clear();
      fpoint.push_back(r*cos(ang)/sqrt(iters));
      fpoint.push_back(r*sin(ang)/sqrt(iters));
      points.push_back(fpoint);
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
