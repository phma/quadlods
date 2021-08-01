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
#include <fftw3.h>
#include "fourier.h"
#include "histogram.h"
#include "discrepancy.h"
#include "ldecimal.h"

#define BUCKETS 2971
/* for graphing. On A4 landscape, this produces more than 10
 * vertical lines per millimeter.
 */

using namespace std;
using namespace quadlods;

map<int,fftw_plan> plans;
map<int,double *> inmem,outmem;

fftw_plan makePlan(int n)
// fftw_plan is a pointer to plan_s, so it can be checked as bool
{
  if (!plans[n])
  {
    inmem[n]=fftw_alloc_real(n);
    outmem[n]=fftw_alloc_real(n);
    plans[n]=fftw_plan_r2r_1d(n,inmem[n],outmem[n],FFTW_R2HC,FFTW_ESTIMATE);
  }
  return plans[n];
}

void destroyPlans()
{
  map<int,fftw_plan>::iterator i;
  map<int,double *>::iterator j;
  for (i=plans.begin();i!=plans.end();i++)
    fftw_destroy_plan(i->second);
  for (j=inmem.begin();j!=inmem.end();j++)
    fftw_free(j->second);
  for (j=outmem.begin();j!=outmem.end();j++)
    fftw_free(j->second);
  plans.clear();
  inmem.clear();
  outmem.clear();
}

vector<double> fft(vector<double> input)
/* The output is calibrated so that the frequency-domain terms are independent
 * of the size of the input.
 */
{
  vector<double> output(input);
  int i,sz=input.size();
  fftw_plan plan=makePlan(sz);
  memcpy(inmem[sz],&input[0],sz*sizeof(double));
  fftw_execute(plan);
  for (i=0;i<sz;i++)
    output[i]=outmem[sz][i]/sz;
  return output;
}

double window(int i,int iters)
{
  return 1-cos((i+0.5)/iters*M_PI*2);
}

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
