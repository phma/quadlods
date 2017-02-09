/* Copyright 2017 Pierre Abbat.
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
#include <vector>
#include <ctime>
#include <cmath>
#include "circletest.h"
#include "ps.h"

using namespace std;

double halfstep[]={1,1.05946309435929526455,1.12246204830937298142};
double minorthird[]={1,1.18920711500272106671,1.41421356237309504878,1.68179283050742908604};
papersize a4land={297,210};

double hstep(int i)
{
  int octave,note;
  octave=i/12;
  if (i<octave*12)
    --octave;
  note=i-12*octave;
  return ldexp(minorthird[note/3]*halfstep[note%3],octave);
}

void circletest(quadlods &quad)
/* Find all pairs of primes (dimensions) which have high discrepancy by
 * counting the pairs of points that fall inside a quadrant of a circle.
 * The number of points in the quadrant should be (π/4)i+O(ln(i)²/i).
 */
{
  int i,j,k,inx,iters=1048576,logi=-1;
  char buf[24];
  time_t now,then;
  bool recordthis;
  vector<int> incircle,ri;
  vector<double> point;
  vector<vector<double> > rrelError;
  PostScript ps;
  double relativeError,maxError,scale;
  incircle.resize(quad.size()*(quad.size()-1)/2);
  rrelError.resize(quad.size()*(quad.size()-1)/2);
  for (i=0;i<iters;i++)
  {
    recordthis=false;
    while (hstep(logi)<i+1)
    {
      logi++;
      recordthis=true;
    }
    point=quad.dgen();
    for (j=0;j<quad.size();j++)
      for (k=0;k<j;k++)
      {
        inx=j*(j-1)/2+k;
        if (point[j]*point[j]+point[k]*point[k]<1)
          incircle[inx]++;
        if (point[j]*point[j]+(1-point[k])*(1-point[k])<1)
          incircle[inx]--;
        if (recordthis)
        {
          relativeError=(incircle[inx])/log(i+1)/log(i+1);
          if (i==0)
            relativeError=0;
          rrelError[inx].push_back(relativeError);
        }
      }
    if (recordthis)
      ri.push_back(i+1);
    now=time(nullptr);
    if (now!=then)
    {
      cout<<rint((double)i/iters*100)<<"% \r";
      cout.flush();
      then=now;
    }
  }
  for (j=0;j<quad.size();j++)
    for (k=0;k<j;k++)
    {
      inx=j*(j-1)/2+k;
      relativeError=(incircle[inx])/log(i)/log(i);
      cout<<j<<' '<<k<<' '<<incircle[inx]<<' '<<relativeError<<endl;
    }
  ps.open("circletest.ps");
  ps.setpaper(a4land,0);
  ps.prolog();
  for (i=0;i<rrelError.size();i++)
  {
    for (maxError=j=0;j<ri.size();j++)
      if (fabs(rrelError[i][j])>maxError)
        maxError=fabs(rrelError[i][j]);
    ps.startpage();
    ps.setscale(0,-1,3,1);
    scale=maxError;
    ps.startline();
    ps.lineto(0,-1);
    ps.lineto(3,-1);
    ps.lineto(3,1);
    ps.lineto(0,1);
    ps.endline(true);
    sprintf(buf,"%g",scale);
    ps.write(3,1,buf);
    ps.startline();
    for (j=0;j<ri.size();j++)
      ps.lineto(log(ri[j])/log(iters)*3,rrelError[i][j]/maxError);
    ps.endline();
    ps.endpage();
  }
}
