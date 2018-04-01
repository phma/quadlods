/* Copyright 2017-2018 Pierre Abbat.
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
#include "circletest.h"
#include "ps.h"
#include "hstep.h"

using namespace std;

papersize a4land={297,210};
const double criterion=21.25;
/* criterion is set to the geometric mean of the relative error of 28901 and
 * 1327, which reaches 7.53 shortly before 1M and then turns back down, and
 * the relative error of 65027 and 28901, which reaches almost -60 at 1M and
 * keeps going down.
 */

void circletest(quadlods &quad)
/* Find all pairs of primes (dimensions) which have high discrepancy by
 * counting the pairs of points that fall inside a quadrant of a circle.
 * The number of points in the quadrant should be (π/4)i+O(ln(i)²/i).
 */
{
  int i,j,k,inx,iters=1048576,allinx;
  char buf[24];
  set<int> halfsteps=hsteps(iters);
  time_t now,then;
  bool recordthis;
  quadlods sel2;
  vector<int> pinx2;
  vector<int> incircle,ri;
  vector<double> point;
  vector<vector<double> > rrelError;
  vector<vector<int> > primepairs;
  vector<errorrec> errorrecs;
  PostScript ps;
  double relativeError,maxError,scale;
  allinx=quad.size()*(quad.size()-1)/2;
  incircle.resize(quad.size()*(quad.size()-1)/2);
  rrelError.resize(quad.size()*(quad.size()-1)/2);
  primepairs.resize(quad.size()*(quad.size()-1)/2);
  for (j=0;j<quad.size();j++)
    for (k=0;k<j;k++)
    {
      inx=j*(j-1)/2+k;
      //primepairs[inx].push_back(quad.getprime(k));
      //primepairs[inx].push_back(quad.getprime(j));
      pinx2.clear();
      pinx2.push_back(j);
      pinx2.push_back(k);
      sel2=select(quad,pinx2);
      for (i=0;i<=iters;i++)
      {
        recordthis=halfsteps.count(i);
        point=sel2.dgen();
        if (inx>=errorrecs.size())
        {
          errorrecs.resize(inx+1);
          errorrecs[inx].primepair[0]=sel2.getprime(0);
          errorrecs[inx].primepair[1]=sel2.getprime(1);
        }
        if (point[0]*point[0]+point[1]*point[1]<1)
          incircle[inx]++;
        if (point[0]*point[0]+(1-point[1])*(1-point[1])<1)
          incircle[inx]--;
        if (recordthis)
        {
          relativeError=(incircle[inx])/log(i+1)/log(i+1);
          if (i==0)
            relativeError=0;
          rrelError[inx].push_back(relativeError);
          errorrecs[inx].relError.push_back(relativeError);
          if (ri.size()==0 || ri.back()<i+1)
            ri.push_back(i+1);
        }
      }
      now=time(nullptr);
      if (now!=then)
      {
        cout<<rint((double)inx/allinx*100)<<"% \r";
        cout.flush();
        then=now;
      }
      errorrecs[inx].maxError=0;
      for (i=0;i<errorrecs[inx].relError.size();i++)
        if (fabs(errorrecs[inx].relError[i])>errorrecs[inx].maxError)
          errorrecs[inx].maxError=fabs(errorrecs[inx].relError[i]);
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
  for (i=0;i<errorrecs.size();i++)
  {
    /*for (maxError=j=0;j<ri.size();j++)
      if (fabs(rrelError[i][j])>maxError)
        maxError=fabs(rrelError[i][j]);*/
    ps.startpage();
    ps.setscale(0,-1,3,1);
    scale=errorrecs[i].maxError;
    ps.startline();
    ps.lineto(0,-1);
    ps.lineto(3,-1);
    ps.lineto(3,1);
    ps.lineto(0,1);
    ps.endline(true);
    sprintf(buf,"%g",scale);
    ps.write(3,1,buf);
    sprintf(buf,"%d %d",errorrecs[i].primepair[0],errorrecs[i].primepair[1]);
    ps.write(0,1,buf);
    ps.startline();
    for (j=0;j<ri.size();j++)
      ps.lineto(log(ri[j])/log(iters)*3,errorrecs[i].relError[j]/scale);
    ps.endline();
    ps.endpage();
  }
}
