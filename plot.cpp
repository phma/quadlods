/******************************************************/
/*                                                    */
/* plot.cpp - plot graphs                             */
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

#include "plot.h"
#include "hstep.h"
using namespace std;

void logLogPlot(PostScript &ps,vector<int> xcoords,vector<double> ycoords)
// ycoords are logarithms.
{
  int i,j,k,l,sz=xcoords.size(),decades,byte,iters=0;
  char buf[24];
  double hi=-INFINITY,lo=INFINITY,scale;
  for (i=0;i<sz;i++)
  {
    if (ycoords[i]>hi)
      hi=ycoords[i];
    if (ycoords[i]<lo)
      lo=ycoords[i];
    if (xcoords[i]>iters)
      iters=xcoords[i];
  }
  hi=ceil (hi/log(10))*log(10);
  lo=floor(lo/log(10))*log(10);
  ps.startpage();
  ps.setscale(0,-1,3,1);
  scale=(hi-lo)/2;
  ps.startline();
  ps.lineto(0,-1);
  ps.lineto(3,-1);
  ps.lineto(3,1);
  ps.lineto(0,1);
  ps.endline(true);
  decades=rint((hi-lo)/log(10));
  for (i=0;i<=decades;i++)
  {
    sprintf(buf,"%g",exp(lo)*pow(10,i));
    ps.write(3.1,i*2./decades-1,buf);
    ps.startline();
    ps.lineto(3,i*2./decades-1);
    ps.lineto(3.1,i*2./decades-1);
    ps.endline();
  }
  xticks(1,iters,ps);
  ps.startline();
  for (i=0;i<sz;i++)
    ps.lineto(log(xcoords[i])/log(iters)*3,(ycoords[i]-lo)/scale-1);
  ps.endline();
  ps.endpage();
}

void logLogPlot(PostScript &ps,set<int> xcoords,vector<double> ycoords)
{
  vector<int> xcoordsv;
  set<int>::iterator it;
  for (it=xcoords.begin();it!=xcoords.end();++it)
    xcoordsv.push_back(*it);
  logLogPlot(ps,xcoordsv,ycoords);
}

