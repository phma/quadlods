/******************************************************/
/*                                                    */
/* filltest.cpp - test how well numbers fill space    */
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

#include <cassert>
#include "filltest.h"
#include "hstep.h"
#include "manysum.h"
#include "matrix.h"
#include "histogram.h"
#include "random.h"
#include "ps.h"
using namespace std;

double distsq(vector<double> a,vector<double> b)
{
  vector<double> d;
  int i;
  assert(a.size()==b.size());
  for (i=0;i<a.size();i++)
    d.push_back(sqr(a[i]-b[i]));
  return pairwisesum(d);
}

/* Filltest works like this: For an n-dimensional generator, pick n random
 * points in n-space. (The points have coordinates equal to (c+0.5)/256, where
 * c is a random byte.) For each of these points, remember the vector to the
 * closest generated point so far. The n vectors form a square matrix. Its
 * determinant should decrease at a known rate; the determinant of the
 * normalized vectors should have a known distribution depending only on n.
 */
void filltest(quadlods &quad)
{
  int i,j,k,sz=quad.size(),iters=1048576,decades;
  char buf[24];
  set<int>::iterator it;
  vector<vector<double> > points,disp;
  vector<double> closedist,point;
  vector<double> detGraph;
  double thisdist,hi=-INFINITY,lo=INFINITY,scale;
  matrix actualSize(sz,sz),normalized(sz,sz);
  set<int> halfsteps=hsteps(iters);
  time_t now,then;
  PostScript ps;
  while (points.size()<sz)
  {
    i=points.size();
    points.resize(i+1);
    disp.resize(i+1);
    for (j=0;j<sz;j++)
    {
      points[i].push_back((rng.ucrandom()+0.5)/256);
      disp[i].push_back(1);
    }
    for (j=0;j<i;j++)
      if (distsq(points[i],points[j])==0)
      {
	points.resize(i);
	disp.resize(i);
	break;
      }
  }
  for (i=0;i<sz;i++)
    closedist.push_back(sz);
  ps.open("filltest.ps");
  ps.setpaper(a4land,0);
  ps.prolog();
  for (i=0;i<=iters;i++)
  {
    now=time(nullptr);
    if (now!=then)
    {
      cout<<rint((double)i/iters*100)<<"% \r";
      cout.flush();
      then=now;
    }
    point=quad.dgen();
    for (j=0;j<sz;j++)
    {
      thisdist=distsq(point,points[j]);
      if (thisdist<closedist[j])
      {
	closedist[j]=thisdist;
	for (k=0;k<sz;k++)
	  disp[j][k]=point[k]-points[j][k];
      }
    }
    if (halfsteps.count(i))
    {
      for (j=0;j<sz;j++)
	for (k=0;k<sz;k++)
	{
	  actualSize[j][k]=disp[j][k];
	  normalized[j][k]=disp[j][k]*sqrt(j/closedist[j]);
	}
      detGraph.push_back(log(fabs(actualSize.determinant())));
    }
  }
  for (i=0;i<detGraph.size();i++)
  {
    if (detGraph[i]>hi)
      hi=detGraph[i];
    if (detGraph[i]<lo)
      lo=detGraph[i];
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
  cout<<"halfsteps size "<<halfsteps.size()<<" detGraph size "<<detGraph.size()<<endl;
  ps.startline();
  for (it=halfsteps.begin(),i=0;it!=halfsteps.end();i++,it++)
    if (*it)
      ps.lineto(log(*it)/log(iters)*3,(detGraph[i]-lo)/scale-1);
  ps.endline();
  ps.endpage();
}
