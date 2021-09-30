/******************************************************/
/*                                                    */
/* filltest.cpp - test how well numbers fill space    */
/*                                                    */
/******************************************************/
/* Copyright 2018-2021 Pierre Abbat.
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
#include "plot.h"
using namespace std;
using namespace quadlods;

double rootBallVolume(int n)
/* Returns the nth root of the volume of a unit n-ball. For n greater than 200
 * or so, computing the volume of a ball between the LD points would naïvely
 * consist of multiplying the nth power of a distance, which could be 14 or more,
 * by the volume of the unit n-ball, which underflows. This avoids the underflow
 * as long as possible.
 */
{
  return sqrt(M_PI)/exp(lgamma(n*0.5+1)/n);
}

double distsq(vector<double> a,vector<double> b)
// Computes the distance with opposite faces identified.
{
  vector<double> d;
  int i;
  double d1,d2;
  assert(a.size()==b.size());
  for (i=0;i<a.size();i++)
  {
    d1=a[i]-b[i];
    if (a[i]>b[i])
      d2=(a[i]-1)-b[i];
    else
      d2=a[i]-(b[i]-1);
    if (fabs(d1)>fabs(d2))
      d1=d2;
    d.push_back(sqr(d1));
  }
  return pairwisesum(d);
}

/* Filltest works like this: For an n-dimensional generator, pick n random
 * points in n-space. (The points have coordinates equal to (c+0.5)/256, where
 * c is a random byte.) For each of these points, remember the vector to the
 * closest generated point so far. The n vectors form a square matrix. Its
 * determinant should decrease at a known rate; the determinant of the
 * normalized vectors should have a known distribution depending only on n.
 *
 * Around the n random points in n-space, put the largest open ball that does
 * not contain any generated point. Its volume should be about 1/i, where i
 * is the number of points generated so far. The second graph is of i times
 * the average of the volumes of the n balls. It should end close to 1.
 *
 * The third graph is the determinant, normalized so that for low-discrepancy
 * sequences, it should end close to 1. It indicates that the vectors from the
 * n random points to the closest generated points are pointing in all directions.
 * If the points form planes, as happens with the unscrambled sequences from
 * the twin primes 3251,3253,13691,13693,21611,21613,59051,59053,65027,65029
 * in both Richtmyer and Halton, the vectors will be close to perpendicular to
 * the planes, and the determinant will be small.
 */
void filltest(Quadlods &quad,int iters,PostScript &ps)
{
  int i,j,k,l,sz=quad.size(),decades,byte;
  char buf[24];
  set<int>::iterator it;
  array<vector<vector<double> >,3> points,disp;
  array<vector<double>,3> closedist;
  vector<double> point;
  vector<double> detGraph,ballGraph,normGraph;
  double hi=-INFINITY,lo=INFINITY,bhi=-INFINITY,blo=INFINITY,nhi=-INFINITY,nlo=INFINITY;
  double thisdist,scale,ballvol,detsqsum,ballsqsum,normsqsum;
  double rbv=rootBallVolume(sz);
  matrix actualSize(sz,sz),normalized(sz,sz);
  set<int> halfsteps=hsteps(1,iters);
  vector<int> halfstepsv;
  manysum weights,reldets,balls;
  time_t now,then;
  for (it=halfsteps.begin();it!=halfsteps.end();++it)
    halfstepsv.push_back(*it);
  for (k=0;k<3;k++)
    while (points[k].size()<sz)
    { // Select n random points in n-space without replacement, three times
      i=points[k].size();
      points[k].resize(i+1);
      disp[k].resize(i+1);
      for (j=0;j<sz;j++)
      {
	byte=rng.ucrandom();
	/* Use binary fractions when testing Richtmyer and quadratic irrationals
	 * when testing Halton. Using binary fractions when testing 1D Halton
	 * with base=2 will result in hitting the random points within
	 * 512 iterations.
	 */
	if (quad.getMode()==QL_MODE_RICHTMYER)
	  points[k][i].push_back((byte+0.5)/256);
	else
	  points[k][i].push_back(nthquad(byte,true));
	disp[k][i].push_back(1);
      }
      for (j=0;j<i;j++)
	if (distsq(points[k][i],points[k][j])==0)
	{
	  points[k].resize(i);
	  disp[k].resize(i);
	  break;
	}
    }
  for (i=0;i<sz;i++)
    for (l=0;l<3;l++)
      closedist[l].push_back(sz);
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
    for (l=0;l<3;l++)
      for (j=0;j<sz;j++)
      {
	thisdist=distsq(point,points[l][j]);
	if (thisdist<closedist[l][j])
	{
	  closedist[l][j]=thisdist;
	  for (k=0;k<sz;k++)
	    disp[l][j][k]=point[k]-points[l][j][k];
	}
      }
    if (halfsteps.count(i))
    {
      detsqsum=ballsqsum=normsqsum=0;
      for (l=0;l<3;l++)
      {
	for (j=0,ballvol=0;j<sz;j++)
	{
	  for (k=0;k<sz;k++)
	  {
	    actualSize[j][k]=disp[l][j][k];
	    normalized[j][k]=disp[l][j][k]*sqrt(sz/(j+1.)/closedist[l][j]);
	  }
	  ballvol+=i*pow(closedist[l][j]*rbv,sz*0.5);
	}
	detsqsum+=sqr(actualSize.determinant());
	ballsqsum+=sqr(ballvol/sz);
	normsqsum+=sqr(normalized.determinant());
      }
      detGraph.push_back(log(detsqsum/3)/2);
      ballGraph.push_back(log(ballsqsum/3)/2);
      normGraph.push_back(log(normsqsum/3)/2);
      weights+=i;
      reldets+=i*sqr(i)*detsqsum/3;
      balls+=i*sqrt(ballsqsum/3);
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
  for (i=0;i<ballGraph.size();i++)
  {
    if (ballGraph[i]>bhi)
      bhi=ballGraph[i];
    if (ballGraph[i]<blo)
      blo=ballGraph[i];
  }
  bhi=ceil (bhi/log(10))*log(10);
  blo=floor(blo/log(10))*log(10);
  for (i=0;i<normGraph.size();i++)
  {
    if (normGraph[i]>nhi)
      nhi=normGraph[i];
    if (normGraph[i]<nlo)
      nlo=normGraph[i];
  }
  if (nhi==nlo)
  {
    nhi+=0.1;
    nlo-=0.1;
  }
  nhi=ceil (nhi/log(10))*log(10);
  nlo=floor(nlo/log(10))*log(10);
  logLogPlot(ps,halfstepsv,detGraph);
  logLogPlot(ps,halfstepsv,ballGraph);
  logLogPlot(ps,halfstepsv,normGraph);
  //cout<<"Average relative determinant "<<reldets.total()/weights.total()<<endl;
  cout<<"Average total ball volume "<<balls.total()/weights.total()<<endl;
}
