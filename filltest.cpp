/******************************************************/
/*                                                    */
/* filltest.cpp - test how well numbers fill space    */
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

#include <cassert>
#include "filltest.h"
#include "hstep.h"
#include "manysum.h"
#include "matrix.h"
#include "histogram.h"
#include "random.h"
using namespace std;

double rootBallVolume(int n)
/* Returns the nth root of the volume of a unit n-ball. For n greater than 200
 * or so, computing the volume of a ball between the LD points would naïvely
 * consist of raising the nth power of a distance, which could be up to 14,
 * by the volume of the unit n-ball, which underflows. This avoids the underflow
 * as long as possible.
 */
{
  return sqrt(M_PI)/exp(lgamma(n*0.5+1)/n);
}

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
void filltest(Quadlods &quad,int iters,PostScript &ps)
{
  int i,j,k,l,sz=quad.size(),decades;
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
  manysum weights,reldets,balls;
  time_t now,then;
  for (k=0;k<3;k++)
    while (points[k].size()<sz)
    {
      i=points[k].size();
      points[k].resize(i+1);
      disp[k].resize(i+1);
      for (j=0;j<sz;j++)
      {
	points[k][i].push_back((rng.ucrandom()+0.5)/256);
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
  cout<<"halfsteps size "<<halfsteps.size()<<" detGraph size "<<detGraph.size()<<endl;
  ps.startline();
  for (it=halfsteps.begin(),i=0;it!=halfsteps.end();i++,it++)
    if (*it)
      ps.lineto(log(*it)/log(iters)*3,(detGraph[i]-lo)/scale-1);
  ps.endline();
  ps.endpage();
  ps.startpage();
  ps.setscale(0,-1,3,1);
  scale=(bhi-blo)/2;
  ps.startline();
  ps.lineto(0,-1);
  ps.lineto(3,-1);
  ps.lineto(3,1);
  ps.lineto(0,1);
  ps.endline(true);
  decades=rint((bhi-blo)/log(10));
  for (i=0;i<=decades;i++)
  {
    sprintf(buf,"%g",exp(blo)*pow(10,i));
    ps.write(3.1,i*2./decades-1,buf);
    ps.startline();
    ps.lineto(3,i*2./decades-1);
    ps.lineto(3.1,i*2./decades-1);
    ps.endline();
  }
  xticks(1,iters,ps);
  ps.startline();
  for (it=halfsteps.begin(),i=0;it!=halfsteps.end();i++,it++)
    if (*it)
      ps.lineto(log(*it)/log(iters)*3,(ballGraph[i]-blo)/scale-1);
  ps.endline();
  ps.endpage();
  ps.startpage();
  ps.setscale(0,-1,3,1);
  scale=(nhi-nlo)/2;
  ps.startline();
  ps.lineto(0,-1);
  ps.lineto(3,-1);
  ps.lineto(3,1);
  ps.lineto(0,1);
  ps.endline(true);
  decades=rint((nhi-nlo)/log(10));
  for (i=0;i<=decades;i++)
  {
    sprintf(buf,"%g",exp(nlo)*pow(10,i));
    ps.write(3.1,i*2./decades-1,buf);
    ps.startline();
    ps.lineto(3,i*2./decades-1);
    ps.lineto(3.1,i*2./decades-1);
    ps.endline();
  }
  xticks(1,iters,ps);
  ps.startline();
  for (it=halfsteps.begin(),i=0;it!=halfsteps.end();i++,it++)
    if (*it)
      ps.lineto(log(*it)/log(iters)*3,(normGraph[i]-nlo)/scale-1);
  ps.endline();
  ps.endpage();
  //cout<<"Average relative determinant "<<reldets.total()/weights.total()<<endl;
  cout<<"Average total ball volume "<<balls.total()/weights.total()<<endl;
}
