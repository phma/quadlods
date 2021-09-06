/******************************************************/
/*                                                    */
/* discrepancy.cpp - compute discrepancy              */
/*                                                    */
/******************************************************/
/* Copyright 2020,2021 Pierre Abbat.
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
#include <cmath>
#include <cassert>
#include <iostream>
#include <chrono>
#include "discrepancy.h"
#include "quadlods.h"
#include "random.h"
#include "threads.h"
#include "dotbaton.h"

using namespace std;
using namespace quadlods;
namespace cr=std::chrono;

vector<Box> emptyPop;
vector<Box> population;
BoxCountBlock boxCountBlock;
double flowerDisc[2]={0,1};
/* When computing the discrepancy of a flower plot, these changes apply:
 * • This array is set to {-1,1}.
 * • The number of dimensions is two.
 * • The lower limit of bounds is -1.
 * • Areas are clipped to the unit circle.
 */

double clipToCircle(double x0,double x1,double y)
{
  if (x0*x0+y*y>1)
  {
    double x2=sqrt(1-y*y);
    if (x2*x0<0)
      x2=-x2;
    if ((x1-x2)*(x2-x0)>=0)
      x0=x2;
    else if (fabs(x0)>fabs(x1))
      x0=x1;
  }
  return x0;
}

double sectorArea(double x0,double y0,double x1,double y1)
{
  if (x0+x1<0) // avoid branch cut of atan2
  {
    x0=-x0;
    y0=-y0;
    x1=-x1;
    y1=-y1;
  }
  return (atan2(y1,x1)-atan2(y0,x0))/2;
}

double areaInCircle(double minx,double miny,double maxx,double maxy)
/* Computes the area of the part of a rectangle that is inside the unit circle.
 * Used when computing the discrepancy of a flower plot.
 */
{
  double maxylef,minylef,minxbot,maxxbot,minyrig,maxyrig,maxxtop,minxtop;
  double leftri,bottri,rigtri,toptri;
  double blsect,brsect,trsect,tlsect;
  double negsect=0,possect=0,ret;
  maxylef=clipToCircle(maxy,miny,minx);
  minylef=clipToCircle(miny,maxy,minx);
  minxbot=clipToCircle(minx,maxx,miny);
  maxxbot=clipToCircle(maxx,minx,miny);
  minyrig=clipToCircle(miny,maxy,maxx);
  maxyrig=clipToCircle(maxy,miny,maxx);
  maxxtop=clipToCircle(maxx,minx,maxy);
  minxtop=clipToCircle(minx,maxx,maxy);
  leftri=minx*(minylef-maxylef)/2;
  bottri=miny*(minxbot-maxxbot)/2;
  rigtri=maxx*(maxyrig-minyrig)/2;
  toptri=maxy*(maxxtop-minxtop)/2;
  blsect=sectorArea(minx,minylef,minxbot,miny);
  brsect=sectorArea(maxxbot,miny,maxx,minyrig);
  trsect=sectorArea(maxx,maxyrig,maxxtop,maxy);
  tlsect=sectorArea(minxtop,maxy,minx,maxylef);
  if (blsect>0)
    possect+=blsect;
  else
    negsect+=blsect;
  if (brsect>0)
    possect+=brsect;
  else
    negsect+=brsect;
  if (trsect>0)
    possect+=trsect;
  else
    negsect+=trsect;
  if (tlsect>0)
    possect+=tlsect;
  else
    negsect+=tlsect;
  ret=(leftri+bottri+rigtri+toptri)+(possect+negsect);
  if (fabs(ret)<=DBL_EPSILON)
    ret=0;
  assert(ret>=0);
  return ret;
}

void setFlowerDisc(bool fd)
{
  if (fd)
    flowerDisc[0]=-1;
  else
    flowerDisc[0]=0;
  flowerDisc[1]=1;
}

Box::Box()
{
  pointsIn=pointsBound=pointsTotal=0;
  volume=NAN;
}

Box::Box(vector<double> pnt0,vector<double> pnt1)
{
  int i;
  assert(pnt0.size()==pnt1.size());
  pointsIn=pointsBound=pointsTotal=0;
  volume=NAN;
  bounds.resize(pnt0.size());
  for (i=0;i<pnt0.size();i++)
  {
    bounds[i][0]=pnt0[i];
    bounds[i][1]=pnt1[i];
    if (bounds[i][0]>bounds[i][1])
      swap(bounds[i][0],bounds[i][1]);
  }
}

Box::Box(Box &mother,Box &father)
{
  int i;
  pointsIn=pointsBound=pointsTotal=0;
  volume=NAN;
  for (i=0;i<mother.bounds.size() && i<father.bounds.size();i++)
    if (rng.brandom())
      bounds.push_back(father.bounds[i]);
    else
      bounds.push_back(mother.bounds[i]);
  bounds.shrink_to_fit();
}

int Box::in(const vector<double> &point)
// Returns 2 if inside, 1 if on the boundary, 0 if outside.
{
  int i,ret=2;
  if (point.size()!=bounds.size())
    throw sizeMismatch;
  for (i=0;i<point.size() && ret;i++)
  {
    if (bounds[i][0]==point[i] || bounds[i][1]==point[i])
      ret=1;
    if (bounds[i][0]>point[i] || bounds[i][1]<point[i])
      ret=0;
  }
  return ret;
}

void Box::countPoints(const vector<vector<double> > &points)
{
  int i,ptin;
  volume=1;
  if (flowerDisc[0])
    volume=areaInCircle(bounds[0][0],bounds[1][0],bounds[0][1],bounds[1][1])/M_PI;
  else
    for (i=0;i<bounds.size();i++)
      volume*=bounds[i][1]-bounds[i][0];
  pointsIn=pointsBound=0;
  pointsTotal=points.size();
  for (i=0;i<points.size();i++)
  {
    ptin=in(points[i]);
    if (ptin==1)
      pointsBound++;
    if (ptin==2)
      pointsIn++;
  }
}

double Box::discrepancy()
/* Returns the signed discrepancy including or excluding the boundary,
 * whichever is larger in absolute value.
 */
{
  double opendisc,closedisc;
  opendisc=(double)pointsIn/pointsTotal-volume;
  closedisc=(double)(pointsIn+pointsBound)/pointsTotal-volume;
  if (fabs(opendisc)>fabs(closedisc))
    return opendisc;
  else
    return closedisc;
}

void Box::mutate(const std::vector<std::vector<double> > &points,int pntnum,int coord)
/* Replaces one of the bounds, at random, with the corresponding coordinate
 * of one of the points, at random.
 */
{
  if (pntnum<0)
    pntnum=rng.rangerandom(points.size()+2);
  if (coord<0)
    coord=rng.rangerandom(bounds.size());
  bounds[coord][rng.brandom()]=(pntnum<points.size())?points[pntnum][coord]:(flowerDisc[pntnum-points.size()]);
  if (bounds[coord][0]>bounds[coord][1])
    swap(bounds[coord][0],bounds[coord][1]);
}

void Box::dump()
{
  int i;
  printf("%5.3f %d,%d/%d",volume,pointsIn,pointsBound,pointsTotal);
  for (i=0;i<bounds.size();i++)
    printf(" [%5.3f,%5.3f]",bounds[i][0],bounds[i][1]);
  printf("\n");
}

bool operator==(const Box &l,const Box &r)
{
  int i;
  bool ret=l.bounds.size()==r.bounds.size();
  for (i=0;ret && i<l.bounds.size();i++)
    ret=l.bounds[i][0]==r.bounds[i][0] && l.bounds[i][1]==r.bounds[i][1];
  return ret;
}

BoxCountBlock::BoxCountBlock()
{
  pop=&emptyPop;
  pts=nullptr;
  b=e=left=0;
}

void BoxCountBlock::load(vector<Box> &population,int begin,int end,const vector<vector<double> > &points)
{
  mtx.lock();
  pop=&population;
  b=begin;
  e=end;
  left=e-b;
  pts=&points;
  mtx.unlock();
}

BoxCountItem BoxCountBlock::getItem()
{
  BoxCountItem ret{&(*pop)[b],*pts};
  mtx.lock();
  if (b<e)
    ret.box=&(*pop)[b++];
  else
    ret.box=nullptr;
  mtx.unlock();
  return ret;
}

void BoxCountBlock::countFinished()
{
  mtx.lock();
  left--;
  mtx.unlock();
}

bool BoxCountBlock::done()
{
  bool ret;
  mtx.lock();
  ret=left==0 && b==e;
  mtx.unlock();
  return ret;
}

void sort(vector<Box> &pop,int begin,int end,int popLimit)
/* Quicksort, but with some recursion omitted.
 * If begin=0 or popLimit is in the interval, sort; else don't bother.
 * The 0th element is the most discrepant, and everything before popLimit
 * is more discrepant than everything after. Other than that, order
 * doesn't matter.
 */
{
  int i,j;
  double pivot;
  if ((begin==0 || (begin<=popLimit && end>=popLimit)) && end-begin>1)
  {
    i=begin;
    j=end-1;
    pivot=fabs(pop[(i+j)/2].discrepancy());;
    while (i<j)
    {
      while (fabs(pop[i].discrepancy())>pivot)
        i++;
      while (fabs(pop[j].discrepancy())<pivot)
        j--;
      if (i<j)
      {
        swap(pop[i],pop[j]);
        if (fabs(pop[i].discrepancy())==fabs(pop[j].discrepancy()))
          if (rng.brandom())
            i++;
          else
            j--;
      }
    }
    i=(i+j)/2;
    sort(pop,begin,i,popLimit);
    sort(pop,i+1,end,popLimit);
  }
}

void shuffle(vector<Box> &pop)
{
  int i;
  for (i=pop.size();i>1;i-=2)
    swap(pop[i-1],pop[rng.rangerandom(i)]);
}

void dumpdisc(vector<Box> &pop)
{
  int i;
  for (i=0;i<pop.size();i++)
  {
    cout<<i<<' '<<pop[i].discrepancy()<<"  ";
    if (i==pop.size()-1 || i%5==4)
      cout<<endl;
  }
}

double prog(int nsteady,int niter)
// Decreases to 0 as progress is made.
{
  int endpt=niter/3+20;
  return ((double)endpt-nsteady)/(endpt+0.5);
}

double discrepancy(const vector<vector<double> > &points,bool keepPop)
/* Computes the discrepancy (or a lower bound) of the points. keepPop is for
 * incrementally computing the discrepancy of a long list of points. The next
 * call will assume that all points up to Box::pointsTotal are the same as
 * in this call.
 */
{
  vector<int> delenda;
  DotBaton dotbaton;
  mpq_class mutationRate(1,points[0].size());
  double lastdisc=-1;
  int i,j,prevsz=0,sz,dim,nParents,popLimit,niter=0,nsteady=0;
  vector<double> all0,all1;
  cr::nanoseconds elapsed;
  cr::time_point<cr::steady_clock> timeStart;
  if (population.size())
    prevsz=population[0].getPointsTotal();
  sz=points.size();
  dim=points[0].size();
  popLimit=3*dim*sz+8192;
  for (i=0;i<dim;i++)
  {
    all0.push_back(flowerDisc[0]);
    all1.push_back(flowerDisc[1]);
  }
  for (i=0;i<sz;i++)
    population.push_back(Box(points[i],points[(i+1)%sz]));
  population.push_back(Box(all0,all1));
  for (i=0;i<sz*2;i++)
    population.push_back(Box(points[i%sz],points[rng.rangerandom(sz)]));
  for (i=0;i<sz;i++)
    for (j=0;j<dim;j++)
    {
      if (prevsz)
	population.push_back(population[rng.rangerandom(prevsz)]);
      else
	population.push_back(Box(all0,all1));
      population.back().mutate(points,i,j);
    }
  boxCountBlock.load(population,0,population.size(),points);
  while (!boxCountBlock.done())
  {
    this_thread::sleep_for(chrono::milliseconds(1));
    dotbaton.update(1e-7,boxCountBlock.getLeft());
  }
  while (prog(nsteady,niter) || population.size()<popLimit)
  {
    timeStart=clk.now();
    shuffle(population);
    nParents=population.size();
    for (i=0;i<nParents;i+=2)
      population.push_back(Box(population[i],population[i+1]));
    for (i=nParents;i<population.size();i++)
      if (rng.frandom(mutationRate))
        population[i].mutate(points);
    elapsed=clk.now()-timeStart;
    //cout<<"Breeding took "<<elapsed.count()/1e6<<" ms\n";
    timeStart=clk.now();
    boxCountBlock.load(population,nParents,population.size(),points);
    while (!boxCountBlock.done())
    {
      this_thread::sleep_for(chrono::milliseconds(1));
      dotbaton.update(prog(nsteady,niter),boxCountBlock.getLeft());
    }
    elapsed=clk.now()-timeStart;
    //cout<<population.size()-nParents<<" new boxes took "<<elapsed.count()/1e6<<" ms\n";
    sort(population,0,population.size(),popLimit);
    delenda.clear();
    sz=population.size();
    for (i=0;i<popLimit-1 && i<sz-1;i++)
      if (population[i]==population[i+1])
	delenda.push_back(i+1); // Eliminate duplicates
    for (i=0;i<delenda.size();i++)
      if (popLimit+i<sz)
	swap(population[delenda[i]],population[i+popLimit]);
    if (sz>popLimit)
      population.resize(popLimit);
    niter++;
    if (lastdisc==population[0].discrepancy())
      nsteady++;
    else
    {
      lastdisc=population[0].discrepancy();
      //cout<<"iter "<<niter<<" disc "<<lastdisc<<endl;
      nsteady=0;
    }
  }
  dotbaton.update(0,0);
  for (i=0;i>3 && i<population.size();i++)
    population[i].dump();
  if (keepPop)
    population.resize(3);
  else
    population.clear();
  return fabs(population[0].discrepancy());
}

bool countAnyBlock()
{
  BoxCountItem item=boxCountBlock.getItem();
  if (item.box)
  {
    item.box->countPoints(item.points);
    boxCountBlock.countFinished();
  }
  return item.box!=nullptr;
}
