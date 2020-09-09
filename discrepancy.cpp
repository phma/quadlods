/******************************************************/
/*                                                    */
/* discrepancy.cpp - compute discrepancy              */
/*                                                    */
/******************************************************/
/* Copyright 2020 Pierre Abbat.
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
#include "discrepancy.h"
#include "random.h"

using namespace std;

const mpq_class mutationRate(1,30);

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

void Box::mutate(const std::vector<std::vector<double> > &points)
/* Replaces one of the bounds, at random, with the corresponding coordinate
 * of one of the points, at random.
 */
{
  int pntnum=rng.rangerandom(points.size()+2);
  int coord=rng.rangerandom(bounds.size());
  bounds[coord][rng.brandom()]=(pntnum<points.size())?points[pntnum][coord]:(pntnum-points.size());
  if (bounds[coord][0]>bounds[coord][1])
    swap(bounds[coord][0],bounds[coord][1]);
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
        if (pop[i].discrepancy()==pop[j].discrepancy())
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

double discrepancy(const vector<vector<double> > &points)
{
  vector<Box> population;
  double lastdisc=-1;
  int i,sz,nParents,popLimit,niter=0,nsteady=0;
  vector<double> all0,all1;
  sz=points.size();
  popLimit=5*sz+256;
  for (i=0;i<points[0].size();i++)
  {
    all0.push_back(0);
    all1.push_back(1);
  }
  for (i=0;i<sz;i++)
    population.push_back(Box(points[i],points[(i+1)%sz]));
  population.push_back(Box(all0,all1));
  for (i=0;i<=sz;i++)
    population[i].countPoints(points);
  while (nsteady<niter/2+10 || population.size()<popLimit)
  {
    shuffle(population);
    nParents=population.size();
    for (i=0;i<nParents;i+=2)
      population.push_back(Box(population[i],population[i+1]));
    for (i=nParents;i<population.size();i++)
    {
      if (rng.frandom(mutationRate))
        population[i].mutate(points);
      population[i].countPoints(points);
    }
    sort(population,0,population.size(),popLimit);
    if (population.size()>popLimit)
      population.resize(popLimit);
    niter++;
    if (lastdisc==population[0].discrepancy())
      nsteady++;
    else
    {
      lastdisc=population[0].discrepancy();
      cout<<"iter "<<niter<<" disc "<<lastdisc<<endl;
    }
  }
  return population[0].discrepancy();
}
