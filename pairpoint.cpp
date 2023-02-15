/******************************************************/
/*                                                    */
/* pairpoint.cpp - pairs of points in scatter plot    */
/*                                                    */
/******************************************************/
/* Copyright 2023 Pierre Abbat.
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
#include <iostream>
#include "pairpoint.h"
using namespace std;

int buck(xy pnt)
/* Returns a number from 0 to 399, which is the bucket that pnt goes in.
 * pnt.x and pnt.y are both in (-1,1) and nonzero. Multiplies by π², which
 * is transcendental, to avoid both quadratic irrationals (Richtmyer) and
 * rationals (Halton) falling on the lines between squares. -pnt is in
 * bucket 399-buck(pnt).
 */
{
  return floor(pnt.getx()*M_PI*M_PI)+20*floor(pnt.gety()*M_PI*M_PI)+210;
}

void Layer::insert(const xy &pnt,int32_t n)
{
  PairDot pairDot;
  map<int64_t,PairDot>::iterator i;
  DotDiff dotDiff;
  pairDot.location=pnt;
  pairDot.inx=n;
  dotDiff.a=next;
  for (i=dots.begin();i!=dots.end();++i)
    if (i->second.inx==n)
    {
      dotDiff.b=i->first;
      dotDiff.diff=pnt-i->second.location;
      diffs[buck(dotDiff.diff)].push_back(dotDiff);
    }
  dots[next++]=pairDot;
}

void Layer::cleanBucket(int n)
{
  vector<DotDiff> newBucket;
  int i;
  for (i=0;i<diffs[n].size();i++)
    if (dots.count(diffs[n][i].a) && dots.count(diffs[n][i].b))
      newBucket.push_back(diffs[n][i]);
  swap(diffs[n],newBucket);
}

PairCompressor::PairCompressor()
{
  PairPoint pairPoint;
  pairPoint.level=0;
  pairPoint.sub=-1;
  pairPoint.sep=xy(NAN,NAN);
  pairPoints.push_back(pairPoint);
  layers.push_back(Layer());
}

bool PairCompressor::findOldPair(int layerNum)
{
  map<int64_t,PairDot>::iterator k,l;
  xy diff;
  xy location;
  int i;
  bool found=false;
  l=layers[layerNum].dots.end();
  --l;
  for (k=layers[layerNum].dots.begin();k!=l;++k)
    for (i=0;k->second.inx==l->second.inx && i<pairPoints.size();++i)
      if (pairPoints[i].sub==l->second.inx)
      {
	diff=l->second.location-k->second.location;
	if ((diff-pairPoints[i].sep).length()<1.5e-6 ||
	    (diff+pairPoints[i].sep).length()<1.5e-6)
	  found=true;
	if (found)
	  goto foundit;
      }
  foundit:
  if (found)
  {
    location=(k->second.location+l->second.location-pairPoints[i].sep)/2;
    layers[layerNum+1].insert(location,i);
    layers[layerNum].dots.erase(k);
    layers[layerNum].dots.erase(l);
  }
  return found;
}

bool PairCompressor::findNewPair(int layerNum)
/* Finds two pairs of points with the same difference in the current layer.
 * The last point added, which l will point to, must not differ by a known
 * difference from another point in the same layer. All points, if pairs,
 * must represent identical patterns (represented by inx). Replaces the
 * pairs found in this layer by single points in the next layer, unless
 * two of them share a point, in which case it replaces one pair.
 */
{
  map<int64_t,PairDot>::iterator iit,jit,kit,lit;
  int i,j,k,l;
  PairPoint newPair;
  xy twidiff,diff;
  PairDot newDot;
  bool found=false;
  lit=layers[layerNum].dots.end();
  --lit;
  l=lit->first;
  for (kit=layers[layerNum].dots.begin();kit!=lit;++kit)
    for (jit=layers[layerNum].dots.begin();kit->second.inx==lit->second.inx && jit!=kit;++jit)
      for (iit=layers[layerNum].dots.begin();jit->second.inx==lit->second.inx && iit!=lit;++iit)
	if (iit->second.inx==lit->second.inx)
	{
	  i=iit->first;
	  j=jit->first;
	  k=kit->first;
	  l=lit->first;
	  diff=layers[layerNum].dots[l].location-layers[layerNum].dots[k].location;
	  twidiff=diff-layers[layerNum].dots[j].location+layers[layerNum].dots[i].location;
	  if (diff.length()>1e-4 && twidiff.length()<1e-6)
	  {
	    found=true;
	    goto foundit;
	  }
	}
  foundit:
  if (found)
  {
    diff=(diff-layers[layerNum].dots[i].location+layers[layerNum].dots[j].location)/2;
    newPair.level=layerNum+1;
    newPair.sub=layers[layerNum].dots[l].inx;
    newPair.sep=diff;
    newPair.lowleft=pairPoints[newPair.sub].lowleft+xy(min(diff.getx(),0.),min(diff.gety(),0.));
    newPair.upright=pairPoints[newPair.sub].upright+xy(max(diff.getx(),0.),max(diff.gety(),0.));
    newDot.inx=pairPoints.size();
    pairPoints.push_back(newPair);
    if (layers.size()<=layerNum+1)
      layers.push_back(Layer());
    if (j==k) //FIXME j should go up to k, but no farther, in the loop
    {
      newDot.location=(layers[layerNum].dots[i].location+layers[layerNum].dots[j].location-diff)/2;
      layers[layerNum+1].insert(newDot.location,newDot.inx);
      layers[layerNum].dots.erase(i);
      layers[layerNum].dots.erase(j);
    }
    else
    {
      newDot.location=(layers[layerNum].dots[i].location+layers[layerNum].dots[j].location-diff)/2;
      layers[layerNum+1].insert(newDot.location,newDot.inx);
      newDot.location=(layers[layerNum].dots[k].location+layers[layerNum].dots[l].location-diff)/2;
      layers[layerNum+1].insert(newDot.location,newDot.inx);
      layers[layerNum].dots.erase(i);
      layers[layerNum].dots.erase(j);
      layers[layerNum].dots.erase(k);
      layers[layerNum].dots.erase(l);
    }
  }
  return found;
}

void PairCompressor::insert(const xy &pnt)
{
  int i=0;
  layers[0].insert(pnt,0);
  while (findOldPair(i) || findNewPair(i))
    i++;
}
