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
#include "pairpoint.h"
using namespace std;

void Layer::insert(const xy &pnt,int32_t n)
{
  PairDot pairDot;
  pairDot.location=pnt;
  pairDot.inx=n;
  dots[next++]=pairDot;
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

void PairCompressor::findNewPair(int layerNum)
/* Finds two pairs of points with the same difference in the current layer.
 * The last point added, which l will point to, must not differ by a known
 * difference from another point in the same layer. All points, if pairs,
 * must represent identical patterns (represented by inx). Replaces the
 * pairs found in this layer by single points in the next layer, unless
 * two of them share a point, in which case it replaces one pair.
 */
{
  map<int64_t,PairDot>::iterator i,j,k,l;
  xy twidiff,diff;
  bool found=false;
  l=layers[layerNum].dots.end();
  --l;
  for (k=layers[layerNum].dots.begin();k!=l;++k)
    for (j=layers[layerNum].dots.begin();k->second.inx==l->second.inx && j!=k;++j)
      for (i=layers[layerNum].dots.begin();j->second.inx==l->second.inx && i!=l;++i)
	if (i->second.inx==l->second.inx)
	{
	  diff=l->second.location-k->second.location;
	  twidiff=diff-j->second.location+i->second.location;
	  if (diff.length()>1e-4 && twidiff.length()<1e-6)
	  {
	    found=true;
	    goto foundit;
	  }
	}
  foundit:
  diff=(diff-i->second.location+j->second.location)/2;
  // ...
}

void PairCompressor::insert(const xy &pnt)
{
  layers[0].insert(pnt,0);
}
