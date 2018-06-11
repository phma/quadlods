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
  int i,j,k,sz=quad.size();
  vector<vector<double> > points,disp;
  vector<double> closedist;
  while (points.size()<sz)
  {
    i=points.size();
    points.resize(i+1);
    for (j=0;j<sz;j++)
      points[i].push_back((rng.ucrandom()+0.5)/256);
    for (j=0;j<i;j++)
      if (distsq(points[i],points[j])==0)
      {
	points.resize(i);
	break;
      }
  }
  for (i=0;i<sz;i++)
    closedist.push_back(sz);
}
