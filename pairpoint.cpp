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

void PairCompressor::insert(const xy &pnt)
{
  layers[0].insert(pnt,0);
}
