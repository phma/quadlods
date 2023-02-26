/******************************************************/
/*                                                    */
/* pairpoint.h - pairs of points in scatter plot      */
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

#ifndef PAIRPOINT_H
#define PAIRPOINT_H

#include <deque>
#include <vector>
#include <map>
#include <cstdint>
#include "xy.h"

#define PBUCKETS 400

struct PairPoint
/* A PairPoint represents a single point, if sep is NaN and sub is -1,
 * or two PairPoints indexed by sub, separated by sep. It is used for
 * compressing scatter plots.
 */
{
  uint32_t level; // 0=single point 3=eight points etc.
  int32_t sub;
  xy sep;
  xy lowleft,upright;
};

struct PairDot
/* A PairDot is an instance of PairPoint at a location in 2-space.
 */
{
  int32_t inx; // index to PairCompressor::pairPoints
  xy location;
};

struct DotDiff
{
  int64_t a,b;
  xy diff; // dots[a]-dots[b]; they must have the same inx
};

class Layer
{
public:
  Layer()
  {
    next=0;
  }
  void insert(const xy &pnt,int32_t n,bool comp=true);
  void cleanBucket(int n);
  std::map<int64_t,PairDot> dots;
  std::vector<DotDiff> diffs[PBUCKETS];
  std::vector<int> pairs[PBUCKETS];
  int64_t next;
};

class PairCompressor
{
public:
  PairCompressor();
  void insert(const xy &pnt);
private:
  bool findOldPair(int layerNum);
  bool findNewPair(int layerNum);
  std::deque<PairPoint> pairPoints;
  std::vector<Layer> layers;
  int compress;
  friend class PostScript;
};

#endif
