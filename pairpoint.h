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
#include <cstdint>
#include "xy.h"

struct PairPoint
/* A PairPoint represents a single point, if sep is NaN and sub is -1,
 * or two PairPoints indexed by sub, separated by sep. It is used for
 * compressing scatter plots.
 */
{
  uint32_t level; // 0=single point 3=eight points etc.
  int32_t sub;
  xy sep;
};

class PairCompressor
{
public:
  PairCompressor();
private:
  std::deque<PairPoint> pairPoints;
};

#endif