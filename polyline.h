/******************************************************/
/*                                                    */
/* polyline.h - polylines                             */
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
#ifndef POLYLINE_H
#define POLYLINE_H

#include <vector>
#include "xy.h"

class polyline
{
protected:
  std::vector<xy> endpoints;
  std::vector<double> lengths,cumLengths;
public:
  bool isopen();
  int size();
  xy getEndpoint(int i);
  virtual void insert(xy newpoint,int pos=-1);
  virtual void setlengths();
  virtual void open();
  virtual void close();
  double length();
};

#endif
