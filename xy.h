/******************************************************/
/*                                                    */
/* xy.h - 2D points                                   */
/*                                                    */
/******************************************************/
/* Copyright 2018,2021 Pierre Abbat.
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
#ifndef XY_H
#define XY_H
#include <cmath>

class xy
{
public:
  xy(double e,double n);
  xy();
  double getx() const;
  double gety() const;
  double length() const;
  bool isfinite() const;
  bool isnan() const;
  friend xy operator+(const xy &l,const xy &r);
  friend xy operator+=(xy &l,const xy &r);
  friend xy operator-=(xy &l,const xy &r);
  friend xy operator-(const xy &l,const xy &r);
  friend xy operator-(const xy &r);
  friend xy operator*(const xy &l,double r);
  friend xy operator*(double l,const xy &r);
  friend xy operator/(const xy &l,double r);
  friend xy operator/=(xy &l,double r);
  friend bool operator!=(const xy &l,const xy &r);
  friend bool operator==(const xy &l,const xy &r);
  friend xy turn90(xy a);
  friend xy turn(xy a,int angle);
  friend double dist(xy a,xy b);
protected:
  double x,y;
};

#endif
