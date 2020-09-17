/******************************************************/
/*                                                    */
/* discrepancy.h - compute discrepancy                */
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
#include <vector>
#include <array>
#include "threads.h"
/* This computes the discrepancy using a genetic algorithm like that invented
 * by Manan Shah. Each individual is a box; its fitness is its discrepancy.
 * In each generation, the least fit boxes die, and the remaining boxes have
 * children, with occasional mutations.
 */
#define sizeMismatch 1

class Box
{
public:
  Box();
  Box(std::vector<double> pnt0,std::vector<double> pnt1);
  Box(Box &mother,Box &father);
  int in(const std::vector<double> &point);
  void countPoints(const std::vector<std::vector<double> > &points);
  double discrepancy();
  void mutate(const std::vector<std::vector<double> > &points);
private:
  std::vector<std::array<double,2> > bounds;
  double volume;
  int pointsIn,pointsBound,pointsTotal;
};

struct BoxCountItem
{
  Box *box;
  const std::vector<std::vector<double> > &points;
};

class BoxCountBlock
{
public:
  BoxCountBlock();
  void load(std::vector<Box> &population,int begin,int end,const std::vector<std::vector<double> > &points);
  BoxCountItem getItem();
  void countFinished();
  bool done();
private:
  std::vector<Box> *pop;
  int b,e,left;
  const std::vector<std::vector<double> > *pts;
  std::mutex mtx;
};

double discrepancy(const std::vector<std::vector<double> > &points);

bool countAnyBlock();
