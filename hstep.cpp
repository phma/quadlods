/******************************************************/
/*                                                    */
/* hstep.cpp - half step sequences                    */
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

#include <cmath>
#include "hstep.h"
using namespace std;

double halfstep[]={1,1.05946309435929526455,1.12246204830937298142};
double minorthird[]={1,1.18920711500272106671,1.41421356237309504878,1.68179283050742908604};

double hstep(int i)
{
  int octave,note;
  octave=i/12;
  if (i<octave*12)
    --octave;
  note=i-12*octave;
  return ldexp(minorthird[note/3]*halfstep[note%3],octave);
}

set<int> hsteps(int start,int iters)
/* Returns a set of numbers which are half steps rounded to integers.
 * The largest possible number is 2026954652.
 */
{
  set<int> ret;
  int i;
  unsigned n;
  for (i=-12;i<372;i++)
  {
    n=rint(iters/hstep(i));
    if (n>=start && n<=iters)
      ret.insert(n);
  }
  return ret;
}

void xticks(int start,int iters,PostScript &ps)
// Draw tickmarks at every power of 10 on the log-scaled x axis.
{
  int i;
  double x;
  char buf[24];
  for (i=1;i<=iters && i%10!=8;i*=10)
  {
    sprintf(buf,"%g",(double)i);
    x=log(i+1-start)/log(iters)*3;
    ps.write(x,-1.1,buf);
    ps.startline();
    ps.lineto(x,-1);
    ps.lineto(x,-1.1);
    ps.endline();
  }
}
