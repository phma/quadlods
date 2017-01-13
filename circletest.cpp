/* Copyright 2017 Pierre Abbat.
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
#include <iostream>
#include <iomanip>
#include <vector>
#include <ctime>
#include <cmath>
#include "circletest.h"

using namespace std;

void circletest(quadlods &quad)
/* Find all pairs of primes (dimensions) which have high discrepancy by
 * counting the pairs of points that fall inside a quadrant of a circle.
 * The number of points in the quadrant should be (π/4)i+O(ln(i)²/i).
 */
{
  int i,j,k,inx,iters=1000000;
  time_t now,then;
  vector<int> incircle;
  vector<double> point;
  double relativeError;
  incircle.resize(quad.size()*(quad.size()-1)/2);
  for (i=0;i<iters;i++)
  {
    point=quad.dgen();
    for (j=0;j<quad.size();j++)
      for (k=0;k<j;k++)
      {
        inx=j*(j-1)/2+k;
        if (point[j]*point[j]+point[k]*point[k]<1)
          incircle[inx]++;
        if (point[j]*point[j]+(1-point[k])*(1-point[k])<1)
          incircle[inx]--;
      }
    now=time(nullptr);
    if (now!=then)
    {
      cout<<rint((double)i/iters*100)<<"% \r";
      cout.flush();
      then=now;
    }
  }
  for (j=0;j<quad.size();j++)
    for (k=0;k<j;k++)
    {
      inx=j*(j-1)/2+k;
      relativeError=(incircle[inx])/log(i)/log(i);
      cout<<j<<' '<<k<<' '<<incircle[inx]<<' '<<relativeError<<endl;
    }
}
