/* Copyright 2016 Pierre Abbat.
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
/* The discrepancy of that part of a sequence that is in a box is the absolute
 * value of the difference between the fraction of the sequence that is in the
 * box (considering both the closed box and the open box) and the volume of
 * the box. The discrepancy of a sequence is the maximum, over all boxes with
 * points in the sequence on all faces (including those where the corners are
 * the same point), of the discrepancy of that part of the sequence that is in
 * the box. The star-discrepancy is the maximum over all boxes with one corner
 * at the origin and points in the sequence on all faces not through the origin.
 * 
 * Computing the star-discrepancy takes O((n**s)ns) time, where n is the length of
 * the sequence and s is the dimension. Computing the discrepancy takes O((n**2s)ns)
 * time. Quadlods therefore computes the star-discrepancy, using random walks
 * from (0,...,0) to (1,...,1) to get an estimate in a reasonable time.
 * 
 * Sequences generated by Quadlods never have a coordinate equal to 0.
 * Coordinates cannot equal 1 either unless due to roundoff error.
 */
#include <cstdlib>
#include <set>
#include <iostream>
#include "discrepancy.h"
using namespace std;

void remove(int n,vector<int> &v)
{
  int i,j;
  for (i=0,j=v.size()-1;i<j;i++)
  {
    while (j>i && v[j]==n)
      j--;
    while (i<j && v[i]!=n)
      i++;
    if (i<j)
      swap(v[i],v[j]);
  }
  v.resize(i);
}

double stardiscrepancy(const vector<vector<double> > &sequence)
{
  int i,j,dbl,n=0;
  vector<set<double> > ordaxes;
  vector<set<double>::iterator> axit;
  vector<double> corner;
  vector<int> dim;
  bool x;
  if (sequence.size())
  {
    ordaxes.resize(sequence[0].size());
    axit.resize(sequence[0].size());
    corner.resize(sequence[0].size());
  }
  for (i=0;i<sequence.size();i++)
    for (j=0;j<sequence[i].size();j++)
      ordaxes[j].insert(sequence[i][j]);
  for (dbl=0;dbl<ordaxes.size();dbl++)
  {
    dim.clear();
    for (i=0;i<ordaxes.size();i++)
    {
      axit[i]=ordaxes[i].begin();
      corner[i]=*axit[i];
      dim.push_back(i);
    }
    dim.push_back(dbl);
    while (dim.size())
    {
      for (i=0;i<corner.size();i++)
      {
	if (i)
	  cout<<',';
	cout<<corner[i];
      }
      cout<<endl;
      x=true;
      while (x)
      {
	n=(n+random()%dim.size())%dim.size();
	axit[dim[n]]++;
	if (axit[dim[n]]==ordaxes[dim[n]].end())
	{
	  remove(dim[n],dim);
	  if (dim.size()==0)
	    x=false;
	}
	else
	{
	  corner[dim[n]]=*axit[dim[n]];
	  x=false;
	}
      }
    }
  }
  return 0;
}
