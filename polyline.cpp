/******************************************************/
/*                                                    */
/* polyline.cpp - polylines                           */
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
#include <iostream>
#include "polyline.h"
#include "manysum.h"
#include "ldecimal.h"
using namespace std;

bool polyline::isopen()
{
  return endpoints.size()>lengths.size();
}

int polyline::size()
{
  return lengths.size();
}

xy polyline::getEndpoint(int i)
{
  i%=endpoints.size();
  if (i<0)
    i+=endpoints.size();
  return endpoints[i];
}

void polyline::insert(xy newpoint,int pos)
/* Inserts newpoint in position pos. E.g. insert(xy(8,5),2) does
 * {(0,0),(1,1),(2,2),(3,3)} -> {(0,0),(1,1),(8,5),(2,2),(3,3)}.
 * If the polyline is open, inserting a point in position 0, -1, or after the last
 * (-1 means after the last) results in adding a line segment.
 * If the polyline is empty (and therefore closed), inserting a point results in
 * adding a line segment from that point to itself.
 * In all other cases, newpoint is inserted between two points and connected to
 * them with line segments.
 */
{
  bool wasopen;
  int i;
  vector<xy>::iterator ptit;
  vector<double>::iterator lenit;
  if (newpoint.isnan())
    cerr<<"Inserting NaN"<<endl;
  wasopen=isopen();
  if (pos<0 || pos>endpoints.size())
    pos=endpoints.size();
  ptit=endpoints.begin()+pos;
  lenit=lengths.begin()+pos;
  endpoints.insert(ptit,newpoint);
  lengths.insert(lenit,0);
  lenit=cumLengths.begin()+pos;
  if (pos<cumLengths.size())
    cumLengths.insert(lenit,cumLengths[pos]);
  else
    cumLengths.insert(lenit,0);
  pos--;
  if (pos<0)
    if (wasopen)
      pos=0;
    else
      pos+=endpoints.size();
  for (i=0;i<2;i++)
  {
    if (pos+1<endpoints.size())
      lengths[pos]=dist(endpoints[pos],endpoints[pos+1]);
    if (pos+1==endpoints.size() && !wasopen)
      lengths[pos]=dist(endpoints[pos],endpoints[0]);
    pos++;
    if (pos>=lengths.size())
      pos=0;
  }
}

/* After inserting, opening, or closing, call setlengths before calling
 * length. If insert called setlengths, manipulation would take
 * too long. So do a lot of inserts, then call setlengths.
 */
void polyline::setlengths()
{
  int i;
  manysum m;
  assert(lengths.size()==cumLengths.size());
  for (i=0;i<lengths.size();i++)
  {
    if (i+1<endpoints.size())
      lengths[i]=dist(endpoints[i],endpoints[i+1]);
    else
      lengths[i]=dist(endpoints[i],endpoints[0]);
    m+=lengths[i];
    cumLengths[i]=m.total();
  }
}

double polyline::length()
{
  if (cumLengths.size())
    return cumLengths.back();
  else
    return 0;
}

void polyline::open()
{
  lengths.resize(endpoints.size()-1);
  cumLengths.resize(endpoints.size()-1);
}

void polyline::close()
{
  lengths.resize(endpoints.size());
  cumLengths.resize(endpoints.size());
  if (lengths.size())
    if (lengths.size()>1)
      cumLengths[lengths.size()-1]=cumLengths[lengths.size()-2]+lengths[lengths.size()-1];
    else
      cumLengths[0]=lengths[0];
}
