/******************************************************/
/*                                                    */
/* dotbaton.cpp - line of dots and twirling baton     */
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

#include <iostream>
#include <cmath>
#include "dotbaton.h"
using namespace std;

char batonChars[4]={'-','/','|','\\'};

DotBaton::DotBaton()
{
  lastDots=lastBaton=0;
}

void DotBaton::update(double progress,int activity)
/* Updates the display to show a number of dots proportional to progress,
 * followed by a baton depending on activity. If progress is 0 or 1,
 * the baton is hidden.
 */
{
  int i;
  int dots=lrint(progress*64);
  int baton=activity&3;
  if (lastProgress==0 || lastProgress==1)
    cout<<batonChars[lastBaton];
  if (baton!=lastBaton || dots!=lastDots)
    cout<<"\b \b";
  for (i=dots;i<lastDots;i++)
    cout<<"\b \b";
  for (i=lastDots;i<dots;i++)
    cout<<'.';
  if (baton!=lastBaton || dots!=lastDots)
    cout<<batonChars[lastBaton];
  if (progress==0 || progress==1)
    cout<<"\b \b";
  cout.flush();
  lastDots=dots;
  lastBaton=baton;
  lastProgress=progress;
}
