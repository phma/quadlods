/******************************************************/
/*                                                    */
/* interact.cpp - interactive mode                    */
/*                                                    */
/******************************************************/
/* Copyright 2019 Pierre Abbat.
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
#include <cassert>
#include "config.h"
#include "interact.h"
using namespace std;

int commandInt(string &command)
/* Removes the first word of command and returns it as an int.
 * The commands are four letters each.
 */
{
  int i,ret=0,ch;
  for (i=0;i<command.length();i++)
  {
    ch=command[i];
    if (isspace(ch))
      break;
    ch=toupper(ch);
    if (i<4)
      ret+=(ch&0xff)<<(8*(3-i));
  }
  command.erase(0,i);
  return ret;
}

void reply(int code,bool done,string text)
{
  assert(code>=100 && code<1000);
  cout<<code<<(done?' ':'-')<<text<<endl;
}

/* Commands for interactive mode, which can be used as a server:
 * init n s res: Initialize generator #n with s dimensions and resolution res.
 * form n dec/hex/flo/rat: Set format to decimal/hexadecimal/floating point/rational.
 * gene n i: Generate i points from generator n.
 * seed n: Seed generator n with random numbers.
 */
void interact()
{
  bool cont=true;
  string command;
  int opcode;
  reply(220,true,string("Quadlods version ")+VERSION+" ready");
  while (cont)
  {
    getline(cin,command);
    opcode=commandInt(command);
    switch (opcode)
    {
      case 0x45584954:
      case 0x51554954:
	cont=false;
	reply(221,true,"Quadlods exiting");
	break;
      default:
	reply(400,true,"Invalid command");
    }
  }
}
