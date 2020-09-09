/******************************************************/
/*                                                    */
/* random.cpp - random numbers                        */
/*                                                    */
/******************************************************/
/* Copyright 2018,2020 Pierre Abbat.
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

#define _CRT_RAND_S
#include <cstdlib>
#include <cstdio>
#include <cassert>
#include <cmath>
#include <iostream>
#include "random.h"
using namespace std;

// Cygwin does not define _WIN32, but does have /dev/urandom
#if defined(_WIN32)
randm::randm()
{
  ucnum=usnum=bitcnt=0;
  bigrange=1;
}

randm::~randm()
{
  ucnum=usnum=-1;
}

unsigned int randm::uirandom()
{
  unsigned int n;
  rand_s(&n);
  return n;
}

unsigned short randm::usrandom()
{
  unsigned short n;
  if (!usnum)
    rand_s(&usbuf);
  n=(usbuf>>usnum)&0xffff;
  usnum=(usnum+16)&31;
  return n;
}

unsigned char randm::ucrandom()
{
  unsigned char n;
  if (!ucnum)
    rand_s(&ucbuf);
  n=(ucbuf>>ucnum)&0xff;
  ucnum=(ucnum+8)&31;
  return n;
}
#else
randm::randm()
{
  randfil=fopen("/dev/urandom","rb");
  bigrange=1;
}

randm::~randm()
{
  fclose(randfil);
}

unsigned int randm::uirandom()
{
  unsigned int n;
  fread(&n,1,4,randfil);
  return n;
}

unsigned short randm::usrandom()
{
  unsigned short n;
  fread(&n,1,2,randfil);
  return n;
}

unsigned char randm::ucrandom()
{
  unsigned char n;
  fread(&n,1,1,randfil);
  return n;
}
#endif

double randm::expirandom()
{
  return -log((uirandom()+0.5)/4294967296.);
}

double randm::expsrandom()
{
  return -log((usrandom()+0.5)/65536.);
}

double randm::expcrandom()
{
  return -log((ucrandom()+0.5)/256.);
}

bool randm::brandom()
{
  bool ret;
  if (bitcnt==0)
  {
    bitbuf=uirandom();
    bitcnt=32;
  }
  ret=bitbuf&1;
  bitbuf>>=1;
  bitcnt--;
  return ret;
}

mpz_class randm::rangerandom(mpz_class range)
{
  mpz_class ret;
  assert(range>0);
  while (bigrange<16777216*range || bigacc/range==bigrange/range)
  {
    bigacc=(bigacc<<32)+uirandom();
    bigrange<<=32;
  }
  ret=bigacc%range;
  bigrange/=range;
  bigacc/=range;
  return ret;
}

bool randm::frandom(mpq_class prob)
{
  bool ret;
  mpz_class upper0,lower1;
  prob.canonicalize();
  assert(prob>=0 && prob<=1);
  while (true)
  {
    upper0=bigrange*prob.get_num()/prob.get_den();
    lower1=(bigrange*prob.get_num()+prob.get_den()-1)/prob.get_den();
    if (bigacc<upper0 || bigacc>=lower1)
      break;
    bigacc=(bigacc<<32)+uirandom();
    bigrange<<=32;
  }
  ret=bigacc<upper0;
  if (ret)
    bigrange=upper0;
  else
  {
    bigacc-=lower1;
    bigrange-=lower1;
  }
  return ret;
}

randm rng;
