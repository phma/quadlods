/******************************************************/
/*                                                    */
/* random.h - random numbers                          */
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
#include <gmpxx.h>
#include "config.h"

class randm
{
public:
  randm();
  unsigned int uirandom();
  unsigned short usrandom();
  unsigned char ucrandom();
  double expirandom();
  double expsrandom();
  double expcrandom();
  bool brandom();
  mpz_class rangerandom(mpz_class range);
  unsigned int rangerandom(unsigned int range)
  {
    return rangerandom(mpz_class(range)).get_ui();
  }
  bool frandom(mpq_class prob);
  ~randm();
private:
#if defined(_WIN32)
  unsigned int usbuf,ucbuf,usnum,ucnum;
#else
  FILE *randfil;
#endif
  unsigned int bitbuf,bitcnt;
  mpz_class bigacc,bigrange;
};

extern randm rng;
