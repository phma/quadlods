/******************************************************/
/*                                                    */
/* flowertest.h - draw flower diagrams of sequence    */
/*                                                    */
/******************************************************/
/* Copyright 2018,2019 Pierre Abbat.
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
#include "quadlods.h"
#include "ps.h"

void flowertest(Quadlods &quad,int iters,PostScript &ps,bool disc2d);
void quadplot(PostScript &ps);
