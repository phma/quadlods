/******************************************************/
/*                                                    */
/* ps.h - PostScript output                           */
/*                                                    */
/******************************************************/
/* Copyright 2014 Pierre Abbat.
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

extern FILE *psfile;
extern int orientation;
void psprolog();
void startpage();
void endpage();
void dot(double x,double y);
void pstrailer();
void psopen(const char * psfname);
void psclose();
void line2p(double x1,double y1,double x2,double y2);
void widen(double factor);
void setcolor(double r,double g,double b);
void setscale(double minx,double miny,double maxx,double maxy);
