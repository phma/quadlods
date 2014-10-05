/******************************************************/
/*                                                    */
/* ps.h - PostScript output                           */
/*                                                    */
/******************************************************/

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
