#include <vector>
#include <gmpxx.h>

class quadlods
{
protected:
  std::vector<mpz_class> num,denom,acc;
public:
  void init(int dimensions,double resolution);
  mpz_class getnum(int n)
  {
    return num[n];
  }
  mpz_class getdenom(int n)
  {
    return denom[n];
  }
  mpz_class getacc(int n)
  {
    return acc[n];
  }
  std::vector<mpq_class> readout();
  void setmiddle();
  void advance(mpz_class n);
  std::vector<mpq_class> gen();
  std::vector<double> dreadout();
  std::vector<double> dgen();
};
