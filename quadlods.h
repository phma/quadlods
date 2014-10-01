#include <vector>
#include <gmpxx.h>

class quadlods
{
public:
  std::vector<mpz_class> num,denom,acc;
  void init(int dimensions,double resolution);
  std::vector<mpq_class> readout();
  void advance(mpz_class n);
  std::vector<mpq_class> gen();
  std::vector<double> dreadout();
  std::vector<double> dgen();
};
