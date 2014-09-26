#include <vector>
#include <gmpxx.h>

class quadlods
{
public:
  std::vector<mpz_class> num,denom,inc;
  void init(int dimensions,double resolution);
  std::vector<mpq_class> gen();
};
