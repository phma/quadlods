#include <vector>
#include <gmpxx.h>

#define QL_JUMBLE_NONE 0
#define QL_JUMBLE_THIRD 1
#define QL_JUMBLE_THUEMORSE 2
#define QL_JUMBLE_GRAY 3
/* The jumbletype controls how to jumble the bits of the accumulator when
 * reading the generator. If acc is 0xc0de and denom is 0x10000, it returns:
 * QL_JUMBLE_NONE:      1100000011011110
 * QL_JUMBLE_THIRD:     1001010110001011 (xor with 5555, 1/3=0.5555... hex)
 * QL_JUMBLE_THUEMORSE: 0101011010110111 (xor with the Thue-Morse word)
 * QL_JUMBLE_GRAY:      1000000010010100 (inverse Gray code).
 * Of these, Gray code is the best, as it leaves no tendency to slope
 * one way or the other. It is therefore the default. However, Gray code
 * takes about 5/3 times as much time as the others.
 */

class quadlods
{
protected:
  std::vector<mpz_class> num,denom,acc;
  int jumbletype;
public:
  void init(int dimensions,double resolution,int j=QL_JUMBLE_GRAY);
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
  void setjumble(int j)
  {
    jumbletype=j;
  }
  void advance(mpz_class n);
  std::vector<mpq_class> gen();
  std::vector<double> dreadout();
  std::vector<double> dgen();
};
