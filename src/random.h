#ifndef METAOPT_RANDOM_H
#define METAOPT_RANDOM_H
#include<random>
namespace MetaOpt{

class RandomClass{
public:
  RandomClass(){
    std::random_device rd;
    engine_.seed(rd());
  };
  std::mt19937_64 engine_;
  // gen.seed(rd());
  /*return a random double in [0,1]*/
  template <typename Real>
  Real gen() { return (Real)engine_() / engine_.max(); };
  /*return a random double in [a,b], numbers can be gened*/
  template <typename T>
  T gen(T a, T b){
    return gen(a, b);
  };
  int    gen(int a, int b) { return engine_() % (b - a + 1) + a; };
  float  gen(float a, float b) { return (b - a) * engine_() / engine_.max() + a; };
  double gen(double a, double b) { return (b - a) * engine_() / engine_.max() + a; };

  template<typename Iter>
  void shuffle( Iter first, Iter last ) {
      std::shuffle( first, last, engine_ );
  };
  template<typename Container>
  void shuffle( Container& container ) {
      shuffle( std::begin( container ), std::end( container ) );
  };
  
};

static RandomClass Random;
}

#endif // METAOPT_RANDOM_H