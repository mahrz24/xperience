#include "Experiment.hpp"

namespace xp
{
  Experiment::Experiment()
    {};

  Experiment::Experiment(Callback exp, Callback in, Callback fin) : 
    experiment(exp), 
    init(in), 
    finish(fin) 
  {};
}
