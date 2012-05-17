#ifndef _EXPERIMENT_H_
#define _EXPERIMENT_H_

#include "Types.hpp"

namespace xp
{
  class Experiment
  {
  public:
    Experiment();
    Experiment(Callback exp, Callback in, Callback fin);

    Callback experiment;
    Callback init;
    Callback finish;
  };
}

#endif /* _EXPERIMENT_H_ */
