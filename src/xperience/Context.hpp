#ifndef _CONTEXT_H_
#define _CONTEXT_H_

#include <map>
#include <vector>

#include "Types.hpp"
#include "Scope.hpp"

namespace xp
{
  class Parameter;
  template <typename Type>
  class TypedParameter;

  template <typename Type>
  struct SweepConverter
  {
    static Type toSweep(double x) { return Type(); };
  };

  template <>
  struct SweepConverter<double>
  {
    static double toSweep(double x) { return x; };
  };

  class Context
  {
  public:
    Context();
    Context(Scope * s);

    template<typename Type>
    Error getParameter(const std::string& parameterName, Type * value)
    {
      if(!parameterName.compare(sweepParameterName))
      {
        *value = SweepConverter<Type>::toSweep(sweepParameterValue);
        return kNoErr;
      }

      if(scope->params()->find(parameterName) == scope->params()->end())
        return kNameNotFound;

      *value=((TypedParameter<Type>*)(*scope->params())[parameterName])->val;

      return kNoErr;
    }

    Error log(const std::string message);
    void sampleValues(std::vector<double> & values);

    Error openDataFileStream(const std::string& name, std::ofstream& stream);

    template<typename Type>
    const std::string parameterName(const std::string& name, Type value)
    {
      std::stringstream filename;
      filename << "_" << name << "_" << value;

      return filename.str();
    }

    Scope * scope;
    std::map<std::string, Parameter> * params;
    std::string sweepParameterName;
    double sweepParameterValue;
    Phase phase;
    void * experimentData;
  };


}

#endif /* _CONTEXT_H_ */
