#include "Context.hpp"
#include "Scope.hpp"
#include "Parameter.hpp"

namespace xp
{
  Context::Context() :
    params(NULL),
    scope(NULL),
    sweepParameterName(""),
    sweepParameterValue(0),
    phase(kInit)
  {}

  Context::Context(Scope * s) :
    scope(s),
    sweepParameterName(""),
    sweepParameterValue(0),
    phase(kInit)
  {}
/*
  Error Context::getNumericParameter(const std::string& parameterName, double * value)
  {
    if(!parameterName.compare(sweepParameterName))
    {
      *value = sweepParameterValue;
      return kNoErr;
    }

    if(scope->params()->find(parameterName)==scope->params()->end())
      return kNameNotFound;

    *value = (*scope->params())[parameterName].doubleVal;

    return kNoErr;
  }

  Error Context::getStringParameter(const std::string& parameterName, std::string * value)
  {
    if(scope->params()->find(parameterName)==scope->params()->end())
      return kNameNotFound;

    *value = (*scope->params())[parameterName].stringVal;

    return kNoErr;
  }

  Error Context::getBoolParameter(const std::string& parameterName, bool * value)
  {
    if(scope->params()->find(parameterName)==scope->params()->end())
      return kNameNotFound;

    *value = (*scope->params())[parameterName].boolVal;

    return kNoErr;
  }*/


  Error Context::log(const std::string message)
  {
    scope->log(message);

    return kNoErr;
  }

  void Context::sampleValues(std::vector<double> & values)
  {
    scope->sampleValues(values);
  }

  Error Context::openDataFileStream(const std::string& name, std::ofstream& stream)
  {
    std::stringstream filename;
    filename << scope->dataFilename() << name;

    switch(phase)
    {
    case kInit:
      filename << "_init";
      break;
    case kFinish:
      filename << "_finish";
      break;
    case kRunning:
      if(sweepParameterName.length())
	filename << "_" << sweepParameterName << "_" << sweepParameterValue;
      break;
    }

    filename << ".txt";

    stream.open(filename.str().c_str());

    if(!stream)
      return kFileSystemError;

    return kNoErr;
  }



}
