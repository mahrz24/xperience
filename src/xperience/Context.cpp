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

  mulog::logger Context::openDataLogger(const std::string& name)
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

    mulog::logger l;

    l.add_transformer<mulog::default_transformer,
      mulog::file_device>(mulog::prefix::severity, filename.str());

    return l;
  }



}
