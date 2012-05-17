#ifndef _SCOPE_H_
#define _SCOPE_H_

#include <sys/stat.h>
#include <sys/time.h>
#include <errno.h>
#include <time.h>
#include <vector>
#include <string>
#include <map>

#include <fstream>
#include <iostream>

#include <tclap/CmdLine.h>
#include <libconfig.h++>

#include <mulog/core>

#include "Types.hpp"
#include "Parameter.hpp"

//#define USE_MPI
//#define USE_OMP

#define BUF_SIZE 128 // This should be long enough for all string sending purposes

#ifdef USE_OMP
#include <omp.h>
#endif

#ifdef USE_MPI
#include <mpi.h>
#endif

#define VERSION "0.1"
#define CMD_STRING "Xperience experiment execution framework"

namespace xp
{
  template <typename Type>
  class TypedParameter;

  class Parameter;
  class Context;
  class Experiment;

  typedef enum
  {
    kSweepChunk=0,
    kLogFile,
    kSize,
    kStats,
    kSamples,
    kDone,
    kDie
  } xpMessage;

  class Scope
  {
  public:
    Scope(const std::string & name, const std::string & sha1);
    ~Scope();

    Error parseCommandLine(int argc, char ** argv);
    Error runExperiment(std::string experimentName,
                        bool sweep=false,
                        std::string sweepParameter="",
                        double begin=0.0,
                        double end=0.0,
                        double step=0.0);

    Error extractParameters(std::string setPath);
    Error registerExperiment(const std::string & experimentName,
                             Callback experiment,
                             Callback init,
                             Callback finish);


    template <typename Type>
    Error registerParameter(const std::string & parameterName, const std::string & parameterDescription)
    {
      if(_params.find(parameterName) != _params.end())
        return kNameAlreadyRegistered;

      // Create a new parameter
      TypedParameter<Type> * newParameter=new TypedParameter<Type > (parameterName, parameterDescription);

      // Store parameter
      _params[parameterName]=newParameter;

      // Register with cmd line options
      _cmd.add(newParameter->arg);

      return kNoErr;
    }

    template <typename Type>
    Error setParameter(const std::string &parameterName, Type paramValue)
    {
      if(_params.find(parameterName) != _params.end())
        return kNameNotFound;

      ((TypedParameter<Type>*)_params[parameterName])->val=paramValue;

      return kNoErr;
    }

    template <typename Type>
    Error getParameter(const std::string &parameterName, Type * paramValue)
    {
      if(_params.find(parameterName) != _params.end())
        return kNameNotFound;

      *paramValue=((TypedParameter<Type>*)_params[parameterName])->val;

      return kNoErr;
    }

    void log(const std::string str);
    void sampleValues(std::vector<double>& values);

    std::string dataFilename();
    std::map<std::string, Parameter*>* params();

    double msecs();

    mulog::logger & logger() { return scopeLogger; };

  private:

    void initCommandLine();
    void listExperiments();
    void listParams();

    std::string mpiLog(int rank, int targetRank= -1, int srcRank= -1);
    std::string dateTimeString();
    std::string timeString();

#ifdef USE_MPI
    int receiveResult();
#endif

    void sampleStats(double lastRuntime,
                     double totalRuntime,
                     bool sweep=false,
                     double sweepParam=0.0,
                     int threadNum=0,
                     int procNum=0);

    void error(std::stringstream& string);
    void stats(std::stringstream& string);
    void msg(std::stringstream& string);
    void verb(std::stringstream& string);

    void error(const std::string& string);
    void stats(const std::string& string);
    void msg(const std::string& string);
    void verb(const std::string& string);

    bool fexists(std::string & filename);

    std::string _name;
    std::string _sha1;

    TCLAP::CmdLine _cmd;
    TCLAP::SwitchArg _listParams;
    TCLAP::SwitchArg _listExps;
    TCLAP::ValueArg<std::string> _executeArg;
    TCLAP::ValueArg<std::string> _paramFileArg;
    TCLAP::ValueArg<std::string> _innerSweepParamArg;
    TCLAP::ValueArg<std::string> _outputDirectory;
    TCLAP::ValueArg<std::string> _tag;
    TCLAP::SwitchArg _stdoutStatistics;
    TCLAP::SwitchArg _verbose;

    TCLAP::ValueArg<int> _threadArg;
    TCLAP::ValueArg<int> _threadChunks;
    TCLAP::ValueArg<int> _processChunks;

    std::map<std::string, Experiment> _experiments;
    std::map<std::string, Parameter*> _params;

    std::string _logdir;
    std::string _dataFilename;
    std::ofstream _samplefile;
    std::ofstream _statsfile;

    libconfig::Config _paramFile;

    mulog::logger scopeLogger;
  };

}

#endif /* _SCOPE_H_ */

