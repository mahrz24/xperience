#include "Scope.hpp"
#include "Experiment.hpp"
#include "Context.hpp"
#include "Parameter.hpp"


namespace xp
{

  Scope::Scope(const std::string & name, const std::string & sha1) :
  _name(name),
  _sha1(sha1),
  _cmd(CMD_STRING, ' ', VERSION),
  _listParams("l", "list-params",
              "List experiment parameters", false),
  _listExps("e", "list-experiments",
            "List experiments", false),
  _executeArg("x", "execute",
              "Experiment to execute", true, "default", "experiment"),
  _paramFileArg("p", "parameters",
                "Filename of parameter file", false, "params.cfg", "filename:paramset"),
  _innerSweepParamArg("i", "inner",
                      "Parameter name of the inner sweep loop", false, "beta", "param"),
  _outputDirectory("o", "outdir",
                   "Directory where experiment outputs are stored", false, ".",
                   "directory name"),
  _tag("g", "tag",
                   "Name of the log directory, if unspecified timestamp is used", false, ".",
                   "directory name"),
  _stdoutStatistics("s", "stdout",
                    "Output statistics to STDOUT", false),
  _verbose("v", "verbose",
           "Output everything to STDOUT", false),
  _threadArg("t", "threads",
             "Number of OpenMP threads to use for sweep", false, 0, "num of threads"),
  _threadChunks("c", "threadChunks",
                "Number of chunks a process range is divided into (OMP)",
                false, 10, "num of chunks"),
  _processChunks("r", "processChunks",
                 "Number of process chunks a sweep range is divided into (MPI)",
                 false, 10, "num of chunks")
  {
    initCommandLine();
  }

  Scope::~Scope()
  {
    for(std::map<std::string, Parameter*>::iterator cur=_params.begin(); cur != _params.end(); cur++)
    {
      delete cur->second;
    }
  }

  Error Scope::parseCommandLine(int argc, char ** argv)
  {
    Error err=kNoErr;

    _cmd.parse(argc, argv);

    if(_verbose.isSet())
      scopeLogger.set_min_severity(mulog::verbose2);
    else
      scopeLogger.set_min_severity(mulog::info);

    scopeLogger.add_transformer<mulog::default_transformer,
      mulog::console_device>(mulog::prefix::severity);

#if defined USE_MPI
    int myRank, numProcs;

    MPI_Init(&argc, &argv); // Check whether this is ok ? Should we maybe put a NULL in here?

    // Get current rank
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    // Get number of running processes
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

    std::cout << "[MPI #" << myRank << "] MPI Process started. Number of processes: " << numProcs << std::endl;
#endif

#if defined USE_OMP
    if(_threadArg.isSet() && _threadArg.getValue() > 0)
      omp_set_num_threads(_threadArg.getValue());
#endif

    if(_listExps.isSet())
    {
      listExperiments();
    }

    if(_listParams.isSet())
    {
      listParams();
    }

    if(_executeArg.isSet())
    {
      // First check whether such an experiment is registered
      if(_experiments.find(_executeArg.getValue()) == _experiments.end())
      {
        error("Experiment not found");
        return kExperimentNotFound;
      }

#ifdef USE_MPI
      if(myRank == 0)
      {
#endif
        // Master process, just do the usual (except for the runexperiment, where only chunks are dispatched
        verb("Checking output directory");

        // Create output dir if neccessary and create logging dir
        if(_outputDirectory.isSet())
        {
          mode_t process_mask=umask(0);
          int result_code=mkdir(_outputDirectory.getValue().c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
          umask(process_mask);

          if(result_code && errno != EEXIST)
          {
            error("Could not create directory");
            return kFileSystemError;
          }

          _logdir=_outputDirectory.getValue();
        }
        else
          _logdir="";

        if(_tag.isSet())
        {
          _logdir+="/" + _tag.getValue();
        }
        else
        {
          _logdir+="/" + dateTimeString();
        }

        std::string odir = _logdir;
        int i=1;
        int result_code = 0;
        do
        {
          mode_t process_mask=umask(0);
          result_code=mkdir(_logdir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
          umask(process_mask);

          if(result_code)
          {
            if(errno == EEXIST)
              msg("A logging directory with the same tag/timestamp exists");
            else
            {
              error("Could not create logging directory");
              return kFileSystemError;
            }
            std::stringstream s;
            s << i;
            i++;
            _logdir = odir + "-" + s.str();
          }
        } while(result_code);
        // Logging directory exists now

        // Create the log files
        std::string filename(_logdir + "/log.txt");

        scopeLogger.add_transformer<mulog::default_transformer,
          mulog::file_device>(mulog::prefix::extended, filename);

        filename=_logdir + "/stats.txt";
        _statsfile.open(filename.c_str());
        filename=_logdir + "/samples.txt";
        _samplefile.open(filename.c_str());

        _dataFilename=_logdir + "/data_";

        // Output running parameters & hash to log
        std::stringstream args;
        for(int i=0; i < argc; i++)
          args << argv[i] << " ";
        msg(args.str());
        msg("SHA1#" + _sha1);

#ifdef USE_MPI
      }
#endif


      // Parameter file handling
      std::string paramFile(_executeArg.getValue() + ".cfg");
      std::string paramSet("");
      std::string setPath("");

      if(_paramFileArg.isSet())
      {
        verb("Loading parameter configuration");
        // Load configuration either specified, <experiment>.cfg or no values given

        paramFile=_paramFileArg.getValue();
        int setStart=0;

        if((setStart=paramFile.find_first_of(':')) != std::string::npos)
        {
          paramSet=paramFile.substr(setStart + 1, paramFile.size() - setStart - 1);
          paramFile=paramFile.substr(0, setStart);
        }

        if(!paramFile.length())
          paramFile=(_executeArg.getValue() + ".cfg");
      }
      else
      {
        if(!fexists(paramFile))
          paramFile=std::string("");
      }

      if(paramFile.length())
      {
        if(paramSet.length())
          verb(std::string("Loading set ") + paramSet + std::string(" of parameter file ") + paramFile);
        else
          verb(std::string("Loading parameter file ") + paramFile);

        // Load parameters
        FILE * paramFD=fopen(paramFile.c_str(), "r");

        if(!paramFD)
        {
          error("Could not find parameter file");
          return kFileNotFound;
        }

        try
        {
          _paramFile.read(paramFD);
        } catch(libconfig::ParseException &e) // catch any exceptions
        {
          fclose(paramFD);
          error("Parameter file invalid: " + std::string(e.what()));
          return kParamParseError;
        }

        fclose(paramFD);

        int rootLength=_paramFile.getRoot().getLength();

        if(rootLength == 0)
          msg("WARNING: Parameter file empty");
        else
        {

          // Check if parameters are in correct format
          if(paramSet.length())
          {
            for(int i=0; i < rootLength; i++)
              if(_paramFile.getRoot()[i].getType() != libconfig::Setting::TypeGroup)
              {
                error("Parameter file in wrong format, no parameter set found");
                return kParamFileFormatError;
              }

            try
            {
              setPath=_paramFile.lookup(paramSet).getPath();
            } catch(std::exception &e)
            {
              error("Parameter set not found: " + std::string(e.what()));
              return kSetNotFound;
            }

          }
          else
          {
            setPath=_paramFile.getRoot().getPath();
          }

          // Found matching parameter set now extract parameters
          Error paramErr;
          if((paramErr=extractParameters(setPath + ".")) != kNoErr)
          {
            return paramErr;
          }
        }

#ifdef USE_MPI
        if(myRank == 0)
        {
#endif
          // Copy the paramfile to log directory
          std::ifstream f1(paramFile.c_str(), std::fstream::binary);
          std::ofstream f2((_logdir + "/params.cfg").c_str(), std::fstream::trunc | std::fstream::binary);
          f2 << f1.rdbuf();

#ifdef USE_MPI
        }
#endif


      }

      // Overwrite parameters when command line parameter is given
      for(std::map<std::string, Parameter*>::iterator cur=_params.begin(); cur != _params.end(); cur++)
      {
        if(cur->second->isSet())
          cur->second->updateInternalValue();
      }

      if(_innerSweepParamArg.isSet())
      {
        // Do a sweep
        double begin, end, step;

        if(!_paramFile.lookupValue(setPath + "." + _innerSweepParamArg.getValue() + ".begin", begin) ||
           !_paramFile.lookupValue(setPath + "." + _innerSweepParamArg.getValue() + ".end", end) ||
           !_paramFile.lookupValue(setPath + "." + _innerSweepParamArg.getValue() + ".step", step))
        {
          error("Sweep range not specified");
          return kSweepRangeNotSpecified;
        }

#ifndef USE_MPI
        err=runExperiment(_executeArg.getValue(), true, _innerSweepParamArg.getValue(), begin, end, step);
#else
        // MPI stuff starts here
        char strBuf[BUF_SIZE];
        std::stringstream l;
        MPI_Status status;
        double chunkRange[3];

        if(myRank == 0)
        {
          // Determine number of steps
          int size=(end - begin) / step;
          int mul=1;

          if(end < begin)
          {
            size*= -1;
            mul= -1;
          }

          int chunk=size / (_processChunks.getValue()) + 1;
          double cur=begin;
          double time=msecs();

          // Output the number of chunks and steps
          l << mpiLog(myRank) << "Chunksize: " << chunk << " Steps: " << size;
          log(l.str());
          l.str("");

          // Dispatch a job to every process running
          for(int rank=1; rank < numProcs; rank++)
          {
            chunkRange[0]=cur;
            chunkRange[1]=cur + mul * chunk*step;
            chunkRange[2]=step;
            cur=chunkRange[1];

            // First send the log directory to every child process
            strcpy(strBuf, _logdir.c_str());
            MPI_Send(strBuf, BUF_SIZE, MPI_CHAR, rank, kLogFile, MPI_COMM_WORLD);

            l << mpiLog(myRank, rank) << "Sending chunk [" << chunkRange[0] << " .. " << chunkRange[1] << "]";
            log(l.str());
            l.str("");

            // Send the actual chunk
            MPI_Send(chunkRange, 3, MPI_DOUBLE, rank, kSweepChunk, MPI_COMM_WORLD);
          }

          bool run=true;

          // Dispatch new chunks until done
          while(run)
          {
            chunkRange[0]=cur;
            chunkRange[1]=cur + mul * chunk*step;
            chunkRange[2]=step;
            cur=chunkRange[1];

            int src=receiveResult();

            l << mpiLog(myRank, -1, src) << "Chunk done";
            log(l.str());
            l.str("");

            // Send the slave a new work unit
            l << mpiLog(myRank, src) << "Sending chunk [" << chunkRange[0] << " .. " << chunkRange[1] << "]";
            log(l.str());
            l.str("");

            // If the next chunk is the last chunk stop sending after that
            if((end < begin && chunkRange[1] < end) ||
               (end >= begin && chunkRange[1] > end))
            {
              run=false;
              chunkRange[1]=end;
            }

            // Send another chunk out
            MPI_Send(chunkRange, 3, MPI_DOUBLE, src, kSweepChunk, MPI_COMM_WORLD);

          }

          for(int rank=1; rank < numProcs; rank++)
          {
            chunkRange[0]=cur;
            chunkRange[1]=cur + mul * chunk*step;
            chunkRange[2]=step;

            cur=chunkRange[1];

            l << mpiLog(myRank) << "Receving remaining results";
            log(l.str());
            l.str("");

            int src=receiveResult();

            l << mpiLog(myRank, -1, src) << "Chunk done";
            log(l.str());
            l.str("");

          }


          for(int rank=1; rank < numProcs; rank++)
          {
            l << mpiLog(myRank, rank) << "Sending die signal";
            log(l.str());
            l.str("");

            MPI_Send(0, 0, MPI_DOUBLE, rank, kDie, MPI_COMM_WORLD);
          }


          // Output the number of chunks and steps
          l << mpiLog(myRank) << "Total time: " << (msecs() - time) << size;
          log(l.str());
          l.str("");

        }
        else
        {
          // Get the logdir and setup logfiles
          MPI_Recv(strBuf, BUF_SIZE, MPI_CHAR, 0, kLogFile, MPI_COMM_WORLD, &status);

          _logdir=std::string(strBuf);

          std::stringstream filename;
          filename << _logdir << "/log-" << myRank << ".txt";
          _logfile.open(filename.str().c_str());

          std::stringstream dfilename;
          dfilename << _logdir << "/data_";
          _dataFilename=dfilename.str();

          // Now we can log messages

          bool running=true;

          // Worker drone
          while(running)
          {

            MPI_Recv(chunkRange, 3, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

            if(status.MPI_TAG == kDie)
            {
              l << mpiLog(myRank, -1, status.MPI_SOURCE) << "Dying";
              log(l.str());
              l.str("");
              running=false;
            }
            else
            {
              l << mpiLog(myRank, -1, status.MPI_SOURCE) << "Received chunk [" << chunkRange[0] << " .. " << chunkRange[1] << "]";
              log(l.str());
              l.str("");

              err=runExperiment(_executeArg.getValue(), true, _innerSweepParamArg.getValue(), chunkRange[0], chunkRange[1], chunkRange[2]);

              if(err != kNoErr)
              {
                l << mpiLog(myRank, -1, status.MPI_SOURCE) << "Chunk error";
                log(l.str());
                l.str("");
              }

              // Send the error back
              MPI_Send(&err, 1, MPI_INT, 0, kDone, MPI_COMM_WORLD);
            }

          }
        }
#endif
      }
      else
      {
        // Run experiment
        err=runExperiment(_executeArg.getValue());
      }

#ifdef USE_MPI
      if(myRank == 0)
      {
#endif
        _statsfile.close();

#ifdef USE_MPI
      }
#endif

    }

#ifdef USE_MPI
    MPI_Finalize();
#endif

    return err;
  }

  Error Scope::runExperiment(std::string experimentName,
                             bool sweep,
                             std::string sweepParameter,
                             double begin,
                             double end,
                             double step)
  {
    std::map<std::string, Experiment>::iterator exp;

    if((exp=_experiments.find(experimentName)) == _experiments.end())
    {
      error("Experiment not found");
      return kExperimentNotFound;
    }

    Context ctx(this);

    if(sweep)
    {
      ctx.sweepParameterName=sweepParameter;
      ctx.sweepParameterValue=begin;
    }

    if((*exp).second.init && (*exp).second.finish)
    {
      Error initError=(*exp).second.init(ctx);

      if(initError != kNoErr)
      {
        error("Error while initializing experiment");
        return initError;
      }
    }

    verb("Experiment initialized");

    ctx.phase=kRunning;

    double totalTime=msecs();
    double runningTime;

    if(sweep)
    {
      double sweepMul=1.0;
      if(begin > end)
      {
        begin= -begin;
        end= -end;
        sweepMul= -1.0;
      }


      int proc=0;

#ifdef USE_MPI
      // Get current rank
      MPI_Comm_rank(MPI_COMM_WORLD, &proc);
#endif

      int max=(end - begin) / step;
      int sweepInt;
      // Chunksize fixed at the moment
      int chunk=max / 10;

      Error expErr=kNoErr;

      int thread=0;

#pragma omp parallel for private(thread, sweepInt, runningTime) firstprivate(ctx) schedule(dynamic,chunk)
      for(sweepInt=0; sweepInt < max; sweepInt++)
      {
        double sweepArg=begin + (sweepInt * step);

        ctx.sweepParameterValue=sweepArg*sweepMul;

        runningTime=msecs();
        Error experimentError=(*exp).second.experiment(ctx);
        runningTime=msecs() - runningTime;

#ifdef USE_OMP
        thread=omp_get_thread_num();
#endif

        sampleStats(runningTime, msecs() - totalTime, true, ctx.sweepParameterValue, thread, proc);

        if(experimentError != kNoErr)
        {
          error("Error while executing experiment");
          expErr=experimentError;
        }
      }

      if(expErr != kNoErr)
        return expErr;

    }
    else
    {
      runningTime=msecs();
      Error experimentError=(*exp).second.experiment(ctx);
      runningTime=msecs() - runningTime;

      if(experimentError != kNoErr)
      {
        error("Error while executing experiment");
        return experimentError;
      }

      totalTime=msecs() - totalTime;
      sampleStats(runningTime, totalTime);

    }

    ctx.phase=kFinish;

    verb("Experiment finishing");
    if((*exp).second.init && (*exp).second.finish)
    {
      return (*exp).second.finish(ctx);
    }

    return kNoErr;
  }

  Error Scope::extractParameters(std::string setPath)
  {
    for(std::map<std::string, Parameter*>::iterator cur=_params.begin();
        cur != _params.end(); cur++)
    {
      switch(cur->second->type)
      {
        case kDouble:
          double dVal;
          if(!_paramFile.lookupValue(setPath + cur->first, dVal))
            if(!_paramFile.lookupValue(setPath + cur->first + ".value", dVal))
              if(!_paramFile.lookupValue(setPath + cur->first + ".begin", dVal))
              {
                int iVal;
                if(!_paramFile.lookupValue(setPath + cur->first, iVal))
                  if(!_paramFile.lookupValue(setPath + cur->first + ".value", iVal))
                    if(!_paramFile.lookupValue(setPath + cur->first + ".begin", iVal))
                      msg("WARNING: Parameter " + cur->first +
                          " not specified in parameter file");

                ((TypedParameter<double>*)cur->second)->val= (double)iVal;
                break;
              }

          ((TypedParameter<double>*)cur->second)->val = dVal;
          break;
        case kInt:
          int iVal;
          if(!_paramFile.lookupValue(setPath + cur->first, iVal))
            if(!_paramFile.lookupValue(setPath + cur->first + ".value", iVal))
              if(!_paramFile.lookupValue(setPath + cur->first + ".begin", iVal))
                msg("WARNING: Parameter " + cur->first +
                    " not specified in parameter file");

          ((TypedParameter<int>*)cur->second)->val= iVal;

          break;
        case kBool:
          bool bVal;
          if(!_paramFile.lookupValue(setPath + cur->first, bVal))
            msg("WARNING: Parameter " + cur->first +
                " not specified in parameter file");

          ((TypedParameter<bool>*)cur->second)->val=bVal;
          break;
        case kUnknown:
          msg("WARNING: Parameter " + cur->first +
                " of unknown type");
          break;
        case kString:
          std::string sVal;
          if(!_paramFile.lookupValue(setPath + cur->first, sVal))
            msg("WARNING: Parameter " + cur->first +
                " not specified in parameter file");

          ((TypedParameter<std::string>*)cur->second)->val=sVal;
          break;
      }
    }

    return kNoErr;
  }

  Error Scope::registerExperiment(const std::string & experimentName,
                                  Callback experiment,
                                  Callback init,
                                  Callback finish)
  {
    if(_experiments.find(std::string(experimentName)) != _experiments.end())
      return kNameAlreadyRegistered;

    Experiment newExperiment(experiment, init, finish);

    _experiments[experimentName]=newExperiment;

    return kNoErr;
  }

  void Scope::log(const std::string str)
  {
    msg(str);
  }

  std::string Scope::dataFilename()
  {
    return _dataFilename;
  }

  void Scope::sampleValues(std::vector<double>& values)
  {
#pragma omp critical (sample)
    {
#ifdef USE_MPI
      int myRank;
      // Get current rank
      MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
      if(myRank == 0)
      {
#endif
        for(std::vector<double>::iterator value=values.begin(); value != values.end(); value++)
          _samplefile << *value << " ";
        _samplefile << std::endl;
#ifdef USE_MPI
      }
      else
      {
        double* valueArray;
        int size=values.size();
        valueArray=new double [size];
        copy(values.begin(), values.end(), valueArray);

        MPI_Send(&size, 1, MPI_INT, 0, kSize, MPI_COMM_WORLD);
        MPI_Send(valueArray, size, MPI_DOUBLE, 0, kSamples, MPI_COMM_WORLD);

        delete [] valueArray;
      }
#endif
    }
  }

  std::map<std::string, Parameter*>* Scope::params()
  {
    return &_params;
  }

  void Scope::initCommandLine()
  {
    std::vector<TCLAP::Arg*> xorlist;
    xorlist.push_back(&_listParams);
    xorlist.push_back(&_listExps);
    xorlist.push_back(&_executeArg);

    _cmd.xorAdd(xorlist);
    _cmd.add(_paramFileArg);

    _cmd.add(_innerSweepParamArg);
    _cmd.add(_outputDirectory);
    _cmd.add(_stdoutStatistics);
    _cmd.add(_verbose);

    _cmd.add(_tag);

#ifdef USE_OMP
    _cmd.add(_threadArg);
    _cmd.add(_threadChunks);
#endif

#ifdef USE_MPI
    _cmd.add(_processChunks);
#endif

  }

  void Scope::listExperiments()
  {
    std::cout << "List of experiments:" << std::endl;
    for(std::map<std::string, Experiment>::iterator cur=_experiments.begin(); cur != _experiments.end(); cur++)
    {
      std::cout << cur->first << std::endl;
    }
  }

  void Scope::listParams()
  {
    std::cout << "List of parameters:" << std::endl;
    for(std::map<std::string, Parameter*>::iterator cur=_params.begin(); cur != _params.end(); cur++)
    {
      std::cout << cur->first << std::endl;
    }
  }

  double Scope::msecs()
  {
    struct timeval tv;
    struct timezone tz;
    gettimeofday(&tv, &tz);

    return tv.tv_sec * 1000 + (tv.tv_usec) / 1000.0;
  }

  std::string Scope::mpiLog(int rank, int targetRank, int srcRank)
  {
    std::stringstream log;
    log << "[MPI #" << rank;
    if(targetRank != -1)
      log << "=>#" << targetRank;
    if(srcRank != -1)
      log << "<=#" << srcRank;
    log << "] ";
    return log.str();
  }

  std::string Scope::dateTimeString()
  {
    time_t seconds;
    struct tm * now;
    char timeStr[48];
    seconds=time(NULL);
    now=localtime(&seconds);

    strftime(timeStr, 48, "%Y-%m-%d-%H-%M-%S", now);

    return std::string(timeStr);
  }

  std::string Scope::timeString()
  {
    time_t seconds;
    struct tm * now;
    char timeStr[48];
    seconds=time(NULL);
    now=localtime(&seconds);

    strftime(timeStr, 48, "%H-%M-%S", now);

    return std::string(timeStr);
  }

#ifdef USE_MPI

  int Scope::receiveResult()
  {
    int result;
    MPI_Status status;

    do
    {
      MPI_Recv(&result, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

      if(status.MPI_TAG == kSize)
      {
        double * dVals=new double[result];
        MPI_Status is;
        MPI_Recv(dVals, result, MPI_DOUBLE, status.MPI_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &is);

        if(is.MPI_TAG == kStats)
        {
          sampleStats(dVals[0], dVals[1], true, dVals[2], dVals[3], dVals[4]);
        }
        else if(is.MPI_TAG == kSamples)
        {
          for(int i=0; i < result; i++)
            _samplefile << dVals[i] << " ";
          _samplefile << std::endl;
        }

        delete [] dVals;
      }
    }
    while(status.MPI_TAG != kDone);

    return status.MPI_SOURCE;
  }
#endif

  void Scope::sampleStats(double lastRuntime, double totalRuntime, bool sweep,
                          double sweepParam, int threadNum, int procNum)
  {
    std::stringstream statsStr("");

#pragma omp critical (stats)
    {
      int numProcs=1;
      int numThreads=1;

#ifdef USE_OMP
      numThreads=omp_get_num_threads();
#endif

#ifdef USE_MPI
      int myRank;

      // Get number of running processes
      MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
      // Get current rank
      MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
      if(myRank == 0)
      {
#endif
        if(sweep)
        {
          _statsfile << sweepParam << " ";
          statsStr << "Sweep param: " << sweepParam << "; ";
        }

        _statsfile << lastRuntime << " " << totalRuntime << " " << threadNum << " " << procNum << std::endl;

        statsStr << "Experiment runtime: " << lastRuntime << "ms on thread " << threadNum << "/" <<
            numThreads << "; Total runtime: " << totalRuntime << "ms on process " << procNum << "/" << numProcs;


        stats(statsStr.str());
#ifdef USE_MPI
      }
      else
      {
        int size=5;
        double statsArray[5];

        statsArray[0]=lastRuntime;
        statsArray[1]=totalRuntime;
        statsArray[2]=sweepParam;
        statsArray[3]=threadNum;
        statsArray[4]=procNum;

        MPI_Send(&size,
                 1,
                 MPI_INT,
                 0,
                 kSize,
                 MPI_COMM_WORLD);

        MPI_Send(statsArray, /* message buffer */
                 5, /* one data item */
                 MPI_DOUBLE, /* data item is an integer */
                 0, /* destination process rank */
                 kStats, /* user chosen message tag */
                 MPI_COMM_WORLD); /* default communicator */
      }
#endif
    }
  }

  void Scope::error(std::stringstream& string)
  {
   LLOG(scopeLogger,error) << mulog::format::red << string.str();
  }

  void Scope::stats(std::stringstream& string)
  {
   LLOG(scopeLogger,verbose2) << string.str();
  }

  void Scope::msg(std::stringstream& string)
  {
   LLOG(scopeLogger,info) << string.str();
  }

  void Scope::verb(std::stringstream& string)
  {
   LLOG(scopeLogger,verbose) << string.str();
  }

  void Scope::error(const std::string& string)
  {
    LLOG(scopeLogger,error) <<  mulog::format::red <<  string;
  }

  void Scope::stats(const std::string& string)
  {
    LLOG(scopeLogger,verbose2) << string;
  }

  void Scope::msg(const std::string& string)
  {
    LLOG(scopeLogger,info) << string;
  }

  void Scope::verb(const std::string& string)
  {
    LLOG(scopeLogger,verbose) << string;
  }


  bool Scope::fexists(std::string & filename)
  {
    std::ifstream ifile(filename.c_str());
    return ifile;
  }
}


