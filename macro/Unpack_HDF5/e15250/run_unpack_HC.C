#define cRED "\033[1;31m"
#define cYELLOW "\033[1;33m"
#define cNORMAL "\033[0m"
#define cGREEN "\033[1;32m"

struct auxchannel
{
  std::string name;
  uint8_t cobo;
  uint8_t asad;
  uint8_t aget;
  uint8_t channel;
};

void run_unpack_HC(std::string dataFile = "/mnt/analysis/e15250_attpc/h5/run_0250.h5",
                   TString parameterFile = "ATTPC.e15250_sim.par", TString mappath = "")
{

  // -----   Timer   --------------------------------------------------------
 TStopwatch timer;
 timer.Start();
 // ------------------------------------------------------------------------

  gSystem->Load("libXMLParser.so");
  // -----------------------------------------------------------------
  // Set file names
  TString scriptfile = "Lookup20150611.xml";
  TString dir = getenv("VMCWORKDIR");
  TString scriptdir = dir + "/scripts/"+ scriptfile;
  TString dataDir = dir + "/macro/data/";
  TString geomDir = dir + "/geometry/";
  gSystem -> Setenv("GEOMPATH", geomDir.Data());

  //TString inputFile   = dataDir + name + ".digi.root";
  //TString outputFile  = dataDir + "output.root";
  TString outputFile = "output.root";
  //TString mcParFile   = dataDir + name + ".params.root";
  TString loggerFile  = dataDir + "ATTPCLog.log";
  TString digiParFile = dir + "/parameters/" + parameterFile;
  TString geoManFile  = dir + "/geometry/ATTPC_v1.1.root";

  TString inimap   = mappath + "inhib.txt";
  TString lowgmap  = mappath + "lowgain.txt";
  TString xtalkmap = mappath + "beampads_e15503b.txt";

  // -----------------------------------------------------------------
  // Logger
  /*FairLogger *fLogger = FairLogger::GetLogger();
  fLogger -> SetLogFileName(loggerFile);
  fLogger -> SetLogToScreen(kTRUE);
  fLogger -> SetLogToFile(kTRUE);
  fLogger -> SetLogVerbosityLevel("LOW");*/

  FairRunAna* run = new FairRunAna();
  run -> SetOutputFile(outputFile);
  //run -> SetGeomFile("../geometry/ATTPC_Proto_v1.0.root");
  run -> SetGeomFile(geoManFile);

  FairRuntimeDb* rtdb = run->GetRuntimeDb();
  FairParAsciiFileIo* parIo1 = new FairParAsciiFileIo();
  parIo1 -> open(digiParFile.Data(), "in");
  //FairParRootFileIo* parIo2 = new FairParRootFileIo();
  //parIo2 -> open("param.dummy_proto.root");
 // rtdb -> setFirstInput(parIo2);
  rtdb -> setSecondInput(parIo1);

  //Auxiliary channels
  //Hash table: cobo, asad, aget, channel
  std::vector<auxchannel> aux_channels;

  auxchannel ch_1{"ion_chamber",4,1,0,65};
  aux_channels.push_back(ch_1);
  auxchannel ch_2{"ion_chamber_downscale",4,1,0,66};
  aux_channels.push_back(ch_2);
  auxchannel ch_3{"unknown_1",4,1,0,67};
  aux_channels.push_back(ch_3);
  auxchannel ch_4{"unknown_2",4,1,0,61};
  aux_channels.push_back(ch_4);
  auxchannel ch_5{"unknown_3",4,1,0,64};
  aux_channels.push_back(ch_5);
  auxchannel ch_6{"unknown_4",4,1,0,59};
  aux_channels.push_back(ch_6);

  AtHDFParserTask *HDFParserTask = new AtHDFParserTask();
  HDFParserTask->SetPersistence(kFALSE);
  HDFParserTask->SetAtTPCMap(scriptdir.Data());
  HDFParserTask->SetFileName(dataFile);
  HDFParserTask->SetBaseLineSubtraction(kTRUE);

  for (auto iaux : aux_channels) {
     auto hash = HDFParserTask->CalculateHash(iaux.cobo, iaux.asad, iaux.aget, iaux.channel);
     auto isaux = HDFParserTask->SetAuxChannel(hash, iaux.name);
  }

  AtPSASimple2 *psa = new AtPSASimple2();
  // psa -> SetPeakFinder(); //NB: Use either peak finder of maximum finder but not both at the same time
  // psa -> SetBaseCorrection(kFALSE);
  // psa -> SetTimeCorrection(kFALSE);

  // AtPSAFilter *psa = new AtPSAFilter();

  AtPSAtask *psaTask = new AtPSAtask(psa);
  psaTask->SetPersistence(kTRUE);
  psa->SetThreshold(30);
  psa->SetMaxFinder();
  // psa->SetMeanK(4);
  // psa->SetStddevMulThresh(0.1);

  // AtPRATask *praTask = new AtPRATask();
  // praTask -> SetPersistence(kTRUE);

  run -> AddTask(HDFParserTask);
  run -> AddTask(psaTask);
  //run -> AddTask(praTask);

  run -> Init();

  auto numEvents = HDFParserTask->GetNumEvents() / 2;

  run->Run(0, numEvents);
  // run->Run(0, 1000);
  // run -> RunOnTBData();

  std::cout << std::endl << std::endl;
  std::cout << "Macro finished succesfully."  << std::endl << std::endl;
  std::cout << "- Output file : " << outputFile << std::endl << std::endl;
  // -----   Finish   -------------------------------------------------------
  timer.Stop();
  Double_t rtime = timer.RealTime();
  Double_t ctime = timer.CpuTime();
  cout << endl << endl;
  cout << "Real time " << rtime << " s, CPU time " << ctime << " s" << endl;
  cout << endl;
  // ------------------------------------------------------------------------

  gApplication->Terminate();

}

