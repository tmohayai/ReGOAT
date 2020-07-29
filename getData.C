/*************************************************************
 * authors: Tanaz Mohayai (mtanaz@fnal.gov) 
 * To run, open ROOT and enter the following commands:
 * root -b
 * root [0] .L getData.C
 * root [1] getData("inputfiles", "outfilename")
 * example of outfilename: signals.txt (see the process.sh bash script)
*************************************************************/

//some standard C++ includes
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <vector>

/*//some ROOT includes
#include "TInterpreter.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TStopwatch.h"
#include "TString.h"
*/
//"art" includes (canvas, and gallery)
#include "canvas/Utilities/InputTag.h"
#include "gallery/Event.h"
#include "gallery/ValidHandle.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindOne.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"

//global struct for hits
struct hit_t {
    Int_t runID;
    Int_t sID;
    Int_t eID;
    Int_t nCh;
    Int_t  nSamp;
    vector<Int_t> chID;
    vector<Float_t> E;
    vector<Float_t> noise;
    vector<Float_t> offset;
    vector<Int_t> maxPos;
    vector<Float_t> maxValue;
    vector<Int_t> minPos;
    vector<Float_t> minValue;
    vector<Int_t> tInit;
    vector<Int_t> tFin;
    vector<Float_t> E0;
    vector<Float_t> minSlope;
    vector<Int_t> tMinSlope;
    
    vector <UShort_t> t;
    vector<vector<Short_t>> signal;
    
    hit_t():runID(-1), sID(-1), eID(-1), nCh(-1), nSamp(-1) {};
    
    void reset() {
        runID=-1; sID=-1; eID=-1; nCh=-1; nSamp=-1;
        chID.clear(); E.clear(); noise.clear(); offset.clear() ;maxPos.clear(); maxValue.clear(); minPos.clear(); minValue.clear(); tInit.clear();tFin.clear(); E0.clear(); minSlope.clear(); tMinSlope.clear(); t.clear();
        for (size_t i = 0; i < signal.size();i++){signal[i].clear();}};
    
    void resize(const size_t n, const size_t nSig){chID.resize(n,0); E.resize(n,0); noise.resize(n,0); offset.resize(n,0); maxPos.resize(n,0); maxValue.resize(n,0); minPos.resize(n,0); minValue.resize(n,0); tInit.resize(n,0); tFin.resize(n,0); E0.resize(n,0); minSlope.resize(n,0); tMinSlope.resize(n,0); signal.resize(n); t.resize(n); for (size_t i = 0; i < n;i++){signal[i].resize(nSig,0);}};
    
    void resizeSignal(const unsigned int n) {signal.resize(n);};
    
    ~hit_t(){};
};

void createBranches(TTree* hitSumm, hit_t &hit, const std::vector <Int_t> &chIDs){
    hitSumm->Branch("runID", &(hit.runID), "runID/I");
    hitSumm->Branch("sID", &(hit.sID), "sID/I");
    hitSumm->Branch("eID", &(hit.eID), "eID/I");
    hitSumm->Branch("nCh",&(hit.nCh),"nCh/I");
    hitSumm->Branch("nSamp", &(hit.nSamp), "nSamp/I");
    hitSumm->Branch("E", &(hit.E[0]), "E[nCh]/F");
    hitSumm->Branch("chID", &(hit.chID[0]), "chID[nCh]/I");
    hitSumm->Branch("noise", &(hit.noise[0]), "noise[nCh]/F");
    hitSumm->Branch("offset", &(hit.offset[0]), "offset[nCh]/F");
    hitSumm->Branch("maxPos", &(hit.maxPos[0]), "maxPos[nCh]/I");
    hitSumm->Branch("maxValue", &(hit.maxValue[0]), "maxValue[nCh]/F");
    hitSumm->Branch("minPos", &(hit.minPos[0]), "minPos[nCh]/I");
    hitSumm->Branch("minValue", &(hit.minValue[0]), "minValue[nCh]/F");
    hitSumm->Branch("tInit", &(hit.tInit[0]), "tInit[nCh]/I");
    hitSumm->Branch("tFin", &(hit.tFin[0]), "tFin[nCh]/I");
    hitSumm->Branch("E0", &(hit.E0[0]), "E0[nCh]/F");
    hitSumm->Branch("minSlope", &(hit.minSlope[0]), "minSlope[nCh]/F");
    hitSumm->Branch("tMinSlope", &(hit.tMinSlope[0]), "tMinSlope[nCh]/I");
    
    int recordLength = 3072;
    stringstream tbranchName;
    tbranchName << "t["<<recordLength<<"]/s";
    hitSumm->Branch("t", &(hit.t[0]), tbranchName.str().c_str());
    for (int i=0; i<chIDs.size(); i++) {
        stringstream branchName, branchSpec;
        branchName <<"signal_" << chIDs[i];
        branchSpec << branchName.str().c_str() << "[" <<recordLength << "]/S";
        hitSumm->Branch(branchName.str().c_str(), &(hit.signal), branchSpec.str().c_str());
        //   hitSumm->Branch("signal", &(hit.signal), "signal[3,3072]/I");
    }
    
}

// Function to get root file directly *or* get root file names from a text file
#include <boost/algorithm/string/predicate.hpp>
std::vector<std::string> GetFileList(std::string input_file)
{
    std::vector<std::string> filenames;
    if(boost::algorithm::ends_with(input_file,".root")){
        filenames.emplace_back(input_file);
        std::cout << "Added file " << filenames.back() << " to be processed." << std::endl;
    }
    else{
        std::ifstream input(input_file);
        for(std::string line; std::getline(input,line);){
            filenames.emplace_back(line);
            std::cout << "Added file " << filenames.back() << " to be processed." << std::endl;
        }
    }
    return filenames;
}

void refreshTreeAddresses(TTree *hitSumm, hit_t &hit, const std::vector <Int_t> &chIDs)
{
    hitSumm->SetBranchAddress("runID",  &(hit.runID));
    hitSumm->SetBranchAddress("sID",  &(hit.sID));
    hitSumm->SetBranchAddress("eID", &(hit.eID));
    hitSumm->SetBranchAddress("nCh",&(hit.nCh));
    hitSumm->SetBranchAddress("nSamp", &(hit.nSamp));
    hitSumm->SetBranchAddress("chID", &(hit.chID[0]));
    hitSumm->SetBranchAddress("E", &(hit.E[0]));
    hitSumm->SetBranchAddress("noise", &(hit.noise[0]));
    hitSumm->SetBranchAddress("offset", &(hit.offset[0]));
    hitSumm->SetBranchAddress("maxPos",&(hit.maxPos[0]));
    hitSumm->SetBranchAddress("maxValue",&(hit.maxValue[0]));
    hitSumm->SetBranchAddress("minPos",&(hit.minPos[0]));
    hitSumm->SetBranchAddress("minValue",&(hit.minValue[0]));
    hitSumm->SetBranchAddress("tInit",&(hit.tInit[0]));
    hitSumm->SetBranchAddress("tFin",&(hit.tFin[0]));
    hitSumm->SetBranchAddress("E0",&(hit.E0[0]));
    hitSumm->SetBranchAddress("minSlope",&(hit.minSlope[0]));
    
    hitSumm->SetBranchAddress("t",&(hit.t[0]));
    for (int i=0; i<hit.signal.size(); i++) {
        stringstream branchName, branchSpec;
        branchName <<"signal_" << chIDs[i];
        hitSumm->SetBranchAddress(branchName.str().c_str(),&(hit.signal[i][0]));
    }
}



void ProcessEvent(gallery::Event *ev, TFile *outfile, std::string filenamebase, const std::vector <Int_t> &chIDs,TTree *hitSumm, TTree *config, hit_t &hit)
{
	std::cout << "Processing "
	<< "Run " << ev->eventAuxiliary().run() << ", "
	<< "Event " << ev->eventAuxiliary().event() << std::endl;
	art::InputTag const fDigitModuleLabel("daq");

	auto const& digitVecHandle = ev->getValidHandle<std::vector<raw::RawDigit>>(fDigitModuleLabel);

 	const int nchannels = (int)digitVecHandle->size();
    	std::cout << "Size of digitVecHandle is " << nchannels << std::endl;
	unsigned int dataSize = 0; // size of raw data vectors
    	raw::ChannelID_t channel = raw::InvalidChannelID; // channel number
    	unsigned int bin(0); // time bin loop variable
	const int recordLength = 3072;
    	std::vector<short> rawadc(recordLength);
 	TH1F *waveforms[nchannels];

	for (int i_ch = 0; i_ch < nchannels; i_ch++) {
        
        auto const& digitVec = (*digitVecHandle)[i_ch];
        
        channel = digitVec.Channel();
        dataSize = digitVec.Samples();
        rawadc.resize(dataSize);
        TString histname;
    	TString histtitle;
        waveforms[i_ch] = new TH1F(histname.Data(),histtitle.Data(),dataSize+1,-0.5,dataSize+0.5);
        
        raw::Uncompress(digitVec.ADCs(), rawadc, digitVec.Compression());
        for (bin = 0; bin < dataSize; bin++){
            	std::cout << "rawadc[" << bin << "] : " << std::endl;
		waveforms[i_ch]->SetBinContent(bin+1,rawadc[bin]);
        }
        
	}

	int evtnum = ev->eventAuxiliary().event();
    	TString dirname = Form("evt_%d/",evtnum);
	
	Double_t evEnergy = -1;
 	Double_t noiseStd = -1;
    	Double_t sigOffset = -1;

    	int bin_ch;
    	double canvmax = 1.3;
    	double canvmin = -1.3;

	hit.reset();
    	hit.resize(chIDs.size(),dataSize);
    	hit.eID = evtnum;
    	hit.nCh = chIDs.size();
    	hit.nSamp = dataSize;
    	//for (UShort_t it = 0; it<dataSize; it++) {
        //	hit.t[it] = it;
	//}


    for (size_t i_ch = 0; i_ch<chIDs.size();i_ch++) {


        hit.chID[i_ch] = chIDs[i_ch];
        canvmax = 2;
        canvmin = -2;

 	Int_t nBins = waveforms[chIDs[i_ch]]->GetNbinsX();
        Int_t finBinNoise = nBins*2/10;
        Int_t initBinNoise = 1;
        sigOffset = waveforms[chIDs[i_ch]]->Integral(initBinNoise,finBinNoise)/(finBinNoise - initBinNoise + 1);
        evEnergy = fabs(waveforms[chIDs[i_ch]]->Integral() - nBins*sigOffset);
        hit.E[i_ch] = evEnergy;;
        hit.offset[i_ch] = sigOffset;

	for (Int_t i = initBinNoise; i < finBinNoise + 1; ++i) {
            noiseStd += pow(waveforms[chIDs[i_ch]]->GetBinContent(i) - sigOffset,2);
        }
	noiseStd = sqrt(noiseStd/double(finBinNoise-initBinNoise+1));
        hit.noise[i_ch] = noiseStd;

	Int_t maxPos = waveforms[chIDs[i_ch]]->GetMaximumBin();
        Double_t maxValue = waveforms[chIDs[i_ch]]->GetBinContent(maxPos);
        hit.maxPos[i_ch] = maxPos;
        hit.maxValue[i_ch] = maxValue;
        Int_t minPos = waveforms[chIDs[i_ch]]->GetMinimumBin();
        Double_t minValue = waveforms[chIDs[i_ch]]->GetBinContent(minPos);
        hit.minPos[i_ch] = minPos;
        hit.minValue[i_ch] = minValue;
	Int_t tInit = waveforms[chIDs[i_ch]]->FindFirstBinAbove(hit.offset[i_ch]-3*hit.noise[i_ch]);
        Int_t tFin = waveforms[chIDs[i_ch]]->FindLastBinAbove(hit.offset[i_ch]-3*hit.noise[i_ch]);

	Bool_t flagAboveThres = false;
        for (Int_t i = 1; i < nBins+1; ++i) {
            if (waveforms[chIDs[i_ch]]->GetBinContent(i)<hit.offset[i_ch]-5*hit.noise[i_ch]){
                if (!flagAboveThres){
                    tInit = i;
                    flagAboveThres = true;
                }
                tFin = i;
            }
        }
	Double_t E0 = fabs(waveforms[chIDs[i_ch]]->Integral(tInit,tFin,"width") - (tFin-tInit+1)*sigOffset);
        hit.tInit[i_ch] = tInit;
        hit.tFin[i_ch] = tFin;
        hit.E0[i_ch] = E0;
	Double_t slope = 0;
        Double_t minSlope = 0;
        Int_t minSlopePos = 0;
        for (Int_t i = 1; i < hit.nSamp-1; ++i) {
            slope = waveforms[chIDs[i_ch]]->GetBinContent(i+1) - waveforms[chIDs[i_ch]]->GetBinContent(i);
            if (slope < minSlope) {
                minSlope = slope;
                minSlopePos = i;
            }
        }
	hit.minSlope[i_ch] = minSlope;
        hit.tMinSlope[i_ch] = minSlopePos;
	for (int bin = 0; bin < hit.nSamp; bin++){
            hit.signal[i_ch][bin] = (Short_t) waveforms[chIDs[i_ch]]->GetBinContent(bin+1);
        }

    }

    refreshTreeAddresses(hitSumm,hit, chIDs);
    hitSumm->Fill();

    for (int i_ch = 0; i_ch < nchannels; i_ch++) {
        delete waveforms[i_ch];
    }



}
void getData(std::string input_files, std::string outfilename="data.root")
{
    gStyle->SetOptStat(0);
    
    // Make output file
    TFile *outfile = new TFile(outfilename.c_str(),"RECREATE");
    std::string filenamebase = outfilename.substr(0, outfilename.find("."));
    outfile->cd();
    TTree *hitSumm = new TTree("hitSumm","hitSumm");
    TTree *config = new TTree("config","config");
    
    //Channel of interest
    std::vector <Int_t> chIDs;
    
    chIDs.push_back(33);
    
    //create event struct and tree branches
    hit_t hit;
    size_t nTemp = 3072;
    hit.resize(chIDs.size(),nTemp);
    createBranches(hitSumm,hit, chIDs);
    
    
    // Format our files list
    std::vector<std::string> filenames = GetFileList(input_files);
    
    // Make a gallery::Event
    gallery::Event ev(filenames);
    
    
    // Now loop through all events and process
    while (!ev.atEnd()){
        ProcessEvent(&ev, outfile, filenamebase, chIDs,hitSumm,config, hit);
        ev.next();
    }
    
    outfile->cd();
    hitSumm->SetBranchStatus("*",1);
    hitSumm->Write("",TObject::kOverwrite);
    outfile->Close();
    
    //delete outfile;
    outfile->Delete();
    cout << "DONE!" << endl;
}

