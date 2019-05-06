#include <iostream>
#include "TROOT.h"
#include "TH1F.h"
#include "TF1.h"
#include "TFile.h"
#include "TMath.h"
#include "TProfile2D.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TChain.h"
#include "TTree.h"
#include "TGraph.h"
#include "TBox.h"

// ==> m.galli added 12jul17, 
#include "TFormula.h"


using namespace std;

void plotEvent(int i, double sigma, TTree *tbar);

void tgf_dead_time_2_1(int   nrdm        = 2, 
                       char *outfilename = "./sim.root", 
                       int   idoffset    = 0, 
                       float minFluence  = 0.05, 
                       float maxFluence  = 0.07, 
                       float minT50      = 0.02, 
                       float maxT50      = 0.05, 
                       int   verbose     = 1, 
                       int   sourcetype  = 3, 
                       float rate        = 20000.0,
                       char *timeprofile = "./uniform_split10.root",
                       float aeff=150.)
{
    //gROOT->Reset();
    //ROOT::v5::TFormula::SetMaxima(10000); // m.galli changed 29nov17,
    TFormula::SetMaxima(10000);             // m.galli changed 29nov17,
    TRandom3 *rnd = new TRandom3();
        
    rnd->SetSeed(); // m.galli for tests 29nov17
    //rnd->SetSeed(4357);   // m.galli for tests 29nov17
    
//    int ncal = 157; // number of expected events in MCAL above threshold for theta 60, phi 0. NB: must use correct value by MG
//    int nevt = 40924; // number of simulated photons
//    double area = 58105; // on this area (cm^2)
//    double simfluence = nevt / area; // simulated fluence
    double dtbar = 0.020; // MCAL bar blanking time window (ms) NB: energy-dependent dead time to be implemented
    double dtpileup = 0.0045; // MCAL pile-up time window provided by PDHU (ms)
    
    // MCAL bars parameters
    int nbars = 26; // active MCAL bars
    float ly = 21000.; // average light yield (e-/MeV)
    float gain = 5.e-7; // amplifier gain in V/e- after PB spice simulations
    float barlen  =   37.5; // bar length
    float hb = 0.5*barlen; // half bar length
    float alpha = -0.028; // average attenuation (cm^-1)
    float auxP = -1./(2.*alpha); // coefficient for position (cm) calculation
    float auxE = exp(-alpha*hb)/(gain*ly); // coefficient for energy (MeV) calculation
    
    // FEE electronics parameters
    float dt = 0.0001; // time sample (ms) for total signal sampling
    float dv = 0.01; // signal drop (V) needed to issue a trigger on peak detected
    float thr = 0.00224; // arming discriminator threshold (for opening and closing a trigger)
    
    // DFE electronics (FPGA) parameters (see email by PB on 31/07/2017)
    float tjitter = 0.002; // time window (default 2 us) to join bar signals in a single count
    float treject = 0.001; // time window (default 1 us) where triggers are not accepted
    float tblank  = 0.020; // time window (default 20 us) where triggers on the same bar are not accepted
    float tckADC  = 0.00114; // BURST chain ADC clock (us)
    float tlatency = 2.5; // BURST chain ADC latency time (clock periods)
    
    struct MCALevt {
        Double_t	time;
        Int_t      	mult;
        Short_t	bar[30];
        Short_t 	pdA[30];
        Short_t 	pdB[30];
        Float_t	pos[30];
        Float_t    	E[30];
        Float_t	Etot;
        Float_t	E0;  // simulated energy (if applicable)
   };
    MCALevt 	inevt;
    int 		nevents;
    
    struct MCALbar {
        Double_t	time;
        Int_t       ievt;
        Int_t	bar;
        Float_t 	pdA;
        Float_t 	pdB;
        Float_t	pos;
        Float_t    	E;
    };
    MCALbar 	inbar, outbar;
    char tcut[100];
    
    // create the chain of input trees
    TChain data("tdata");
    int nfile = 0;
    
    switch (sourcetype) {
        case 2:
            nfile += data.Add("/Users/dub002/dev/MCAL-front-end-simulator/fake.root");
            break;
 
        case 3:
            nfile += data.Add(timeprofile);            
            break;
           
        default:
            nfile += data.Add("/Users/dub002/data/TGF/triggers_with_WWLLN/MC_simulations/th047_ph290/TGF41313-25_cutpl02-55.root");
            nfile += data.Add("/Users/dub002/data/TGF/triggers_with_WWLLN/MC_simulations/th047_ph290/TGF41313-25_cutpl02-55bis.root");
            break;
    }
    
    // connect the data input tree
    data.SetBranchAddress("time", &inevt.time);
    data.SetBranchAddress("mult", &inevt.mult);
    data.SetBranchAddress("bar", inevt.bar);
    data.SetBranchAddress("pdA", inevt.pdA);
    data.SetBranchAddress("pdB", inevt.pdB);
    data.SetBranchAddress("pos", inevt.pos);
    data.SetBranchAddress("E", inevt.E);
    data.SetBranchAddress("Etot", &inevt.Etot);
    data.SetBranchAddress("E0", &inevt.E0); // simulated energy (if applicable)
    nevents = data.GetEntries();
    if (verbose) cout << nevents << " events to be processed in " << nfile << " files" << endl;
    
    // open file for output and build output trees
    TFile fout(outfilename, "RECREATE"); // "UPDATE"
    int id; // simulation id
    int nc; // number of counts
    int ncout; // detected number of counts
    
    TTree *tin = new TTree("tin", "input events");
    tin->Branch("id", &id, "id/I");
    tin->Branch("time", &inevt.time, "time/D");
    tin->Branch("mult", &inevt.mult, "mult/I");
    tin->Branch("bar", inevt.bar, "bar[30]/S");
    tin->Branch("pdA", inevt.pdA, "pdA[30]/S");
    tin->Branch("pdB", inevt.pdB, "pdB[30]/S");
    tin->Branch("pos", inevt.pos, "pos[30]/F");
    tin->Branch("E", inevt.E, "E[30]/F");
    tin->Branch("Etot", &inevt.Etot, "Etot/F");
    tin->Branch("E0", &inevt.E0, "E0/F"); // simulated energy (if applicable)
   
    TTree *tdet = new TTree("tdet", "output measured events");
    tdet->Branch("id", &id, "id/I");
    tdet->Branch("time", &inevt.time, "time/D");
    tdet->Branch("mult", &inevt.mult, "mult/I");
    tdet->Branch("bar", inevt.bar, "bar[30]/S");
    tdet->Branch("pdA", inevt.pdA, "pdA[30]/S");
    tdet->Branch("pdB", inevt.pdB, "pdB[30]/S");
    tdet->Branch("pos", inevt.pos, "pos[30]/F");
    tdet->Branch("E", inevt.E, "E[30]/F");
    tdet->Branch("Etot", &inevt.Etot, "Etot/F");
    
    float inMedianE, inMaxE, outMedianE, outMaxE;
    float inRMS = 0.;
    float outRMS = 0.;
    double fluence =0.;
    double t50=0.;
    double sigma=0.;
    TTree *tsum = new TTree("tsum", "summary");
    tsum->Branch("id", &id, "id/I");
    tsum->Branch("fluence", &fluence, "fluence/D");
    tsum->Branch("t50", &t50, "t50/D");
    tsum->Branch("nc", &nc, "nc/I");
    tsum->Branch("ncout", &ncout, "ncout/I");
    tsum->Branch("inMedianE", &inMedianE, "inMedianE/F");
    tsum->Branch("outMedianE", &outMedianE, "outMedianE/F");
    tsum->Branch("inMaxE", &inMaxE, "inMaxE/F");
    tsum->Branch("outMaxE", &outMaxE, "outMaxE/F");
    tsum->Branch("inRMS", &inRMS, "inRMS/F");
    tsum->Branch("outRMS", &outRMS, "outRMS/F");
    
    TObjArray *fcoll = new TObjArray(1000); // collection to store functions, to ease memory handling
    
    // main loop on random extractions
    
    for (int i=0; i<nrdm; i++) {
        //if (i%5000 == 0) cout << "Processed " << i << " events" << endl;
        id = idoffset + i;
        double *atcal = NULL; // store events time
        int *atsim = NULL; // store simulated event entry number
        
        switch (sourcetype) {
                
            case 1: // constant rate, uniform distribution in time, for validation of the DFE model by PB
            {
                float twindow = 10.; // time interval per test (ms)
                nc = (int) rnd->PoissonD(0.001*rate*twindow);
                atcal = new double[nc]; // MCAL time array, prior to any dead time correction
                for (int h=0; h<nc; h++) atcal[h] = rnd->Uniform(0., twindow); // change here the TGF time profile as liked
                printf("\n ---- Event: %5d  Id = %5d Uniform distribution with rate = %.3f Hz, MCAL counts in %f ms interval = %d ----- \n", i,id, rate, twindow, nc);
                break;
            }
            
            case 2: // fake events on single bar, for validation, for validation of the FEE model by PB
            {
                nc = 2;
                atcal = new double[nc]; // MCAL time array, prior to any dead time correction
                float tsep = 0.002;
                atcal[0] = 0.;
                atcal[1] = (i+1)*tsep;
                printf("\n ---- Event: %5d  Id = %5d Fake test for validation: n. of counts %d - t separation = %f ms ----- \n", i,id, nc, tsep);
                break;
            }
                
            case 3: // TGF time profile as simulated by D. Sarria and NØ, extract uniform random fluence, extract uniform entries for (energy,time) correlated values
            {
                fluence = rnd->Uniform(minFluence, maxFluence); //(0.01, 0.1)
                t50 = rnd->Uniform(minT50, maxT50); // (0.01, 2.)
                sigma = t50*0.74;
                nc = (int) (fluence*aeff); // expected number of MCAL counts
                atcal = new double[nc]; // MCAL time array, prior to any dead time correction
                atsim = new int[nc]; // simulated event entry number
                for (int h=0; h<nc; h++) {
                    atsim[h] = rnd->Integer(nevents);
                    data.GetEntry(atsim[h]); // extract random event from simulations
                    atcal[h] = 0.001*inevt.time; // this is because in root file is in microsec, while in this program I use ms
               
                    // introduce a gaussian broadening for times
                    if (maxT50>0.) {
                        double tbroad=rnd->Gaus(0.,sigma);
                        atcal[h] += tbroad ;
                    }
                }
                printf("\n ---- Event: %5d  Id = %5d TGF Fluence = %.3f, MCAL counts = %d ----- \n", i,id, fluence, nc); //if (verbose)
                break;
            }
           
            default: // Gaussian shaped TGF, extract uniform random fluence and duration
            {
                fluence = rnd->Uniform(minFluence, maxFluence); //(0.01, 0.1)
                t50 = rnd->Uniform(minT50, maxT50); // (0.01, 2.)
                sigma = t50*0.74;
                nc = (int) (fluence*aeff); // expected number of MCAL counts
                atcal = new double[nc]; // MCAL time array, prior to any dead time correction
                for (int h=0; h<nc; h++) atcal[h] = rnd->Gaus(0., sigma);
                printf("\n ---- Event: %5d  Id = %5d TGF Fluence = %.3f, T50 = %.3f , MCAL counts = %d ----- \n", i,id, fluence, t50, nc); //if (verbose)
                break;
            }
        }
        if (nc>1) inRMS = TMath::RMS(nc, atcal);
        
            
        // create the simulated TGF as seen by MCAL using simulations
        
        //TTree *tgf  = data.CloneTree(0); // input tgf
        TTree *tgf = new TTree("tgf", "input events");
        tgf->Branch("id", &id, "id/I");
        tgf->Branch("time", &inevt.time, "time/D");
        tgf->Branch("mult", &inevt.mult, "mult/I");
        tgf->Branch("bar", inevt.bar, "bar[30]/S");
        tgf->Branch("pdA", inevt.pdA, "pdA[30]/S");
        tgf->Branch("pdB", inevt.pdB, "pdB[30]/S");
        tgf->Branch("pos", inevt.pos, "pos[30]/F");
        tgf->Branch("E", inevt.E, "E[30]/F");
        tgf->Branch("Etot", &inevt.Etot, "Etot/F");
        tgf->Branch("E0", &inevt.E0, "E0/F");
      
        // sort MCAL time array
        int *idc = new int[nc];
        TMath::Sort(nc, atcal, idc, 0);
        if (verbose) for (int h=0; h<nc; h++) printf("%10.6f\n", atcal[idc[h]]);
        float tmin = atcal[idc[0]];  // time of first count
        float tmax = atcal[idc[nc-1]]; // time of last count
        
        TTree *tbar = new TTree("tbar","");
        tbar->Branch("time", &inbar.time, "time/D");
        tbar->Branch("ievt", &inbar.ievt, "ievt/I");
        tbar->Branch("bar", &inbar.bar, "bar/I");
        tbar->Branch("pdA", &inbar.pdA, "pdA/F");
        tbar->Branch("pdB", &inbar.pdB, "pdB/F");
        tbar->Branch("pos", &inbar.pos, "pos/F");
        tbar->Branch("E", &inbar.E, "E/F");
        
        for (int h=0; h<nc; h++) {
            switch (sourcetype) {
                case 2:
                    data.GetEntry(150);
                    break;
                    
                case 3: 
                {
                    int ievent = atsim[idc[h]];
                    data.GetEntry(ievent); // extract random event from simulations
                    break;
                }  
                default: 
                {
                    int ievent = rnd->Integer(nevents);
                    data.GetEntry(ievent); // extract random event from simulations
                    break;
                }
            }
            inevt.time = atcal[idc[h]]; // update event time
            tgf->Fill(); // save to new tree
            tin->Fill(); // save to input tree
            
            // split the event for each fired bar, fill the bar tree
            
            int ibar;
            inbar.time = inevt.time;
            inbar.ievt = h;
            for (int j=0; j<30; j++) {
                ibar = inevt.bar[j];
                if (ibar>=0) {
                    inbar.bar = j;
                    inbar.pos = inevt.pos[ibar];  // check reference frame with MG
                    inbar.E = inevt.E[ibar];
                    inbar.pdA = gain*ly*inbar.E * exp(alpha*(inbar.pos + hb)); //  output in V, includes ly and e/V conversion
                    inbar.pdB = gain*ly*inbar.E * exp(alpha*(hb - inbar.pos));
                    tbar->Fill();
                }
            }
        }
        tbar->SetScanField(0);
        
        if (verbose) tbar->Scan("time:bar:E");
        if (verbose) {
            for (int h=0; h<nc; h++) {
                tgf->Show(h, 30);
                tgf->GetEntry(h);
                sprintf(tcut, "ievt == %d", h);
                tbar->Scan("time:bar:pdA:pdB:pos:E", tcut);
            }
        }
        //plotEvent(i, sigma, tbar);
        
        // build output tree for bar entries after trigger generation
        TTree *toutbar = new TTree("toutbar","");
        toutbar->Branch("time", &outbar.time, "time/D");
        toutbar->Branch("ievt", &outbar.ievt, "ievt/I");
        toutbar->Branch("bar", &outbar.bar, "bar/I");
        toutbar->Branch("pdA", &outbar.pdA, "pdA/F");
        toutbar->Branch("pdB", &outbar.pdB, "pdB/F");
        toutbar->Branch("pos", &outbar.pos, "pos/F");
        toutbar->Branch("E", &outbar.E, "E/F");
        toutbar->SetScanField(0);
        
        char cname[100], cut[100], hName[100], fName[100];
        for (int j=0; j<30; j++) {
            
            sprintf(cname, "c%06d_%02d.png", i, j);
            TCanvas *c1 = new TCanvas(cname, cname, 1200, 800);
            c1->Divide(1,2);
            
            TF1 *fAtot = NULL;
            TF1 *fBtot = NULL;
            TF1 *fABtot = NULL;
            int n;
            double *t = NULL;

            
            // work on pdA
            c1->cd(1);
            sprintf(hName, "h2_%06d_bar%02d_A", i, j);
            TH2F *h2A = new TH2F(hName, hName, 1000, tmin -0.1, tmax +0.1, 1000, 0., 1.);
            h2A->SetStats(0);
            h2A->Draw();
            sprintf(cut, "bar == %d", j);
            tbar->Draw("pdA:time", cut, "psame");
            if (gPad->GetPrimitive("Graph") != NULL) {
                TGraph *gA = (TGraph *) gPad->GetPrimitive("Graph");
                gA->SetMarkerStyle(5);
                
                n = gA->GetN();
                t = gA->GetX();
                double *v = gA->GetY();
                TString fAtotDef("");
                
                // build the signal functions
                for (int k=0; k<n; k++) {
                    sprintf(fName, "fA%02d%04d", j, k);
                    TF1 *f = new TF1(fName, "x< [2] ? [0] : [0] + [1]*((x-[2])/[3])^([4])*exp(-(x-[2])/[3])", t[0], t[n-1]+0.1);
                    f->SetParameters(0., v[k]/2.75, t[k], 3.32432e-03, 3.61901); // no offset
                    f->SetLineColor(2);
                    f->Draw("lsame");
                    fcoll->Add(f);
                    
                    // add function name to total signal definition
                    if (k==0) fAtotDef += fName;
                    else {
                        fAtotDef += " + ";
                        fAtotDef += fName;
                    }
                }
                
                // build the sum signal function for side A
                fAtot = new TF1("fAtot",fAtotDef.Data(), t[0], t[n-1]+0.1);
                fAtot->SetLineColor(kBlue);
                //fAtot->SetLineStyle(2);
                fAtot->SetNpx(2000);
                fAtot->Draw("lsame");
                fcoll->Add(fAtot);
            }
            
            // work on pdB
            c1->cd(2);
            sprintf(hName, "h2_%06d_bar%02d_B", i, j);
            TH2F *h2B = new TH2F(hName, hName, 1000,  tmin -0.1, tmax +0.1, 1000, 0., 1.);
            h2B->SetStats(0);
            h2B->Draw();
            tbar->Draw("pdB:time", cut, "psame");
            if (gPad->GetPrimitive("Graph") != NULL) {
                TGraph *gB = (TGraph *) gPad->GetPrimitive("Graph");
                gB->SetMarkerStyle(5);
                
                n = gB->GetN();
                t = gB->GetX();
                double *v = gB->GetY();
                TString fBtotDef("");
                
                // build the signal functions
                for (int k=0; k<n; k++) {
                    sprintf(fName, "fB%02d%04d", j, k);
                    TF1 *f = new TF1(fName, "x< [2] ? [0] : [0] + [1]*((x-[2])/[3])^([4])*exp(-(x-[2])/[3])", t[0], t[n-1]+0.1);
                    f->SetParameters(0., v[k]/2.75, t[k], 3.32432e-03, 3.61901); // no offset
                    f->SetLineColor(2);
                    f->Draw("lsame");
                    fcoll->Add(f);
                    
                    // add function name to total signal definition
                    if (k==0) fBtotDef += fName;
                    else {
                        fBtotDef += " + ";
                        fBtotDef += fName;
                    }
                }
                
                // build the sum signal function for side B
                fBtot = new TF1("fBtot",fBtotDef.Data(), t[0], t[n-1]+0.1);
                fBtot->SetLineColor(kBlue);
                //fBtot->SetLineStyle(2);
                fBtot->SetNpx(2000);
                fBtot->Draw("lsame");
                fcoll->Add(fBtot);
            }
            
            // build the sum signal function for side A+B
            if(fAtot && fBtot) {
                fABtot = new TF1("fABtot", "fAtot + fBtot", t[0], t[n-1]+0.1);
                fABtot->SetLineColor(kGreen);
                fABtot->SetNpx(2000);
                fcoll->Add(fABtot);
                
                // scan the sum signal and get peaks and trigger start and stop
                int arm=0;
                int trg=0;
                int ntrg = 0;
                float *ttrg = new float[n]; // trigger start time
                float *tstop = new float[n]; // trigger start time
                float *vAB = new float[n]; // A+B sum signal amplitude, discrimination is done against this signal
                float *vA = new float[n]; // A signal amplitude
                float *vB = new float[n]; // B signal amplitude
                for (int h=0; h<n; h++) vAB[h] = -1.; // init amplitude array
                
                for (float tt = t[0]; tt<t[n-1]+0.1; tt+=dt) {
                    //cout << tt << "  " << fABtot->Eval(tt)  << "  " << ntrg << "  " << arm  << "  " << trg << endl;
                    if (!arm && fABtot->Eval(tt)>=thr) { // discriminator is not armed, signal rise above threshold
                        arm = 1; // arm the discriminator
                    }
                    if (arm && !trg) { // if armed but trigger is not issued
                        if (fABtot->Eval(tt) >= vAB[ntrg]) { // follow the rising signal
                            vAB[ntrg] = fABtot->Eval(tt);
                            vA[ntrg]  = fAtot->Eval(tt);
                            vB[ntrg]  = fBtot->Eval(tt);
                            ttrg[ntrg] = tt;
                        } else {
                            if (vAB[ntrg] - fABtot->Eval(tt) > dv || fABtot->Eval(tt)<thr) { // valley found: issue a trigger on the peak
                                trg = 1;
                            }
                        }
                    }	else { // if trigger is already issued
                        if (arm && fABtot->Eval(tt) < thr) { // reset trigger and arming discriminator when signal goes below threshold
                            tstop[ntrg] = tt;
                            ntrg++;
                            arm=0;
                            trg=0;
                        }
                    }
                }
                
                if (verbose) {
                    printf("\nSim %d - bar %d:\nNumber of hits %d - number of triggers %d\n\n  n.    tstart    tstop      vA       vB\n", i, j, n, ntrg);
                    for (int h=0; h<ntrg; h++) printf("%3d  %10.3f  %10.3f  %10.3f  %10.3f\n", h, ttrg[h], tstop[h], vA[h], vB[h]);
                }
                
                // draw triggers, format plot and axis
                c1->cd(1);
                for (int h=0; h<ntrg; h++) {
                    TBox *tbA = new TBox(ttrg[h], 0., tstop[h], vA[h]);
                    tbA->SetFillStyle(3005);
                    tbA->SetFillColor(kBlue);
                    tbA->Draw("same");
                }
                h2A->SetAxisRange(0., 1.1*fAtot->GetMaximum(t[0], t[n-1]+0.1),"Y");
                h2A->GetXaxis()->SetLabelSize(0.05);
                h2A->GetXaxis()->SetTitle("Time (ms)");
                h2A->GetXaxis()->SetTitleSize(0.05);
                h2A->GetXaxis()->CenterTitle();
                h2A->GetYaxis()->SetLabelSize(0.05);
                h2A->GetYaxis()->SetTitle("Signal (V)");
                h2A->GetYaxis()->SetTitleSize(0.05);
                h2A->GetYaxis()->CenterTitle();
                gPad->SetTicks(1,1);
                gPad->Update();
                
                c1->cd(2);
                for (int h=0; h<ntrg; h++) {
                    TBox *tbB = new TBox(ttrg[h], 0., tstop[h], vB[h]);
                    tbB->SetFillStyle(3005);
                    tbB->SetFillColor(kBlue);
                    tbB->Draw("same");
                }
                h2B->SetAxisRange(0., 1.1*fBtot->GetMaximum(t[0], t[n-1]+0.1),"Y");
                h2B->GetXaxis()->SetLabelSize(0.05);
                h2B->GetXaxis()->SetTitle("Time (ms)");
                h2B->GetXaxis()->SetTitleSize(0.05);
                h2B->GetXaxis()->CenterTitle();
                h2B->GetYaxis()->SetLabelSize(0.05);
                h2B->GetYaxis()->SetTitle("Signal (V)");
                h2B->GetYaxis()->SetTitleSize(0.05);
                h2B->GetYaxis()->CenterTitle();
                gPad->SetTicks(1,1);
                gPad->Update();
                
                //c1->Print(cname);  // NB: this is the only way to see axis labels, otherwise they do not appear in the canvas on screen. Don't know why!
                
                // fill output bar tree
                for (int h=0; h<ntrg; h++) {
                    outbar.time = ttrg[h];
                    outbar.ievt = h; // should be modified, this is not relevant
                    outbar.bar = j; // current bar
                    outbar.pdA = vA[h];
                    outbar.pdB = vB[h];
                    outbar.pos = auxP*log(vB[h]/vA[h]);
                    outbar.E = auxE*sqrt(vB[h]*vA[h]);
                    toutbar->Fill();
                }
                
                // clean it up (after bar processing)!
                delete[] ttrg;
                delete[] tstop;
                delete[] vAB;
                delete[] vA;
                delete[] vB;
            }                     // end if(fAtot && fBtot)
            delete h2A;
            delete h2B;
            delete c1;
        }                  // bar j loop end
        
        // sort output bar tree by increasing time (for sorting a tree see https://root.cern.ch/root/roottalk/roottalk01/3646.html )
        int noutbar = toutbar->GetEntries();
        toutbar->Draw("time","","goff");
        int *tidx = new int[noutbar];
        TMath::Sort(noutbar, toutbar->GetV1(), tidx, 0); // sort increasing time
        TTree *toutbarS = (TTree*) toutbar->CloneTree(0); // output bar tree sorted with time
        for (int h=0; h<noutbar; h++) {
            toutbar->GetEntry(tidx[h]);
            toutbarS->Fill();
        }
        
        // implement DFE FPGA logic and build the resulting tree
        TTree *tout = data.CloneTree(0); // output tgf (after application of readout electronics model)
        tout->SetScanField(0);
        
        float tprev[30]; // array with previous trigger time. To be checked against conversion time and blanking time.
        for (int k=0; k<30; k++) tprev[k] = -1.;
        float tconv = 0.;
        int barid;
        
        barid=0 ; // m.galli added 29nov17 
        
        // init with first bar event
        toutbarS->GetEntry(0);
        inevt.time = outbar.time;
        inevt.mult = 1;
        inevt.Etot = outbar.E;
        for (int k=0; k<30; k++) inevt.bar[k] = -1.;
        inevt.bar[outbar.bar] = 0;
        inevt.pos[0] = outbar.pos;
        inevt.E[0] = outbar.E;
        tprev[outbar.bar] = outbar.time;
        
        // scan all bar events
        for (int h=1; h<noutbar; h++) {
            toutbarS->GetEntry(h);
            
            if (outbar.time - tprev[outbar.bar] > TMath::Max(tblank, tconv)) { // was bar free to trigger? (check conversion time and blanking time)
                
                if (outbar.time - inevt.time <= tjitter) { // join current bar in the same event
                    barid++;
                    inevt.mult += 1;
                    inevt.Etot += outbar.E;
                    inevt.bar[outbar.bar] = barid;
                    //inevt.pdA[barid] = outbar.pdA; // pdA and pdB are supposed to be in channels, while for outbar they are in V. Not relevant now.
                    //inevt.pdB[barid] = outbar.pdB;
                    inevt.pos[barid] = outbar.pos;
                    inevt.E[barid] = outbar.E;
                    
                    tprev[outbar.bar] = outbar.time; // update array with previous trigger time
                    
                } else {
                    
                    if (outbar.time - inevt.time <= tjitter + treject) {
                        continue; // bar was hit during blanking time window: skip trigger
                        
                    } else { // a new event is starting
                        
                        tout->Fill(); // save current event
                        tdet->Fill();
                        
                        tconv = (tlatency + 1. + 2.*(inevt.mult - 1))*tckADC; // conversion time for current event
                        if (verbose) printf("Event multiplicity = %d  -  total energy  %f  MeV  - conversion time  %f ms\n", inevt.mult, inevt.Etot, tconv);
                        
                        barid = 0; // init new event
                        inevt.time = outbar.time;
                        inevt.mult = 1;
                        inevt.Etot = outbar.E;
                        for (int k=0; k<30; k++) inevt.bar[k] = -1.;
                        inevt.bar[outbar.bar] = 0;
                        inevt.pos[0] = outbar.pos;
                        inevt.E[0] = outbar.E;
                        tprev[outbar.bar] = outbar.time;
                    }
                    
                }
            } else {
                continue; // bar was busy: skip trigger
            }
        }
        tout->Fill(); // save current event
        tdet->Fill();
        
        ncout = tout->GetEntries(); // detected number of counts
        
        if (ncout>1) {
            tout->Draw("time","","goff");
            outRMS = TMath::RMS(ncout, tout->GetV1());
        }
        
        // do some statistics on energy
        tgf->Draw("Etot","Etot>=0.5","goff");  // median energy above 0.5 MeV (to avoid detection threshold issue)
        inMedianE = TMath::Median(nc, tgf->GetV1());
        inMaxE = TMath::MaxElement(nc, tgf->GetV1());
        tout->Draw("Etot","Etot>=0.5","goff");
        outMedianE = TMath::Median(ncout, tout->GetV1());
        outMaxE = TMath::MaxElement(ncout, tout->GetV1());
        tsum->Fill();
        
        // final printouts for debug
        if (verbose) {
            cout << endl << "Initial bar array:" << endl;
            tbar->Scan("time:bar:E");
            cout << endl << "Final bar array:" << endl;
            //toutbar->Scan("time:bar:E");
            toutbarS->Scan("time:bar:E");
            cout << endl << "Final counts array:" << endl;
            tout->Scan("time:mult:Etot");
            
            for (int h=0; h<ncout; h++) tout->Show(h, 30);
        }
        printf("\nSummary: \nfluence = %f  -  duration = %f ms  -  n. physical hits = %d  -  n. detected counts = %d  -  initial time RMS = %f ms  -  final time RMS = %f ms \n", fluence, t50, nc, ncout, inRMS, outRMS);
        printf("Initial median energy (E>=0.5MeV) = %f  -  final median energy (E>=0.5MeV) = %f   -  n. initial max energy = %f  -  final max energy = %f\n", inMedianE, outMedianE, inMaxE, outMaxE);
        
        //----------------------
        
        
        // clean it up!   
        delete[] atcal;
        if (atsim) delete[] atsim;
        delete[] idc;
        delete tgf;
        delete tbar;
        delete toutbar;
        delete toutbarS;
        fcoll->Delete();
        
        // save data after each iteration
        tin->Write("", TObject::kOverwrite); // save only the new version of the tree
        tdet->Write("", TObject::kOverwrite);
        tsum->Write("", TObject::kOverwrite);
        //    fout.Write();
    }
    
    
    fout.Close();
    return;
}

//------------------------------------

void plotEvent(int i, double sigma, TTree *tbar) {
    
    char cname[100], cut[100], hName[100];
    
    // work on plane X
    sprintf(cname, "c%06d_X.eps", i);
    TCanvas *c1 = new TCanvas(cname, cname, 2400, 1200);
    c1->Divide(4,8);
    
    // loop on all bars and build scatter plots
    int irow, icol;
    for (int j=0; j<15; j++) {
        
        irow = j%8;
        icol = (j/8)*2;
        
        // plot pdA, bar j
        c1->cd(irow*4 + icol +1);
        sprintf(hName, "h2_%06d_bar%02d_A", i, j);
        TH2F *h2A = new TH2F(hName, hName, 1000, -3.*sigma, 3.*sigma, 200, 0., 200.);
        h2A->SetStats(0);
        h2A->GetXaxis()->SetLabelSize(0.1);
        h2A->GetYaxis()->SetLabelSize(0.1);
        h2A->Draw();
        
        sprintf(cut, "bar == %d", j);
        if (0) {
            cout << cut << endl;
            tbar->Scan("time:bar:pdA:pdB:pos:E", cut);
        }
        
        tbar->Draw("pdA:time", cut);
        if (gPad->GetPrimitive("Graph") != NULL) {
            TGraph *gA = (TGraph *) gPad->GetPrimitive("Graph");
            gA->SetMarkerStyle(5);
        }
        h2A->GetXaxis()->SetRangeUser(-3.*sigma, 3.*sigma);
        h2A->GetYaxis()->SetRangeUser(0.1, 100.);
        
        // plot pdB, bar j
        c1->cd(irow*4 + icol +2);
        sprintf(hName, "h2_%06d_bar%02d_B", i, j);
        TH2F *h2B = new TH2F(hName, hName, 1000, -3.*sigma, 3.*sigma, 200, 0., 200.);
        h2B->SetStats(0);
        h2B->GetXaxis()->SetLabelSize(0.1);
        h2B->GetYaxis()->SetLabelSize(0.1);
        h2B->Draw();
        sprintf(cut, "bar == %d", j);
        tbar->Draw("pdB:time", cut);
        if (gPad->GetPrimitive("Graph") != NULL) {
            TGraph *gB = (TGraph *) gPad->GetPrimitive("Graph");
            gB->SetMarkerStyle(5);
        }
        h2B->GetXaxis()->SetRangeUser(-3.*sigma, 3.*sigma);
        h2B->GetYaxis()->SetRangeUser(0.1, 100.);
    }
    //c1->Print(cname);
    delete c1;
    
    // work on plane Z
    sprintf(cname, "c%06d_Z.eps", i);
    TCanvas *c2 = new TCanvas(cname, cname, 2400, 1200);
    c2->Divide(4,8);
    
    // loop on all bars and build scatter plots
    for (int j=15; j<30; j++) {
        irow = (j-15)%8;
        icol = ((j-15)/8)*2;
        
        // plot pdA, bar j
        c2->cd(irow*4 + icol +1);
        sprintf(hName, "h2_%06d_bar%02d_A", i, j);
        TH2F *h2A = new TH2F(hName, hName, 1000, -3.*sigma, 3.*sigma, 200, 0., 200.);
        h2A->SetStats(0);
        h2A->GetXaxis()->SetLabelSize(0.1);
        h2A->GetYaxis()->SetLabelSize(0.1);
        h2A->Draw();
        
        sprintf(cut, "bar == %d", j);
        if (0) {
            cout << cut << endl;
            tbar->Scan("time:bar:pdA:pdB:pos:E", cut);
        }
        
        tbar->Draw("pdA:time", cut);
        if (gPad->GetPrimitive("Graph") != NULL) {
            TGraph *gA = (TGraph *) gPad->GetPrimitive("Graph");
            gA->SetMarkerStyle(5);
        }
        h2A->GetXaxis()->SetRangeUser(-3.*sigma, 3.*sigma);
        h2A->GetYaxis()->SetRangeUser(0.1, 100.);
        
        // plot pdB, bar j
        c2->cd(irow*4 + icol +2);
        sprintf(hName, "h2_%06d_bar%02d_B", i, j);
        TH2F *h2B = new TH2F(hName, hName, 1000, -3.*sigma, 3.*sigma, 200, 0., 200.);
        h2B->SetStats(0);
        h2B->GetXaxis()->SetLabelSize(0.1);
        h2B->GetYaxis()->SetLabelSize(0.1);
        h2B->Draw();
        sprintf(cut, "bar == %d", j);
        tbar->Draw("pdB:time", cut);
        if (gPad->GetPrimitive("Graph") != NULL) {
            TGraph *gB = (TGraph *) gPad->GetPrimitive("Graph");
            gB->SetMarkerStyle(5);
        }
        h2B->GetXaxis()->SetRangeUser(-3.*sigma, 3.*sigma);
        h2B->GetYaxis()->SetRangeUser(0.1, 100.);
    }
    //c2->Print(cname);
    delete c2;
    
}

// m.galli added 23aug17 revised 29nov17
// ------------------- main program --------------

int main( int  argc, char * argv[]) 
{ 
   std::cerr << "Beginning tgf_dead_time .... " << std::endl ; 
  
   int   knrdm = 2; 
   char *koutfilename="./sim.root"; 
   int   kidoffset    = 0  ; 
   float kminFluence  = 0.05; 
   float kmaxFluence  = 0.07; 
   float kminT50      = 0.02; 
   float kmaxT50      = 0.05; 
   int   kverbose     = 1; 
   int   ksourcetype  = 3; 
   float krate        = 20000.;
   char *ktimeprofile="./uniform_split10.root";

   if ( argc < 10 )
   {
     std::cerr << " ERROR: Missing arguments: " << std::endl ; 
     std::cerr << " Arguments are: "                << std::endl ;
     std::cerr << "  int   nrdm = 2 ;             " << std::endl ; 
     std::cerr << "  char *outfilename=./sim.root;" << std::endl ; 
     std::cerr << "  int   idoffset   = 0;    "     << std::endl ; 
     std::cerr << "  float minFluence = 0.05; "     << std::endl ; 
     std::cerr << "  float maxFluence = 0.07; "     << std::endl ; 
     std::cerr << "  float minT50     = 0.02; "     << std::endl ; 
     std::cerr << "  float maxT50     = 0.05; "     << std::endl ; 
     std::cerr << "  int   verbose    = 0;    "     << std::endl ; 
     std::cerr << "  int   sourcetype = 3;    "     << std::endl ; 
     std::cerr << "  float rate       = 20000.; "   << std::endl ;
     std::cerr << "  char *timeprofile=./uniform_split10.root;" << std::endl ; 
     
     exit(1) ; 
   }      
   knrdm        = atoi(argv[1])  ; 
   koutfilename = argv[2]        ; 
   kidoffset    = atoi(argv[3])  ; 
   kminFluence  = atof(argv[4])  ; 
   kmaxFluence  = atof(argv[5])  ; 
   kminT50      = atof(argv[6])  ;
   kmaxT50      = atof(argv[7])  ;
   kverbose     = atoi(argv[8])  ;
   ksourcetype  = atoi(argv[9])  ;
   krate        = atof(argv[10]) ;
   ktimeprofile = argv[11]        ; 
           
   std::cerr << "Run Parameters: "                << std::endl 
               << " nrdm:        " << knrdm       << std::endl  
               << " outfilename: " << koutfilename<< std::endl  
               << " idoffset:    " << kidoffset   << std::endl  
               << " minFluence:  " << kminFluence << std::endl  
               << " maxFluence:  " << kmaxFluence << std::endl  
               << " minT50:      " << kminT50     << std::endl  
               << " maxT50:      " << kmaxT50     << std::endl               
               << " verbose:     " << kverbose    << std::endl     
               << " sourcetype:  " << ksourcetype << std::endl     
               << " rate:        " << krate       << std::endl 
               << " timeprofile: " << ktimeprofile << std::endl ;
   
   tgf_dead_time_2_1( 
       knrdm, 
       koutfilename,kidoffset,
       kminFluence,kmaxFluence,
       kminT50,kmaxT50,
       kverbose,ksourcetype,
       krate,ktimeprofile) ; 
   
      
   std::cerr << "ended tgf_dead_time .... " << std::endl ; 
}



