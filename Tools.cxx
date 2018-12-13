#include "vector"
#include "Tools.h"
#include "Layer.h"
#include "Hit.h"
#include "Particle.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TClonesArray.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TString.h"
#include "TPaveText.h"
#include "THStack.h"
#include "Riostream.h"

//using namespace TMath;

void verbosities(bool b_verbose, bool b_multiscatter, bool b_noise, int kExp){//general information about the simulation

  printf("\n\n+++ TRACKS START - tracks generation +++\n\n");

  if (b_verbose) {
    printf("Printing vertex and hit coordinates: ON\n\n");
  }else{printf("Printing vertex and hit coordinates: OFF - only ten event displayed\n\n");}

  if (b_multiscatter) {
    printf("Applying multiple scattering: ON\n\n");
  }else{printf("Applying multiple scattering: OFF\n\n");}

  if (b_noise) {
    printf("Applying noise: ON\n\nAll distances are in cm, all angles are in rad.\n\n\n");
  }else{printf("Applying noise: OFF\n\nAll distances are in cm, all angles are in rad.\n\n\n");}

  printf("Ready to generate %d events! The detection info will be stored in a .root file containing a TTree with two branches, respectively the hits' coordinates on the first and the second Si layer.\n\n\n",kExp);

}

void graphstyler(TGraph &graph, int divide){//TGraph makeup
  if(divide==1){
    graph.GetXaxis()->SetTitleSize(0.04);
    graph.GetYaxis()->SetTitleSize(0.035);
    graph.GetXaxis()->SetTitleOffset(1.1);
  }
  if(divide==2){
    graph.GetXaxis()->SetLabelSize(0.04);
    graph.GetYaxis()->SetLabelSize(0.04);
    graph.GetXaxis()->SetTitleSize(0.04);
    graph.GetYaxis()->SetTitleSize(0.04);
    graph.GetXaxis()->SetTitleOffset(1.15);
    graph.GetYaxis()->SetTitleOffset(1.3);
  }
  if(divide==4){
    graph.GetXaxis()->SetLabelSize(0.045);
    graph.GetYaxis()->SetLabelSize(0.045);
    graph.GetXaxis()->SetTitleSize(0.05);
    graph.GetYaxis()->SetTitleSize(0.05);
    graph.GetXaxis()->SetTitleOffset(0.9);
    graph.GetYaxis()->SetTitleOffset(0.7);
  }
}

void histostyler(TH1 &histo, int divide){//TH1 makeup
  if(divide==1){
    histo.GetXaxis()->SetTitleSize(0.04);
    histo.GetYaxis()->SetTitleSize(0.035);
    histo.GetXaxis()->SetTitleOffset(1.1);
  }
  if(divide==2){
    histo.GetXaxis()->SetLabelSize(0.04);
    histo.GetYaxis()->SetLabelSize(0.04);
    histo.GetXaxis()->SetTitleSize(0.04);
    histo.GetYaxis()->SetTitleSize(0.04);
    histo.GetXaxis()->SetTitleOffset(1.15);
    histo.GetYaxis()->SetTitleOffset(1.3);
  }
  if(divide==4){
    histo.GetXaxis()->SetLabelSize(0.045);
    histo.GetYaxis()->SetLabelSize(0.045);
    histo.GetXaxis()->SetTitleSize(0.05);
    histo.GetYaxis()->SetTitleSize(0.05);
    histo.GetXaxis()->SetTitleOffset(0.9);
    histo.GetYaxis()->SetTitleOffset(0.7);
  }
}

void stackstyler(THStack &stack){//THStack makeup

    stack.GetXaxis()->SetLabelSize(0.045);
    stack.GetYaxis()->SetLabelSize(0.045);
    stack.GetXaxis()->SetTitleSize(0.05);
    stack.GetYaxis()->SetTitleSize(0.05);
    stack.GetXaxis()->SetTitleOffset(0.9);
    stack.GetYaxis()->SetTitleOffset(0.7);
}

void pavestyler(TPaveText &pave, double textsize){//TPave makeup
  pave.SetTextAlign(13);
  pave.SetFillStyle(0);
  pave.SetShadowColor(0);
  pave.SetLineColor(0);
  pave.SetBorderSize(0);
  pave.SetMargin(0);
  pave.SetTextSize(textsize);
}

void MSaveBigPNG(TString filename, double scale) {//to scale up png resolution
    TCanvas* old_canv = gPad->GetCanvas();

    gROOT->SetBatch(kTRUE);
    gROOT->ForceStyle(kTRUE);

    Int_t orig_msz = gStyle->GetMarkerSize();
    Int_t orig_mst = gStyle->GetMarkerStyle();
    Int_t orig_lt  = gStyle->GetLineWidth();

    gStyle->SetMarkerSize(1.0+scale/5);
    gStyle->SetMarkerStyle(20);
    gStyle->SetLineWidth(orig_lt*scale);

    if(filename == "zio") {
        filename = old_canv->GetName();
        filename += ".png";
    }

    Int_t old_width  = old_canv->GetWindowWidth();
    Int_t old_height = old_canv->GetWindowHeight();

    Int_t new_width = old_width * scale;
    Int_t new_height= old_height* scale;

    TCanvas* temp_canvas = new TCanvas("temp", "", new_width, new_height);
    old_canv->DrawClonePad();

    temp_canvas->Draw();
    temp_canvas->SaveAs(filename);
    temp_canvas->Close();

    gStyle->SetMarkerSize(orig_msz);
    gStyle->SetMarkerStyle(orig_mst);
    gStyle->SetLineWidth(orig_lt);

    gROOT->ForceStyle(kFALSE);
    gROOT->SetBatch(kFALSE);

    return;
}

double *hit_point(double x0, double y0, double z0, double theta, double phi, double R) {//evaluate hit coordinates

  static double hit[3];
  double c1=TMath::Sin(theta)*TMath::Cos(phi), c2=TMath::Sin(theta)*TMath::Sin(phi), c3=TMath::Cos(theta);//direction cosines
  double delta=2*x0*y0*c1*c2-c1*c1*y0*y0+c1*c1*R*R-c2*c2*x0*x0+c2*c2*R*R;//delta of II degree equation (>= 0 by construction)
  double t_p=(-(x0*c1-y0*c2)+TMath::Sqrt(delta))/(c1*c1+c2*c2);//+ solution
  double t_m=(-(x0*c1-y0*c2)-TMath::Sqrt(delta))/(c1*c1+c2*c2);//- solution

  //calculate the intersection coordinates (x,y,z)
  if (t_p >= 0) {
    hit[0] = x0 + c1*t_p;
    hit[1] = y0 + c2*t_p;
    hit[2] = z0 + c3*t_p;
  }

  else {
    hit[0] = x0 + c1*t_m;
    hit[1] = y0 + c2*t_m;
    hit[2] = z0 + c3*t_m;
  }

  return hit;//return a pointer to the first element of hit[]
}

bool detect(Hit* vtx, Layer* L, Particle &part, TClonesArray &cross, bool b_verbose, bool b_multiscatter, int &counter, TH1D** histo){

  double *hit_buffer;
  bool b_cross=false;

  hit_buffer=hit_point(vtx->GetX(),vtx->GetY(),vtx->GetZ(),part.GetTheta(),part.GetPhi(),L->GetRadius());//evaluate hit coordinates

  if(*(hit_buffer+2) >= -(L->GetWidth()/2.) && *(hit_buffer+2) <= (L->GetWidth()/2.)) {

    b_cross = true;//yes we have detection
    //gSystem->Beep(440,1);

    new(cross[counter])Hit(*(hit_buffer+0),*(hit_buffer+1),*(hit_buffer+2));//fill the TCA with Hit object

    for(int i=0;i<=2;i++){histo[i]->Fill(*(hit_buffer+i));}//fill the histogram with hist coordinates

    if (b_cross&&b_multiscatter) {
      part.Rotate(L->GetRMS());//if multiscattering ON set new angles
    }

    if (b_verbose) {
      cout<<"Hit with "<<L->GetLayerName();
      printf(" at (%f, %f, %f)\nAngles after: theta %f - phi %f\n\n",((Hit*)cross[counter])->GetX(),((Hit*)cross[counter])->GetY(),((Hit*)cross[counter])->GetZ(),part.GetTheta(),part.GetPhi());
    }

    counter++;//this counter will be taken by "detect" during the next call

  }else{
    if(b_verbose){
      cout<<"Does not hit "<<L->GetLayerName()<<endl<<endl;
    }
  }

  return b_cross;

}

void noise(bool b_verbose, int Noise, int Mult, TClonesArray &cross, Layer* L){//noise hit generator

  int index_noise=Mult;

  for(int i=0; i<Noise; i++){
    new(cross[index_noise])Hit(L->GetRadius(),L->GetWidth());//fill the TCA with a random noise hit
    if(b_verbose){
      cout<<"> Noise hit with "<<L->GetLayerName();
      printf(" at (%f, %f,%f) <\n\n",((Hit*)cross[index_noise])->GetX(),((Hit*)cross[index_noise])->GetY(),((Hit*)cross[index_noise])->GetZ());
    }
    index_noise++;
  }

}

void smear(int index, double sigmaz, double sigmarf, double R, TClonesArray &cross){//gaussian smearing

  double dphi=(gRandom->Gaus(0,sigmarf))/R;
  double dz=gRandom->Gaus(0,sigmaz);

  Hit *hit_buffer=(Hit*)cross.At(index);//alias to retrieve the hit coordinates

  double x=hit_buffer->GetX();
  double y=hit_buffer->GetY();
  double z=hit_buffer->GetZ();

  double phi=TMath::ACos(x/R);
  double theta=TMath::ACos(z/TMath::Sqrt(x*x+y*y+z*z));//angles before the smearing

  if (gRandom->Rndm()<0.5) {
    phi+=dphi;
    x=R*TMath::Cos(phi);
    y=R*TMath::Sin(phi);
    z+=dz;
    theta=TMath::ACos(z/TMath::Sqrt(x*x+y*y+z*z));
  } else {
    phi+=dphi;
    x=R*TMath::Cos(phi);
    y=R*TMath::Sin(phi);
    z-=dz;
    theta=TMath::ACos(z/TMath::Sqrt(x*x+y*y+z*z));
  }//new angles and coordinates after gaussian smearing

  hit_buffer->SetX(x);
  hit_buffer->SetY(y);
  hit_buffer->SetZ(z);

}

bool peakfinder(TH1D* histo, double ampli, int width){//tracklet ambiguity check

  //the idea of this check is to evaluate whether
  //the tracklet histogram contains more than 1 region
  //with a certain amount of counts - in this case
  //we cannot say which region could contain the true vertex
  //and discard the event, otherwise we accept it
  //and evaluate the reconstructed z

  //ampli is the scale factor that determines the "amount of counts"
  //width in the number of bins in which we look for it

  int fisrtbin=(int)(ampli/2+0.5);
  int lastbin=(int)(ampli/2+0.5)-1;//define the check range
  bool peakit=false;
  int kBin=histo->GetSize()-2;//number of bin excluding under- and overflow
  int binC=histo->GetMaximumBin();//bin index of the first global maximum
  double Max=histo->GetBinContent(binC);//value of the first global maximum
  double ClusterSize=0;

  for(int i=fisrtbin;i<kBin-lastbin;i++){//loop over the check range
    for(int j=-lastbin;j<=lastbin;j++){//loop over width
      ClusterSize+=histo->GetBinContent(i-j);
    }//end loop over width
    if(ClusterSize>=ampli*Max&&((i<=binC-width)||(i>=binC+width))){//this is the case in which we have ambiguity
      peakit=false;
      ClusterSize=0;
      break;
    }else{//this is the case in which we don't have ambiguity - the event is good
      ClusterSize=0;
      peakit=true;
    }
  }//end loop over the check range

  return peakit;
}
