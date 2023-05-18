from ROOT import gInterpreter
def def_cpp():
        return gInterpreter.Declare('''
#include <math.h>
#include <TLorentzVector.h>
#include <algorithm> 

TH2F *bTagEff = 0;

using namespace ROOT::VecOps;
        
        
float CalcMinPhi(float phi1, float phi2){
   float diff = fmod((phi1-phi2 +3.1415), (2 *3.1415)) - 3.1415;
   if (diff < -3.1415) {
       diff = diff + 3.1415*2;
   }
   return diff;
}

float CalcDeltaR(float eta1, float  eta2, float  phi1, float  phi2) {
    return pow(pow(eta1-eta2, 2) + pow(CalcMinPhi(phi1, phi2), 2), .5);
}        
        
bool hard_process(int value){
    return ((value & (1 << 7)) == (1 << 7)) || ((value & (1 << 8)) == (1 << 8));
}      
        

int GetGenMultiplicity(const RVec<int> &GenPart_statusFlags, const RVec<int> &GenPart_pdgId,
                    const RVec<float> &GenPart_pt, const RVec<float> &GenPart_eta, const RVec<float> &GenPart_phi){
    int multiplicity = 0;
    int nGen = GenPart_pdgId.size();
    bool _ishard;
    float _pt;
    float _eta;
    float _phi;
    float _pdgid;
    
    //find initial state information
    int inBs=0;
    int inSs=0;
    int inOQs=0;
    int inGs=0;
    for (unsigned i = 0; i < nGen; ++i) {
            _ishard = hard_process(GenPart_statusFlags.at(i));
            _pt = GenPart_pt.at(i);
            _pdgid = GenPart_pdgId.at(i);
            //only take hard process particles
            if (!_ishard) continue;
            // only initial state has pt of 0
            if (!(_pt==0)) continue;
            if (abs(_pdgid)==5) inBs+=1;
            if (abs(_pdgid)==3) inSs+=1;
            if ((abs(_pdgid)<=9) && (abs(_pdgid)!=3) && (abs(_pdgid)!=5)) inOQs+=1;
            if (abs(_pdgid)==21) inGs+=1; 
    }
    
    //final state information
    //leading quark stats
    float leadOutB = 0;
    float leadOutB_id = 0;
    float leadOutB_eta = 0;
    float leadOutB_phi = 0;
    float leadOutBMinus = 0;
    float leadOutBPlus = 0;
    float leadOutS = 0;
    for (unsigned i = 0; i < nGen; ++i) {
        _ishard = hard_process(GenPart_statusFlags.at(i));
        _pt = GenPart_pt.at(i);
        _eta = GenPart_eta.at(i);
        _phi = GenPart_phi.at(i);
        _pdgid = GenPart_pdgId.at(i);
        
        //only take hard process particles
        if (!_ishard) continue;
        if (_pt==0) continue;
        if ((abs(_pdgid)==5) && (_pt > leadOutB)){
            leadOutB = _pt;
            leadOutB_id = _pdgid;
            leadOutB_eta = _eta;
            leadOutB_phi = _phi;
        }
        if ((_pdgid==5) && (_pt > leadOutBPlus)) leadOutBPlus = _pt;
        if ((_pdgid==-5) && (_pt > leadOutBMinus)) leadOutBMinus = _pt;
        if ((abs(_pdgid)==3) && (_pt > leadOutS)) leadOutS = _pt;
    }
    //subleading b quark stats
    float subLeadOutB = 0;
    float deltaR = 0;
    for (unsigned i = 0; i < nGen; ++i) {
        _ishard = hard_process(GenPart_statusFlags.at(i));
        _pt = GenPart_pt.at(i);
        _eta = GenPart_eta.at(i);
        _phi = GenPart_phi.at(i);
        _pdgid = GenPart_pdgId.at(i);
        //only take hard process particles
        if (!_ishard) continue;
        if (_pt==0) continue;
        //only b
        if (abs(_pdgid)!=5) continue;
        deltaR = CalcDeltaR(_eta, leadOutB_eta, _phi, leadOutB_phi);
        if ((deltaR < .4) && (_pdgid == leadOutB_id)) continue;
        if ((_pt < leadOutB) && (_pt > subLeadOutB)) subLeadOutB = _pt;
    }
        
    // multiplicity logic, based on combined initial state and final state
    // 0j: 0,1, 1b: 2,3,4, 1s: 5,6, 1b+1s: 7,8,9,10, 2b: 11,12,13,14
    //5 5
    if ((inBs == 2) && (inSs == 0) && (inOQs == 0) && (inGs == 0)){
        //0b
        multiplicity += pow(2,0);
        //2b'
        if ((leadOutB>0) && (subLeadOutB>0)){ 
            multiplicity += pow(2,11);
        //1b1s'
        } else if ((leadOutB>0) && (leadOutS>0)) multiplicity += pow(2,7);
    }
    //5 3 
    if ((inBs == 1) && (inSs == 1) && (inOQs == 0) && (inGs == 0)) {
        //0b
        multiplicity += pow(2,0);
        //2b'
        if ((leadOutB>0) && (subLeadOutB>0)) multiplicity += pow(2,12);
    }
    //5 21
    if ((inBs == 1) && (inSs == 0) && (inOQs == 0) && (inGs == 1)){
    
        //1b
        if (leadOutB > leadOutS) multiplicity += pow(2,2);
        //1s
        if (leadOutS > leadOutB) multiplicity += pow(2,5);
    }
    //3 21
    if ((inBs == 0) && (inSs == 1) && (inOQs == 0) && (inGs == 1)) multiplicity += pow(2,3);
    //21 21
    if ((inBs == 0) && (inSs == 0) && (inOQs == 0) && (inGs == 2)){
        //2b
        if ((leadOutBPlus > leadOutS) && (leadOutBMinus > leadOutS)){
            multiplicity += pow(2,13);
        //1b+1s
        } else if ((leadOutS > 0) && (leadOutB > 0)) multiplicity += pow(2,8);
    }
    //5 + x        
    if ((inBs == 1) && (inSs == 0) && (inOQs == 1) && (inGs ==0)){
        //1b+1s'
        if ((leadOutB>0) && (leadOutS>0)) {
            multiplicity += pow(2,9);
        //1b
        } else if (leadOutB>0) multiplicity += pow(2,4);
        //1s
        if (leadOutS>0) multiplicity += pow(2,6);
    }
    //3 + x        
    if ((inBs == 0) && (inSs == 1) && (inOQs == 1) && (inGs ==0)){
        //1b
        multiplicity += pow(2,10);
    }
    //2 q       
    if ((inBs == 0) && (inSs == 0) && (inOQs == 2) && (inGs ==0)){
        //2b
        multiplicity += pow(2,14);
    }

    return multiplicity;
}
       
        
        
      
        
        
RVec<float> DivideVec(const RVec<float> &one, const RVec<float> &two){
    RVec<float> divided = {};
    for (unsigned j = 0; j < one.size(); ++j) {
        divided.push_back(one[j]/two[j]);
    }
    return 1.0-divided;
}
         

        
int GetNbJets(const RVec<int> &GenPart_statusFlags, const RVec<int> &GenPart_pdgId){
    int nBJets =  0;
    std::cout << "-----" << std::endl;
    for (unsigned i = 0; i < GenPart_statusFlags.size(); ++i) {
    
        nBJets += int(hard_process(GenPart_statusFlags.at(i)) && (abs(GenPart_pdgId.at(i))==5));

    }
    return nBJets;
}

int GetNsJets(const RVec<int> &GenPart_statusFlags, const RVec<int> &GenPart_pdgId){
    int nSJets =  0;
    for (unsigned i = 0; i < GenPart_statusFlags.size(); ++i) {
        nSJets += int(hard_process(GenPart_statusFlags.at(i)) && (abs(GenPart_pdgId.at(i))==3));
    }
    return nSJets;
}


RVec<float> GetPDFWeight(const RVec<float> &LHEPdfWeight){
    // based on https://arxiv.org/pdf/2203.05506.pdf
    auto LHEPdfWeight_sorted = Sort(LHEPdfWeight);
    ///take 68% envelope
    int nVariations = LHEPdfWeight_sorted.size();
    float percentOffEnds = (1.-.68)/2.;
    int numberOffEnds = floor(percentOffEnds*(float)nVariations);
    float delta = 0;
    if (nVariations > 0){
        float up = LHEPdfWeight_sorted.at(nVariations-numberOffEnds);
        float down = LHEPdfWeight_sorted.at(numberOffEnds);
        delta = (float)abs(up-down)/2;
    }
    RVec<float> weight={(float)1.-delta,(float)1.+delta};
    return weight;
}     
       
        
 RVec<float> GetScaleUncertainty( const RVec<float> &LHEScaleWeight){
          RVec<float> weight={1.,1.};
       
          //Scale weight usage based on conjecture that this works by choosing the maximal variation in each direction for renormalization and factorization, independently. Contents as of https://cms-nanoaod-integration.web.cern.ch/integration/master/mc94X_doc.html#LHE
          float RenUncertaintyUp=99;
          float RenUncertaintyDown=-99;
          float FactUncertaintyUp=99;
          float FactUncertaintyDown=-99;
          if(LHEScaleWeight.size()==44){ 
              RenUncertaintyUp=TMath::Max(LHEScaleWeight[9]-1., TMath::Max(LHEScaleWeight[23]-1., TMath::Max(LHEScaleWeight[38]-1., 0.)));
              RenUncertaintyDown=TMath::Max(1.-LHEScaleWeight[9], TMath::Max(1.-LHEScaleWeight[23], TMath::Max(1.-LHEScaleWeight[38], 0.)));
              FactUncertaintyUp=TMath::Max(LHEScaleWeight[19]-1., TMath::Max(LHEScaleWeight[23]-1., TMath::Max(LHEScaleWeight[28]-1., 0.)));
              FactUncertaintyDown=TMath::Max(1.-LHEScaleWeight[19], TMath::Max(1.-LHEScaleWeight[23], TMath::Max(1.-LHEScaleWeight[28], 0.)));
          } else if (LHEScaleWeight.size() == 9) {
              RenUncertaintyUp=TMath::Max(LHEScaleWeight[1]-1., TMath::Max(LHEScaleWeight[4]-1., TMath::Max(LHEScaleWeight[7]-1., 0.)));
              RenUncertaintyDown=TMath::Max(1.-LHEScaleWeight[1], TMath::Max(1.-LHEScaleWeight[4], TMath::Max(1.-LHEScaleWeight[7], 0.)));
              FactUncertaintyUp=TMath::Max(LHEScaleWeight[3]-1., TMath::Max(LHEScaleWeight[4]-1., TMath::Max(LHEScaleWeight[5]-1., 0.)));
              FactUncertaintyDown=TMath::Max(1.-LHEScaleWeight[3], TMath::Max(1.-LHEScaleWeight[4], TMath::Max(1.-LHEScaleWeight[5], 0.)));
          }
          //according to Eq. 8 or 10 in https://arxiv.org/pdf/1101.0536.pdf, modified by Eq. 6 to a factor of 1/N_rep, which given each variation is up and down of one representation, should be 51-1
          weight[0]+=TMath::Sqrt(pow(RenUncertaintyUp  ,2)+pow(FactUncertaintyUp  ,2));
          weight[1]-=TMath::Sqrt(pow(RenUncertaintyDown,2)+pow(FactUncertaintyDown,2));
          return weight;
        }        
        
        
        

        
        RVec<float> min_delta_r(const RVec<float> eta1, const RVec<float> eta2, const RVec<float> phi1, const RVec<float> phi2){
            TLorentzVector v1; 
            TLorentzVector v2;
            RVec<float> dr_vec = {};
            for (unsigned i = 0; i < eta1.size(); ++i) {
                v1.SetPtEtaPhiM(1, eta1[i], phi1[i], 1);
                double min_dr = 999;
                for (unsigned j = 0; j < eta2.size(); ++j) {
                    v2.SetPtEtaPhiM(1, eta2[j], phi2[j], 1);
                    min_dr = min(min_dr, v2.DeltaR(v1));
                }
                dr_vec.push_back(min_dr);
            }
            return dr_vec;
        }
        
        float min_vec(const RVec<float> x){
            float min_val = 999999;
            for (unsigned i = 0; i < x.size(); ++i) {
                min_val = min(min_val, x[i]);
            }
            return min_val;
        }
        float max_vec(const RVec<float> x){
            float max_vec = -999999;
            for (unsigned i = 0; i < x.size(); ++i) {
                max_vec = max(max_vec, x[i]);
            }
            return max_vec;
        }
        float CalcAverage(const RVec<float> x){
            float sum = 0;
            for (unsigned i = 0; i < x.size(); ++i) {
                sum+=x[i];
            }
            return sum/float(x.size());
        }

        float GetBTagWeight(const RVec<int> &isGoodBJet, const RVec<int> &isGoodJet, const RVec<int> &JetHadronFlav, const RVec<float> &JetPt, const RVec<float> &bTagSF) {
          float weight = 1.;
          for (unsigned j = 0; j < isGoodBJet.size(); ++j) {
            if (isGoodBJet[j]) {
              weight *= bTagSF[j];
            }  else if (isGoodJet[j]) {
            int HadronFlav = 2.;
            if (JetHadronFlav[j] == 5)
              HadronFlav=0.;
            else if (JetHadronFlav[j] == 4)
              HadronFlav=1.;
              int effBin = bTagEff->FindBin(HadronFlav, JetPt[j]);
              float eff  = bTagEff->GetBinContent(effBin);
              weight *= (1. - eff * bTagSF[j]) / (1. - eff);
            }
          }
          return (float) weight;
        }

        RVec<float> GetBTagWeightPerJet(const RVec<int> &isGoodBJet, const RVec<int> &isGoodJet, const RVec<int> &JetHadronFlav, const RVec<float> &JetPt, const RVec<float> &bTagSF) {
           RVec<float> weights = {};
          for (unsigned j = 0; j < isGoodBJet.size(); ++j) {
            if (isGoodBJet[j]) {
              weights.push_back( bTagSF[j]);
            }  else if (isGoodJet[j]) {
            int HadronFlav = 2.;
            if (JetHadronFlav[j] == 5)
              HadronFlav=0.;
            else if (JetHadronFlav[j] == 4)
              HadronFlav=1.;
              int effBin = bTagEff->FindBin(HadronFlav, JetPt[j]);
              float eff  = bTagEff->GetBinContent(effBin);
              weights.push_back(  (1. - eff * bTagSF[j]) / (1. - eff) );
            }
          }
          return weights;
        }
        


        TH2F *PUIDSF_true  = 0;
        TH2F *PUIDSF_false = 0;
        TH2F *PUIDUnc_true = 0;
        TH2F *PUIDUnc_false= 0;
        TH2F *PUIDEff_true = 0;
        TH2F *PUIDEff_false= 0;
        
        float clip(float x){
            return (float)std::min((float)x,(float)5.);
        }
        RVec<float> GetPUIDweight(const RVec<float> &JetPt, const RVec<float> &JetEta, const RVec<int> &JetGenJetIdx, const RVec<int> &PUID,  float verbose) {
          if (verbose>0) std::cout << "------" << std::endl;
          RVec<float> weights = {1.,1.,1.};
          for(unsigned j = 0; j < JetGenJetIdx.size(); ++j) {
            if (verbose>0) std::cout << JetPt[j] << " " <<  JetEta[j] << std::endl;
            const int SFbin = PUIDSF_true->FindBin(JetPt[j], JetEta[j]);
            if(PUID[j] & 1){
              if(JetGenJetIdx[j]>=0){
                const float SF = PUIDSF_true->GetBinContent(SFbin);
                const float unc= PUIDUnc_true->GetBinContent(SFbin);
                weights[0] *= clip(SF);
                weights[1] *= clip(SF+unc);
                weights[2] *= clip(SF-unc);
                if (verbose>0) std::cout << "PUID gen " <<  SF << " " << unc << std::endl;
              }
              else{
                const float SF = PUIDSF_false->GetBinContent(SFbin);
                const float unc= PUIDUnc_false->GetBinContent(SFbin);
                weights[0] *= clip(SF);
                weights[1] *= clip(SF+unc);
                weights[2] *= clip(SF-unc);
                if (verbose>0) std::cout << "PUID no gen " <<  SF << " " << unc << std::endl;
              }
            }
            else{
              if(JetGenJetIdx[j]>=0){
                const float SF = PUIDSF_true->GetBinContent(SFbin);
                const float unc= PUIDUnc_true->GetBinContent(SFbin);
                const float eff= PUIDEff_true->GetBinContent(SFbin);
                weights[0] *= clip((1-eff*SF)/(1-eff));
                weights[1] *= clip((1-eff*(SF+unc))/(1-eff));
                weights[2] *= clip((1-eff*(SF-unc))/(1-eff));
                if (verbose>0) std::cout << "no PUID gen " <<  SF << " " << unc << " " << eff << std::endl;
              }
              else{
                const float SF = PUIDSF_false->GetBinContent(SFbin);
                const float unc= PUIDUnc_false->GetBinContent(SFbin);
                const float eff= PUIDEff_false->GetBinContent(SFbin);
                weights[0] *= clip((1-eff*SF)/(1-eff));
                weights[1] *= clip((1-eff*(SF+unc))/(1-eff));
                weights[2] *= clip((1-eff*(SF-unc))/(1-eff));
                if (verbose>0) std::cout << "no PUID no gen " <<  SF << " " << unc << " " << eff << std::endl;
              }
            }
          }
          if (verbose>0) std::cout <<"weights: " <<  weights[0] << " " << weights[1] << " " << weights[2] << " " << std::endl;

            if (verbose>0) std::cout <<"ratio: " << (weights[1]-weights[2])/weights[0] << std::endl;
            return weights;
        }

        RVec<float> GetPUIDweightPerJet(const RVec<float> &JetPt, const RVec<float> &JetEta, const RVec<int> &JetGenJetIdx, const RVec<int> &PUID,  float verbose) {
          if (verbose>0) std::cout << "------" << std::endl;
          RVec<float> weights = {};
          for(unsigned j = 0; j < JetGenJetIdx.size(); ++j) {
            if (verbose>0) std::cout << JetPt[j] << " " <<  JetEta[j] << std::endl;
            const int SFbin = PUIDSF_true->FindBin(JetPt[j], JetEta[j]);
            if(PUID[j] & 1){
              if(JetGenJetIdx[j]>=0){
                const float SF = PUIDSF_true->GetBinContent(SFbin);
                const float unc= PUIDUnc_true->GetBinContent(SFbin);
                weights.push_back(clip(SF+unc)-clip(SF));
              }
              else{
                const float SF = PUIDSF_false->GetBinContent(SFbin);
                const float unc= PUIDUnc_false->GetBinContent(SFbin);
                weights.push_back(clip(SF+unc)-clip(SF));
              }
            }
            else{
              if(JetGenJetIdx[j]>=0){
                const float SF = PUIDSF_true->GetBinContent(SFbin);
                const float unc= PUIDUnc_true->GetBinContent(SFbin);
                const float eff= PUIDEff_true->GetBinContent(SFbin);
                weights.push_back( clip((1-eff*(SF+unc))/(1-eff)) - clip((1-eff*SF)/(1-eff)));
              }
              else{
                const float SF = PUIDSF_false->GetBinContent(SFbin);
                const float unc= PUIDUnc_false->GetBinContent(SFbin);
                const float eff= PUIDEff_false->GetBinContent(SFbin);
                weights.push_back( clip((1-eff*(SF+unc))/(1-eff)) - clip((1-eff*SF)/(1-eff)));
              }
            }
          }
            return weights;
        }

        
        
        
        RVec<float> GetPDFandScaleUncertainty(const RVec<float> &LHEPdfWeight, const RVec<float> &LHEScaleWeight){
          RVec<float> weight={1.,1.};
          if(LHEPdfWeight.size()==0 || LHEScaleWeight.size()==0) return weight;
          //according to Eq. 3 https://arxiv.org/pdf/1101.0536.pdf, 103 members expected according to https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopSystematics#PDF, 102 variation members stored, relative to nominal (1)
          float PDFuncertaintyUp=0.;
          float PDFuncertaintyDown=0.;
          for(unsigned i=0; i<LHEPdfWeight.size()/2; ++i){
            PDFuncertaintyUp  += pow(TMath::Max(TMath::Max(LHEPdfWeight[2*i]-1.,LHEPdfWeight[2*i+1]-1.),0.),2);
            PDFuncertaintyDown+= pow(TMath::Max(TMath::Max(1.-LHEPdfWeight[2*i],1.-LHEPdfWeight[2*i+1]),0.),2);
          }
          //Scale weight usage based on conjecture that this works by choosing the maximal variation in each direction for renormalization and factorization, independently. Contents as of https://cms-nanoaod-integration.web.cern.ch/integration/master/mc94X_doc.html#LHE
          float RenUncertaintyUp=99;
          float RenUncertaintyDown=-99;
          float FactUncertaintyUp=99;
          float FactUncertaintyDown=-99;
          if(LHEScaleWeight.size()==44){ 
              RenUncertaintyUp=TMath::Max(LHEScaleWeight[9]-1., TMath::Max(LHEScaleWeight[23]-1., TMath::Max(LHEScaleWeight[38]-1., 0.)));
              RenUncertaintyDown=TMath::Max(1.-LHEScaleWeight[9], TMath::Max(1.-LHEScaleWeight[23], TMath::Max(1.-LHEScaleWeight[38], 0.)));
              FactUncertaintyUp=TMath::Max(LHEScaleWeight[19]-1., TMath::Max(LHEScaleWeight[23]-1., TMath::Max(LHEScaleWeight[28]-1., 0.)));
              FactUncertaintyDown=TMath::Max(1.-LHEScaleWeight[19], TMath::Max(1.-LHEScaleWeight[23], TMath::Max(1.-LHEScaleWeight[28], 0.)));
          } else if (LHEScaleWeight.size() == 9) {
              RenUncertaintyUp=TMath::Max(LHEScaleWeight[1]-1., TMath::Max(LHEScaleWeight[4]-1., TMath::Max(LHEScaleWeight[7]-1., 0.)));
              RenUncertaintyDown=TMath::Max(1.-LHEScaleWeight[1], TMath::Max(1.-LHEScaleWeight[4], TMath::Max(1.-LHEScaleWeight[7], 0.)));
              FactUncertaintyUp=TMath::Max(LHEScaleWeight[3]-1., TMath::Max(LHEScaleWeight[4]-1., TMath::Max(LHEScaleWeight[5]-1., 0.)));
              FactUncertaintyDown=TMath::Max(1.-LHEScaleWeight[3], TMath::Max(1.-LHEScaleWeight[4], TMath::Max(1.-LHEScaleWeight[5], 0.)));
          }
          //according to Eq. 8 or 10 in https://arxiv.org/pdf/1101.0536.pdf, modified by Eq. 6 to a factor of 1/N_rep, which given each variation is up and down of one representation, should be 51-1
          float Nrep=LHEPdfWeight.size()/2.;
          weight[0]+=TMath::Sqrt(PDFuncertaintyUp  /(Nrep-1.)+pow(RenUncertaintyUp  ,2)+pow(FactUncertaintyUp  ,2));
          weight[1]-=TMath::Sqrt(PDFuncertaintyDown/(Nrep-1.)+pow(RenUncertaintyDown,2)+pow(FactUncertaintyDown,2));
          return weight;
        }
        
        RVec<float> CalculateMuonScaleFactor(const RVec<float> &Trigger, const RVec<float> &Trigger_stat, const RVec<float> &ID, const RVec<float> &ID_sys, const RVec<float> &ID_stat, const RVec<float> &ISO, const RVec<float> &ISO_sys, const RVec<float> &ISO_stat){
          RVec<float> weights = {1.,1.,1.};
          for(unsigned m = 0; m < Trigger.size(); ++m){
            weights[0]*=ID[m]*ISO[m];
            weights[1]*=(ID[m]+TMath::Sqrt(pow(ID_stat[m],2)+pow(ID_sys[m],2)))*(ISO[m]+TMath::Sqrt(pow(ISO_stat[m],2)+pow(ISO_sys[m],2)));
            weights[2]*=(ID[m]-TMath::Sqrt(pow(ID_stat[m],2)+pow(ID_sys[m],2)))*(ISO[m]-TMath::Sqrt(pow(ISO_stat[m],2)+pow(ISO_sys[m],2)));
          }
          return weights;
        }

        RVec<float> CalcMuonRecoIdIsoSFPerMuon(const RVec<float> &Trigger, const RVec<float> &Trigger_stat, const RVec<float> &ID, const RVec<float> &ID_sys, const RVec<float> &ID_stat, const RVec<float> &ISO, const RVec<float> &ISO_sys, const RVec<float> &ISO_stat){
          RVec<float> weights = {};
          for(unsigned m = 0; m < Trigger.size(); ++m){
            weights.push_back(
            ((ID[m]+TMath::Sqrt(pow(ID_stat[m],2)+pow(ID_sys[m],2)))*
            (ISO[m]+TMath::Sqrt(pow(ISO_stat[m],2)+pow(ISO_sys[m],2))))
            - 
            ID[m]*ISO[m]
            );
          }
          return weights;
        }
        
        RVec<float> CalculateMuonTriggerEff(const RVec<float> &Trigger, const RVec<float> &Trigger_stat, const RVec<float> &Trigger_sys){
          RVec<float> weights = {1.,1.,1.};
          for(unsigned m = 0; m < Trigger.size(); ++m){
            weights[0]*=1+(Trigger[m]-1);
            weights[1]*=weights[0]*(1-Trigger_stat[m])*(1-Trigger_sys[m]);
            weights[2]*=weights[0]*(1-Trigger_stat[m])*(1-Trigger_sys[m]);
          }
          weights[1] = 1- weights[1];
          weights[2] = 1 - weights[2];
          return weights;
        }
        
        RVec<float> SysPercPerObj(const RVec<float> nom, const RVec<float> up, const RVec<float> down, int verbose){
          RVec<float> weights = {};
          for(unsigned m = 0; m < nom.size(); ++m){
              if (verbose>0) std::cout << abs(up[m]-down[m])/(2*nom[m]) << " " << up[m] << " " << down[m] << " " << nom[m] << std::endl;
            weights.push_back( abs(up[m]-down[m])/(2*nom[m]));
          }
          return weights;
        }        
        
        RVec<float> CalculateElectronScaleFactor(const RVec<float> &centralSF, const RVec<float> &stat, const RVec<float> &syst){
          RVec<float> weights = {1.,1.,1.};
          for(unsigned e = 0; e < centralSF.size(); ++e){
            float weight = centralSF[e];
            weights[0]*=centralSF[e];
            weights[1]*=centralSF[e]+TMath::Sqrt(pow(stat[e],2)+pow(syst[e],2));
            weights[2]*=centralSF[e]-TMath::Sqrt(pow(stat[e],2)+pow(syst[e],2));
          }
          return weights;
        }
        
        int oppositeSign(const RVec<int> &charge){
          int oppositSignFloat = 0;
          if (charge.size()==2){
            if (charge[0]*charge[1]<0) oppositSignFloat=1;
          }
          return oppositSignFloat;
        }
        
        int oppositeSign(const RVec<int> &charge1,const RVec<int> &charge2){
          int oppositSignFloat = 0;
          if (charge1.size()==1 and charge2.size()==1){
            if (charge1[0]*charge2[0]<0) oppositSignFloat=1;
          }
          return oppositSignFloat;
        }
        
        float map_zero_to_one(float value){
          if (value==0){
            return (float) 1.0;
          } else {
            return value;
          }
        }
        
        
        RVec<float> find_replace(const RVec<float> &find_replace_vec, float find, float replace){  
            RVec<float> new_vec = {};
             for(unsigned i = 0; i < find_replace_vec.size(); ++i){
                 if (find_replace_vec[i] == find){
                      new_vec.push_back(replace);
                 } else {
                     new_vec.push_back(find_replace_vec[i]);
                 }
             }
             return new_vec;
        }
        
        float k_factor(int era, float mass){
          float a = 1.067;
          float b = -0.000112;
          float c = 3.176*pow(10,-8);
          float d = -4.068*pow(10,-12);
          float k_fac = 1;
          float X = mass-400;
          k_fac = a+b*X+c*pow(X,2)+d*pow(X,3);
          return k_fac;
        }
        
        ''')
