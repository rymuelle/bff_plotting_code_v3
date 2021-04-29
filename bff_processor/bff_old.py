##!/usr/bin/env python
#
#from __future__ import print_function
#from glob import glob
#from ROOT import vector, RDataFrame, RDF, TFile, TH1F, TH2F, gInterpreter, TMath
#import ROOT
#import sys
#
#gInterpreter.Declare('''
##include <math.h>
#
#TH2F *bTagEff = 0;
#
#using namespace ROOT::VecOps;
#float GetBTagWeight(const RVec<int> &isBJet, const RVec<int> &JetHadronFlav, const RVec<float> &JetPt, const RVec<float> &bTagSF) {
#  float weight = 1.;
#  for (unsigned j = 0; j < isBJet.size(); ++j) {
#    if (isBJet[j]) {
#      weight *= bTagSF[j];
#    }  else {
#    int HadronFlav = 2.;
#    if (JetHadronFlav[j] == 5)
#      HadronFlav=0.;
#    else if (JetHadronFlav[j] == 4)
#      HadronFlav=1.;
#      int effBin = bTagEff->FindBin(HadronFlav, JetPt[j]);
#      float eff  = bTagEff->GetBinContent(effBin);
#      weight *= (1. - eff * bTagSF[j]) / (1. - eff);
#    }
#  }
#  return weight;
#}
#
#TH2F *PUIDSF_true  = 0;
#TH2F *PUIDSF_false = 0;
#TH2F *PUIDUnc_true = 0;
#TH2F *PUIDUnc_false= 0;
#TH2F *PUIDEff_true = 0;
#TH2F *PUIDEff_false= 0;
#
#RVec<float> GetPUIDweight(const RVec<float> &JetPt, const RVec<float> &JetEta, const RVec<int> &JetGenJetIdx, const RVec<int> &PUID) {
#  RVec<float> weights = {1.,1.,1.};
#  for(unsigned j = 0; j < JetGenJetIdx.size(); ++j) {
#    const int SFbin = PUIDSF_true->FindBin(JetPt[j], JetEta[j]);
#    if(PUID[j] & 1){
#      if(JetGenJetIdx[j]>=0){
#        const float SF = PUIDSF_true->GetBinContent(SFbin);
#        const float unc= PUIDUnc_true->GetBinContent(SFbin);
#        weights[0] *= SF;
#        weights[1] *= SF+unc;
#        weights[2] *= SF-unc;
#      }
#      else{
#        const float SF = PUIDSF_false->GetBinContent(SFbin);
#        const float unc= PUIDUnc_false->GetBinContent(SFbin);
#        weights[0] *= SF;
#        weights[1] *= SF+unc;
#        weights[2] *= SF-unc;
#      }
#    }
#    else{
#      if(JetGenJetIdx[j]>=0){
#        const float SF = PUIDSF_true->GetBinContent(SFbin);
#        const float unc= PUIDUnc_true->GetBinContent(SFbin);
#        const float eff= PUIDEff_true->GetBinContent(SFbin);
#        weights[0] *= (1-eff*SF)/(1-eff);
#        weights[1] *= (1-eff*(SF+unc))/(1-eff);
#        weights[2] *= (1-eff*(SF-unc))/(1-eff);
#      }
#      else{
#        const float SF = PUIDSF_false->GetBinContent(SFbin);
#        const float unc= PUIDUnc_false->GetBinContent(SFbin);
#        const float eff= PUIDEff_false->GetBinContent(SFbin);
#        weights[0] *= (1-eff*SF)/(1-eff);
#        weights[1] *= (1-eff*(SF+unc))/(1-eff);
#        weights[2] *= (1-eff*(SF-unc))/(1-eff);
#      }
#    }
#  }
#  return weights;
#}
#
#RVec<float> GetPDFandScaleUncertainty(const RVec<float> &LHEPdfWeight, const RVec<float> &LHEScaleWeight){
#  RVec<float> weight={1.,1.};
#  if(LHEPdfWeight.size()==0 || LHEScaleWeight.size()==0) return weight;
#  //according to Eq. 3 https://arxiv.org/pdf/1101.0536.pdf, 103 members expected according to https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopSystematics#PDF, 102 variation members stored, relative to nominal (1)
#  float PDFuncertaintyUp=0.;
#  float PDFuncertaintyDown=0.;
#  for(unsigned i=0; i<LHEPdfWeight.size()/2; ++i){
#    PDFuncertaintyUp  += pow(TMath::Max(TMath::Max(LHEPdfWeight[2*i]-1.,LHEPdfWeight[2*i+1]-1.),0.),2);
#    PDFuncertaintyDown+= pow(TMath::Max(TMath::Max(1.-LHEPdfWeight[2*i],1.-LHEPdfWeight[2*i+1]),0.),2);
#  }
#  //Scale weight usage based on conjecture that this works by choosing the maximal variation in each direction for renormalization and factorization, independently. Contents as of https://cms-nanoaod-integration.web.cern.ch/integration/master/mc94X_doc.html#LHE
#  float RenUncertaintyUp=TMath::Max(LHEScaleWeight[1]-1., TMath::Max(LHEScaleWeight[4]-1., TMath::Max(LHEScaleWeight[7]-1., 0.)));
#  float RenUncertaintyDown=TMath::Max(1.-LHEScaleWeight[1], TMath::Max(1.-LHEScaleWeight[4], TMath::Max(1.-LHEScaleWeight[7], 0.)));
#  float FactUncertaintyUp=TMath::Max(LHEScaleWeight[3]-1., TMath::Max(LHEScaleWeight[4]-1., TMath::Max(LHEScaleWeight[5]-1., 0.)));
#  float FactUncertaintyDown=TMath::Max(1.-LHEScaleWeight[3], TMath::Max(1.-LHEScaleWeight[4], TMath::Max(1.-LHEScaleWeight[5], 0.)));
#  //according to Eq. 8 or 10 in https://arxiv.org/pdf/1101.0536.pdf, modified by Eq. 6 to a factor of 1/N_rep, which given each variation is up and down of one representation, should be 51-1
#  float Nrep=LHEPdfWeight.size()/2.;
#  weight[0]+=TMath::Sqrt(PDFuncertaintyUp  /(Nrep-1.)+pow(RenUncertaintyUp  ,2)+pow(FactUncertaintyUp  ,2));
#  weight[1]-=TMath::Sqrt(PDFuncertaintyDown/(Nrep-1.)+pow(RenUncertaintyDown,2)+pow(FactUncertaintyDown,2));
#  return weight;
#}
#
#RVec<float> CalculateMuonScaleFactor(const RVec<float> &Trigger, const RVec<float> &Trigger_stat, const RVec<float> &ID, const RVec<float> &ID_sys, const RVec<float> &ID_stat, const RVec<float> &ISO, const RVec<float> &ISO_sys, const RVec<float> &ISO_stat){
#  RVec<float> weights = {1.,1.,1.};
#  for(unsigned m = 0; m < Trigger.size(); ++m){
#    weights[0]*=Trigger[m]*ID[m]*ISO[m];
#    weights[1]*=(Trigger[m]+Trigger_stat[m])*(ID[m]+TMath::Sqrt(pow(ID_stat[m],2)+pow(ID_sys[m],2)))*(ISO[m]+TMath::Sqrt(pow(ISO_stat[m],2)+pow(ISO_sys[m],2)));
#    weights[2]*=(Trigger[m]-Trigger_stat[m])*(ID[m]-TMath::Sqrt(pow(ID_stat[m],2)+pow(ID_sys[m],2)))*(ISO[m]-TMath::Sqrt(pow(ISO_stat[m],2)+pow(ISO_sys[m],2)));
#  }
#  return weights;
#}
#
#RVec<float> CalculateElectronScaleFactor(const RVec<float> &centralSF, const RVec<float> &stat, const RVec<float> &syst){
#  RVec<float> weights = {1.,1.,1.};
#  for(unsigned e = 0; e < centralSF.size(); ++e){
#    float weight = centralSF[e];
#    weights[0]*=centralSF[e];
#    weights[1]*=centralSF[e]+TMath::Sqrt(pow(stat[e],2)+pow(syst[e],2));
#    weights[2]*=centralSF[e]-TMath::Sqrt(pow(stat[e],2)+pow(syst[e],2));
#  }
#  return weights;
#}
#
#int oppositeSign(const RVec<int> &charge){
#  int oppositSignFloat = 0;
#  if (charge.size()==2){
#    if (charge[0]*charge[1]<0) oppositSignFloat=1;
#  }
#  return oppositSignFloat;
#}
#
#int oppositeSign(const RVec<int> &charge1,const RVec<int> &charge2){
#  int oppositSignFloat = 0;
#  if (charge1.size()==1 and charge2.size()==1){
#    if (charge1[0]*charge2[0]<0) oppositSignFloat=1;
#  }
#  return oppositSignFloat;
#}
#
#float map_zero_to_one(float value){
#  if (value==0){
#    return 1;
#  } else {
#    return value;
#  }
#}
#
#float k_factor(int era, float mass){
#  float a = 1.067;
#  float b = -0.000112;
#  float c = 3.176*pow(10,-8);
#  float d = -4.068*pow(10,-12);
#  float k_fac = 1;
#  float X = mass-400;
#  k_fac = a+b*X+c*pow(X,2)+d*pow(X,3);
#  return k_fac;
#}
#
#''')
#
#
#def toVector(vtype, list):
#    v = vector(vtype)()
#    for o in list:
#        v.push_back(o)
#    return v
#
#
#def process_samples(samples='samplesCR_2016_Apr2020ReReco_lo.py', outname='bff.root', is_inclusive=False, useSpark=False, npartitions=32, is_quick=False):
#    if useSpark:
#        import PyRDF
#        PyRDF.use('spark', {'npartitions': npartitions})
#        RDFrame = PyRDF.RDataFrame
#    else:
#        ROOT.ROOT.EnableImplicitMT()
#        RDFrame = RDataFrame
#
#    #load python object with samples
#    print("opening: {}".format(samples))
#    exec(open(samples).read(),locals(),globals())    
#    assert "lumi" in sample_dict, "lumi dict must have lumi key."
#    lumi = sample_dict['lumi']
#
#    lumi = float(lumi)
#
#    out = TFile(outname, 'recreate')
#    hlumi = TH1F("lumi", "lumi", 1, 0, 1)
#    hlumi.SetDirectory(out)
#    hlumi.SetBinContent(1, lumi)
#    hlumi.Write()
#
#    sampleswitch = 0
#    if   '2016' in samples:
#        bDiscValue = 0.6321
#        sampleswitch = 2016
#    elif '2017' in samples:
#        bDiscValue = 0.4941
#        sampleswitch = 2017
#    elif '2018' in samples:
#        bDiscValue = 0.4184
#        sampleswitch = 2018
#    else:
#        sys.exit("no era specified in sample name")
#
#    samples = {}
#
#    for sample in sample_dict['samples']:
#
#        name = sample["name"]
#        xsec = sample["xsec"]
#        nevts = sample["nevts"]
#        ismc = sample["ismc"]
#        fileglob = sample["fileglob"]
#        bTagEff= sample["bTagEff"]
#
#        samples[name] = int(ismc)
#        outdir = out.mkdir(name)
#        print(name, xsec, nevts, ismc)
#        files = toVector('string', glob(fileglob))
#        if int(ismc):
#             bTagFile = TFile(bTagEff, 'read')
#             ROOT.bTagEff = bTagFile.Get("bTagEff")
#             PUIDSFfile = TFile("../data/PUID_SF/PUID_80XTraining_EffSFandUncties.root", 'read')
#             if sampleswitch == 2016:
#                 ROOT.PUIDEff_true  = PUIDSFfile.Get("h2_eff_mc2016_T")
#                 ROOT.PUIDEff_false = PUIDSFfile.Get("h2_mistag_mc2016_T")
#                 ROOT.PUIDSF_true   = PUIDSFfile.Get("h2_eff_sf2016_T")
#                 ROOT.PUIDSF_false  = PUIDSFfile.Get("h2_mistag_sf2016_T")
#                 ROOT.PUIDUnc_true  = PUIDSFfile.Get("h2_eff_sf2016_T_Systuncty")
#                 ROOT.PUIDUnc_false = PUIDSFfile.Get("h2_mistag_sf2016_T_Systuncty")
#             elif sampleswitch == 2017:
#                 ROOT.PUIDEff_true  = PUIDSFfile.Get("h2_eff_mc2017_T")
#                 ROOT.PUIDEff_false = PUIDSFfile.Get("h2_mistag_mc2017_T")
#                 ROOT.PUIDSF_true   = PUIDSFfile.Get("h2_eff_sf2017_T")
#                 ROOT.PUIDSF_false  = PUIDSFfile.Get("h2_mistag_sf2017_T")
#                 ROOT.PUIDUnc_true  = PUIDSFfile.Get("h2_eff_sf2017_T_Systuncty")
#                 ROOT.PUIDUnc_false = PUIDSFfile.Get("h2_mistag_sf2017_T_Systuncty")
#             elif sampleswitch == 2018:
#                 ROOT.PUIDEff_true  = PUIDSFfile.Get("h2_eff_mc2018_T")
#                 ROOT.PUIDEff_false = PUIDSFfile.Get("h2_mistag_mc2018_T")
#                 ROOT.PUIDSF_true   = PUIDSFfile.Get("h2_eff_sf2018_T")
#                 ROOT.PUIDSF_false  = PUIDSFfile.Get("h2_mistag_sf2018_T")
#                 ROOT.PUIDUnc_true  = PUIDSFfile.Get("h2_eff_sf2018_T_Systuncty")
#                 ROOT.PUIDUnc_false = PUIDSFfile.Get("h2_mistag_sf2018_T_Systuncty")
        # Uncomment this line to run on swan with spark
        # files = [f.replace("/eos/cms","root://eoscms.cern.ch//") for f in files]
        df = RDFrame('Events', files)
        df = df.Define("GoodJet", "Jet_pt>20 && abs(Jet_eta)<2.4 && (Jet_jetId & (1 << 1)) && ((Jet_puId & 1) || Jet_pt>50.) && !(Jet_btagDeepB<={} && Jet_pt<=30.)".format(bDiscValue))\
               .Define("GoodJetPt", "Jet_pt[GoodJet]")\
               .Define("GoodJetEta", "Jet_eta[GoodJet]")\
               .Define("nJets", "GoodJetPt.size()")\
               .Define("BJet", "Jet_pt>20 && abs(Jet_eta)<2.4 && (Jet_jetId & (1 << 1)) && (Jet_puId & 1) && Jet_btagDeepB>{}".format(bDiscValue))\
               .Define("GoodBJet", "BJet[GoodJet]")
        df = df.Define("leading_b_jet_pt", "(GoodBJet.size()>0) ? Jet_pt[GoodBJet[0]] : 0")
        df = df.Define("leading_jet_pt", "(GoodJet.size()>0) ? Jet_pt[GoodJet[0]] : 0")

        df = df.Define("NoPUID_GoodJet", "Jet_pt>20 && abs(Jet_eta)<2.4 && (Jet_jetId & (1 << 1)) && !(Jet_btagDeepB<={} && Jet_pt<=30.) && Jet_pt<50.".format(bDiscValue))\
               .Define("NoPUID_GoodJetPt", "Jet_pt[NoPUID_GoodJet]")\
               .Define("NoPUID_GoodJetEta", "Jet_eta[NoPUID_GoodJet]")\
               .Define("NoPUID_PUID", "Jet_puId[NoPUID_GoodJet]")
        if int(ismc):
          df = df.Define("NoPUID_GoodJetGenJetIdx", "Jet_genJetIdx[NoPUID_GoodJet]")
        df = df.Define("GoodMuon", "Muon_pt_corrected > 53 && abs(Muon_eta) < 2.4 && Muon_tightId > 0 && Muon_pfRelIso04_all < 0.25")\
               .Define("GoodMuonPt", "Muon_pt[GoodMuon]")\
               .Define("GoodMuonEta", "Muon_eta[GoodMuon]")\
               .Define("GoodMuonPhi", "Muon_phi[GoodMuon]")\
               .Define("GoodMuonCharge", "Muon_charge[GoodMuon]")\
               .Define("nMuons", "GoodMuonPt.size()")
        df = df.Define("GoodMuonLowPt", "Muon_pt_corrected > 24 && abs(Muon_eta) < 2.4 && Muon_tightId > 0 && Muon_pfRelIso04_all < 0.25")\
               .Define("GoodMuonPtLow", "Muon_pt[GoodMuonLowPt]")\
               .Define("GoodMuonEtaLow", "Muon_eta[GoodMuonLowPt]")\
               .Define("nMuonsLowPt", "GoodMuonPtLow.size()")
        df = df.Define("GoodElectron", "Electron_pt > 53 && abs(Electron_eta) < 2.4 && Electron_cutBased_HEEP > 0")\
               .Define("GoodElePt", "Electron_pt[GoodElectron]")\
               .Define("GoodEleEta", "Electron_eta[GoodElectron]")\
               .Define("GoodElePhi", "Electron_phi[GoodElectron]")\
               .Define("GoodEleCharge", "Electron_charge[GoodElectron]")\
               .Define("nEle", "GoodElePt.size()")
        df = df.Define("GoodElectronLowPt", "Electron_pt > 24 && abs(Electron_eta) < 2.4 && Electron_cutBased_HEEP > 0")\
              .Define("GoodElePtLow", "Electron_pt[GoodElectronLowPt]")\
              .Define("GoodEleEtaLow", "Electron_eta[GoodElectronLowPt]")\
               .Define("nEleLowPt", "GoodElePtLow.size()")
        if int(ismc):
          df = df.Define("GoodMuon_effSF_trigger", "Muon_effSF_trigger[GoodMuon]")\
               .Define("GoodMuon_effSF_stat_trigger", "Muon_effSF_stat_trigger[GoodMuon]")\
               .Define("GoodMuon_effSF_ID", "Muon_effSF_ID[GoodMuon]")\
               .Define("GoodMuon_effSF_sys_ID", "Muon_effSF_sys_ID[GoodMuon]")\
               .Define("GoodMuon_effSF_stat_ID", "Muon_effSF_stat_ID[GoodMuon]")\
               .Define("GoodMuon_effSF_ISO", "Muon_effSF_ISO[GoodMuon]")\
               .Define("GoodMuon_effSF_sys_ISO", "Muon_effSF_sys_ISO[GoodMuon]")\
               .Define("GoodMuon_effSF_stat_ISO", "Muon_effSF_stat_ISO[GoodMuon]")\
               .Define("GoodElectron_effSF", "Electron_effSF[GoodElectron]")\
               .Define("GoodElectron_effSF_stat", "Electron_effSF_stat[GoodElectron]")\
               .Define("GoodElectron_effSF_sys", "Electron_effSF_sys[GoodElectron]")

        if int(ismc):
            if sampleswitch == 2016:
              df = df.Define("TriggerRegion1", "HLT_Mu50>0 or HLT_DoubleEle33_CaloIdL_MW>0")
              df = df.Define("TriggerRegion2", "HLT_Mu50>0 or HLT_TkMu50>0 or HLT_DoubleEle33_CaloIdL_MW>0")
              df = df.Define("TriggerRegion3", "HLT_Mu50>0 or HLT_TkMu50>0 or HLT_DoubleEle33_CaloIdL_MW>0 or HLT_DoubleEle33_CaloIdL_GsfTrkIdVL>0")
              df = df.Define("TriggerWeight", "TriggerRegion1*0.073+TriggerRegion2*(1.0-0.073-0.3087428721)+TriggerRegion3*0.3087428721")
            elif sampleswitch == 2017:
              df = df.Define("TriggerRegion1", "HLT_Mu50>0 or HLT_OldMu100>0 or HLT_TkMu100>0 or HLT_DoubleEle33_CaloIdL_MW>0")
              df = df.Define("TriggerRegion2", "HLT_Mu50>0 or HLT_OldMu100>0 or HLT_TkMu100>0 or HLT_DoubleEle33_CaloIdL_MW>0 or HLT_DoubleEle25_CaloIdL_MW>0")
              df = df.Define("TriggerWeight", "TriggerRegion1*(1-0.6528268793)+TriggerRegion2*0.6528268793")
            elif sampleswitch == 2018:
              df = df.Define("TriggerRegion1", "HLT_Mu50>0 or HLT_OldMu100>0 or HLT_TkMu100>0 or HLT_DoubleEle25_CaloIdL_MW>0")
              df = df.Define("TriggerWeight", "TriggerRegion1")
            sample_weight = float(xsec)*lumi/float(nevts)
            print(name)
            if name not in ['mc_zz', 'mc_wz']:
              df = df.Define("PDF_ISRFSR_uncertainty","GetPDFandScaleUncertainty(LHEPdfWeight, LHEScaleWeight)")
              df = df.Define("PDF_ISRFSR_uncertaintyUp","PDF_ISRFSR_uncertainty[0]")
              df = df.Define("PDF_ISRFSR_uncertaintyDown","PDF_ISRFSR_uncertainty[1]")
            else:
              df = df.Define("PDF_ISRFSR_uncertainty","1")
              df = df.Define("PDF_ISRFSR_uncertaintyUp","1")
              df = df.Define("PDF_ISRFSR_uncertaintyDown","1")

            df = df.Define("MuonSFweights","CalculateMuonScaleFactor(GoodMuon_effSF_trigger, GoodMuon_effSF_stat_trigger, GoodMuon_effSF_ID, GoodMuon_effSF_sys_ID, GoodMuon_effSF_stat_ID, GoodMuon_effSF_ISO, GoodMuon_effSF_sys_ISO, GoodMuon_effSF_stat_ISO)")
            df = df.Define("MuonSFweight","MuonSFweights[0]")
            df = df.Define("MuonSFweightUp","MuonSFweights[1]")
            df = df.Define("MuonSFweightDown","MuonSFweights[2]")
            df = df.Define("ElectronSFweights","CalculateElectronScaleFactor(GoodElectron_effSF, GoodElectron_effSF_stat, GoodElectron_effSF_sys)")
            df = df.Define("ElectronSFweight","map_zero_to_one(ElectronSFweights[0])")
            df = df.Define("ElectronSFweightUp","map_zero_to_one(ElectronSFweights[1])")
            df = df.Define("ElectronSFweightDown","map_zero_to_one(ElectronSFweights[2])")
            if "ZTo" in name: 
              print("applying k_factor")
              df = df.Define("k_factor","k_factor(2016, DiLepMass)")
            else: df = df.Define("k_factor","1")
            if not is_inclusive:
                df = df.Define("PUIDWeights", "GetPUIDweight(NoPUID_GoodJetPt, NoPUID_GoodJetEta, NoPUID_GoodJetGenJetIdx, NoPUID_PUID)")
                df = df.Define("PUIDWeight", "map_zero_to_one(PUIDWeights[0])")
                df = df.Define("PUIDWeightUp", "map_zero_to_one(PUIDWeights[1])")
                df = df.Define("PUIDWeightDown", "map_zero_to_one(PUIDWeights[2])")
                df = df.Define("GoodJetBTagSF", "(Jet_btagSF[GoodJet])")
                df = df.Define("GoodJetBTagSFUp", "(Jet_btagSF_up[GoodJet])")
                df = df.Define("GoodJetBTagSFDown", "(Jet_btagSF_down[GoodJet])")
                df = df.Define("GoodJetHadronFlav","Jet_hadronFlavour[GoodJet]")
                df = df.Define("BTagWeight", "map_zero_to_one(GetBTagWeight(GoodBJet, GoodJetHadronFlav, GoodJetPt, GoodJetBTagSF))")
                df = df.Define("BTagWeightUp", "map_zero_to_one(GetBTagWeight(GoodBJet, GoodJetHadronFlav, GoodJetPt, GoodJetBTagSFUp))")
                df = df.Define("BTagWeightDown", "map_zero_to_one(GetBTagWeight(GoodBJet, GoodJetHadronFlav, GoodJetPt, GoodJetBTagSFDown))")

                df = df.Define("Weight", "{}*copysign(1.,genWeight)*k_factor*puWeight*BTagWeight*PUIDWeight*MuonSFweight*ElectronSFweight*TriggerWeight".format(sample_weight))
                # Systematic weights
                df = df.Define("Weight_PuUp", "Weight/puWeight*puWeightUp")
                df = df.Define("Weight_PuDown", "Weight/puWeight*puWeightDown")
                df = df.Define("Weight_BTagUp", "Weight/BTagWeight*BTagWeightUp")
                df = df.Define("Weight_BTagDown", "Weight/BTagWeight*BTagWeightDown")
                df = df.Define("Weight_PUIDUp", "Weight/PUIDWeight*PUIDWeightUp")
                df = df.Define("Weight_PUIDDown", "Weight/PUIDWeight*PUIDWeightDown")
                df = df.Define("Weight_PDF_ISRFSR_Up", "Weight*PDF_ISRFSR_uncertaintyUp")
                df = df.Define("Weight_PDF_ISRFSR_Down", "Weight*PDF_ISRFSR_uncertaintyDown")
                df = df.Define("Weight_MuonSFUp", "Weight/MuonSFweight*MuonSFweightUp")
                df = df.Define("Weight_MuonSFDown", "Weight/MuonSFweight*MuonSFweightDown")
                df = df.Define("Weight_ElectronSFUp", "Weight/ElectronSFweight*ElectronSFweightUp")
                df = df.Define("Weight_ElectronSFDown", "Weight/ElectronSFweight*ElectronSFweightDown")
            else: # Inclusive sample doesn't cut on btags, so no btag weight needed; also no PUID needed
                df = df.Define("Weight", "{}*copysign(1.,genWeight)*k_factor*puWeight*MuonSFweight*ElectronSFweight*TriggerWeight".format(sample_weight))
                df = df.Define("Weight_PuUp", "Weight/puWeight*puWeightUp")
                df = df.Define("Weight_PuDown", "Weight/puWeight*puWeightDown")
                df = df.Define("Weight_PDF_ISRFSR_Up", "Weight*PDF_ISRFSR_uncertaintyUp")
                df = df.Define("Weight_PDF_ISRFSR_Down", "Weight*PDF_ISRFSR_uncertaintyDown")
                df = df.Define("Weight_MuonSFUp", "Weight/MuonSFweight*MuonSFweightUp")
                df = df.Define("Weight_MuonSFDown", "Weight/MuonSFweight*MuonSFweightDown")
                df = df.Define("Weight_ElectronSFUp", "Weight/ElectronSFweight*ElectronSFweightUp")
                df = df.Define("Weight_ElectronSFDown", "Weight/ElectronSFweight*ElectronSFweightDown")
        else:
            df = df.Define("Weight", "1")
            df = df.Define("Weight_PuDown", "1")
            df = df.Define("Weight_PuUp", "1")
            df = df.Define("Weight_BTagUp", "1")
            df = df.Define("Weight_BTagDown", "1")
            df = df.Define("Weight_PUIDUp", "1")
            df = df.Define("Weight_PUIDDown", "1")
            df = df.Define("Weight_PDF_ISRFSR_Up", "1")
            df = df.Define("Weight_PDF_ISRFSR_Down", "1")
            df = df.Define("Weight_MuonSFUp", "1")
            df = df.Define("Weight_MuonSFDown", "1")
            df = df.Define("Weight_ElectronSFUp", "1")
            df = df.Define("Weight_ElectronSFDown", "1")
            df = df.Define("k_factor","1")

        df = df.Define("oppositeSignMuon", "oppositeSign(GoodMuonCharge)")
        df = df.Define("oppositeSignElectron", "oppositeSign(GoodEleCharge)")
        df = df.Define("oppositeSignEmu", "oppositeSign(GoodMuonCharge,GoodEleCharge)")

        df = df.Define("IncMumuLocal", "nEle==0 && nMuons==2 && oppositeSignMuon==1")
        df = df.Define("IncMumuLocalLow", "IncMumuLocal && nMuonsLowPt==2 && nEleLowPt==0")
        df = df.Define("IncEeLocal", "nEle==2 && nMuons==0 && oppositeSignElectron==1")
        df = df.Define("IncEeLocalLow", "IncEeLocal && nEleLowPt==2 && nMuonsLowPt==0")
        df = df.Define("IncEmuLocal", "nEle==1 && nMuons==1 && oppositeSignEmu==1")
        df = df.Define("IncEmuLocalLow", "IncEmuLocal && nEleLowPt==1 && nMuonsLowPt==1")

        regions = []
        if is_inclusive:
            mu_zPeak = df.Filter("IncMumuLocalLow && DiLepMass>60 && DiLepMass<=120", "IncMuMu_zPeak")
            el_zPeak = df.Filter("IncEeLocalLow && DiLepMass>60 && DiLepMass<=120", "IncEe_zPeak")
            emu_zPeak = df.Filter("IncEmuLocalLow && DiLepMass>60 && DiLepMass<=120", "IncEmu_zPeak")

            mu = df.Filter("IncMumuLocalLow && DiLepMass>120", "IncMuMu")
            el = df.Filter("IncEeLocalLow && DiLepMass>120", "IncEe")
            emu = df.Filter("IncEmuLocalLow && DiLepMass>120", "IncEmu")

            #mu = df.Filter("IncMumuLocalLow && DiLepMass>50", "IncMuMu")
            #el = df.Filter("IncEeLocalLow && DiLepMass>50", "IncEe")
            #emu = df.Filter("IncEmuLocalLow && DiLepMass>50", "IncEmu")

            regions = zip([mu, el, emu, mu_zPeak, el_zPeak, emu_zPeak], ["IncMuMu", "IncEe", "IncEmu", "IncMuMu_zPeak", "IncEe_zPeak", "IncEmu_zPeak"])
        else:
            print(int(ismc))
            if int(ismc):
              JERC_var = ['nominal','jerUp','jerDown','jesUp','jesDown']
            else:
              JERC_var = ['nominal']
            rs = ["CR10", "CR11", "CR12", "CR13", "CR14", "CR20", "CR21", "CR22", "CR23", "CR24", "SR1", "SR2"]
            rs = [(r,var) for r in rs for var in JERC_var]
            for reg,var in rs:
                r = '{}_{}'.format(reg,var)
                HTLT_string = 'HTLT_{}'.format(var)
                RelMET_string = 'RelMET_{}'.format(var)
                SBM_string = 'SBM_{}'.format(var)

                HTLT,RelMET,SBM = -120,0.22,0
                if r[2] == "2":
                    HTLT,RelMET,SBM = -60,0.22,150
                
                region_string = "{}_no_zPeak_no_zkin_cut".format(r)
                format_dict = {"region_string":region_string, "r": r,"HTLT": HTLT,"RelMET": RelMET,"SBM": SBM,"HTLT_string": HTLT_string,"RelMET_string": RelMET_string,"SBM_string": SBM_string}
                regions.append((df.Filter("{r} && DiLepMass>120".format(**format_dict), region_string), format_dict))

                region_string = "{}_no_zPeak".format(r)
                format_dict = {"region_string":region_string, "r": r,"HTLT": HTLT,"RelMET": RelMET,"SBM": SBM,"HTLT_string": HTLT_string,"RelMET_string": RelMET_string,"SBM_string": SBM_string}
                regions.append((df.Filter("{r} && DiLepMass>120 && {HTLT_string}<{HTLT} && {RelMET_string}<{RelMET} && {SBM_string}>{SBM}".format(**format_dict), region_string), format_dict))
                if not is_quick:
                  region_string = "{}".format(r)
                  format_dict = {"region_string":region_string, "r": r,"HTLT": HTLT,"RelMET": RelMET,"SBM": SBM,"HTLT_string": HTLT_string,"RelMET_string": RelMET_string,"SBM_string": SBM_string}
                  regions.append((df.Filter("{r} && DiLepMass>54 && {HTLT_string}<{HTLT} && {RelMET_string}<{RelMET} && {SBM_string}>{SBM}".format(**format_dict), region_string), format_dict))
                  region_string = "{}_zPeak".format(r)
                  format_dict = {"region_string":region_string, "r": r,"HTLT": HTLT,"RelMET": RelMET,"SBM": SBM,"HTLT_string": HTLT_string,"RelMET_string": RelMET_string,"SBM_string": SBM_string}
                  regions.append((df.Filter("{r} && DiLepMass<=120 && DiLepMass>60 && {HTLT_string}<{HTLT} && {RelMET_string}<{RelMET} && {SBM_string}>{SBM}".format(**format_dict), region_string), format_dict))
        hists = []
        hists.append(df.Histo1D(RDF.TH1DModel("total"+"_Weight", ";Event Weight;Events", 100, -10, 10), "Weight"))

        for d, format_dict in regions:
            region = format_dict['region_string']
            nBins = 100
            min_mass, max_mass = 0, 600
            if "zPeak" in region: min_mass, max_mass, nBins = 60, 120, 100
            if "no_zPeak" in region: min_mass, max_mass, nBins = 120, 800, (800-120)/5

            HTLT_string = format_dict['HTLT_string']
            RelMET_string = format_dict['RelMET_string']
            SBM_string = format_dict['SBM_string']
            hists.append(d.Histo1D(RDF.TH1DModel(region+"_HTLT", ";H_{T}-L_{T};Events / 5 GeV", 100, -500, 500),HTLT_string, "Weight"))
            hists.append(d.Histo1D(RDF.TH1DModel(region+"_RelMET", "E_{T}^{miss}/M_{\ell^{+}\ell^{-}};Events / 5 GeV", 100, 0, 1),RelMET_string, "Weight"))
            hists.append(d.Histo1D(RDF.TH1DModel(region+"_SBM", ";max(\{ M_{jet_{i}\ell_{j}}  , M_{jet_{k}\ell_{l}}\}_{min\Delta M} );Events / 5 GeV", 100, 0, 300),SBM_string, "Weight"))

            hists.append(d.Histo1D(RDF.TH1DModel(region+"_DiLepMass", ";Dilepton Mass [GeV];Events / 5 GeV", nBins, min_mass, max_mass), "DiLepMass", "Weight"))
            if int(ismc):
              hists.append(d.Histo1D(RDF.TH1DModel(region+"_W_PU_central", ";PU weight; Events / 0.01 weight", 200, 0., 2.), "puWeight"))         
              hists.append(d.Histo1D(RDF.TH1DModel(region+"_W_PU_up", ";PU weight up-variation; Events / 0.01 weight", 200, 0., 2.), "puWeightUp"))
              hists.append(d.Histo1D(RDF.TH1DModel(region+"_W_PU_down", ";PU weight down-variation; Events / 0.01 weight", 200, 0., 2.), "puWeightDown"))
              if not is_inclusive:
                hists.append(d.Histo1D(RDF.TH1DModel(region+"_W_BTagSF_central", ";BTag weight; Events / 0.01 weight", 200, 0., 2.), "BTagWeight"))
                hists.append(d.Histo1D(RDF.TH1DModel(region+"_W_BTagSF_up", ";BTag weight up-variation; Events / 0.01 weight", 200, 0., 2.), "BTagWeightUp"))
                hists.append(d.Histo1D(RDF.TH1DModel(region+"_W_BTagSF_down", ";BTag weight down-variation; Events / 0.01 weight", 200, 0., 2.), "BTagWeightDown"))
                hists.append(d.Histo1D(RDF.TH1DModel(region+"_W_PUIDSF_central", ";PUID weight; Events / 0.01 weight", 200, 0., 2.), "PUIDWeight"))
                hists.append(d.Histo1D(RDF.TH1DModel(region+"_W_PUIDSF_up", ";PUID weight up-variation; Events / 0.01 weight", 200, 0., 2.), "PUIDWeightUp"))
                hists.append(d.Histo1D(RDF.TH1DModel(region+"_W_PUIDSF_down", ";PUID weight down-variation; Events / 0.01 weight", 200, 0., 2.), "PUIDWeightDown"))
              hists.append(d.Histo1D(RDF.TH1DModel(region+"_W_PDF_ISRFSR_uncertainty", ";PDF+ISR/FSR uncertainty; Events / 0.01 weight", 250, -2., 3.), "PDF_ISRFSR_uncertainty"))
              hists.append(d.Histo1D(RDF.TH1DModel(region+"_W_PDF_ISRFSR_uncertainty_up", ";PDF+ISR/FSR uncertainty up-variation; Events / 0.01 weight", 250, -2., 3.), "PDF_ISRFSR_uncertaintyUp"))
              hists.append(d.Histo1D(RDF.TH1DModel(region+"_W_PDF_ISRFSR_uncertainty_down", ";PDF+ISR/FSR uncertainty down-variation; Events / 0.01 weight", 250, -2., 3.), "PDF_ISRFSR_uncertaintyDown"))
              hists.append(d.Histo1D(RDF.TH1DModel(region+"_W_MuonSF_central", ";MuonSF weight; Events / 0.01 weight", 200, 0., 2.), "MuonSFweight"))
              hists.append(d.Histo1D(RDF.TH1DModel(region+"_W_MuonSF_up", ";MuonSF weight up-variation; Events / 0.01 weight", 200, 0., 2.), "MuonSFweightUp"))
              hists.append(d.Histo1D(RDF.TH1DModel(region+"_W_MuonSF_down", ";MuonSF weight down-variation; Events / 0.01 weight", 200, 0., 2.), "MuonSFweightDown"))
              hists.append(d.Histo1D(RDF.TH1DModel(region+"_W_ElectronSF_central", ";ElectronSF weight; Events / 0.01 weight", 200, 0., 2.), "ElectronSFweight"))
              hists.append(d.Histo1D(RDF.TH1DModel(region+"_W_ElectronSF_up", ";ElectronSF weight up-variation; Events / 0.01 weight", 200, 0., 2.), "ElectronSFweightUp"))
              hists.append(d.Histo1D(RDF.TH1DModel(region+"_W_ElectronSF_down", ";ElectronSF weight down-variation; Events / 0.01 weight", 200, 0., 2.), "ElectronSFweightDown"))

            hists.append(d.Histo1D(RDF.TH1DModel(region+"_kFac", ";kFac;Events / 5 GeV", 100, 0,2), "k_factor", "Weight"))

            hists.append(d.Histo1D(RDF.TH1DModel(region+"_DiLepMass_Weight_PuDown", ";Dilepton Mass [GeV];Events / 5 GeV", nBins, min_mass, max_mass), "DiLepMass", "Weight_PuDown"))
            hists.append(d.Histo1D(RDF.TH1DModel(region+"_DiLepMass_Weight_PDF_ISRFSR_Down", ";Dilepton Mass [GeV];Events / 5 GeV", nBins, min_mass, max_mass), "DiLepMass", "Weight_PDF_ISRFSR_Down"))
            hists.append(d.Histo1D(RDF.TH1DModel(region+"_DiLepMass_Weight_MuonSFDown", ";Dilepton Mass [GeV];Events / 5 GeV", nBins, min_mass, max_mass), "DiLepMass", "Weight_MuonSFDown"))
            hists.append(d.Histo1D(RDF.TH1DModel(region+"_DiLepMass_Weight_ElectronSFDown", ";Dilepton Mass [GeV];Events / 5 GeV", nBins, min_mass, max_mass), "DiLepMass", "Weight_ElectronSFDown"))

            hists.append(d.Histo1D(RDF.TH1DModel(region+"_DiLepMass_Weight_PuUp", ";Dilepton Mass [GeV];Events / 5 GeV", nBins, min_mass, max_mass), "DiLepMass", "Weight_PuUp"))
            hists.append(d.Histo1D(RDF.TH1DModel(region+"_DiLepMass_Weight_PDF_ISRFSR_Up", ";Dilepton Mass [GeV];Events / 5 GeV", nBins, min_mass, max_mass), "DiLepMass", "Weight_PDF_ISRFSR_Up"))
            hists.append(d.Histo1D(RDF.TH1DModel(region+"_DiLepMass_Weight_MuonSFUp", ";Dilepton Mass [GeV];Events / 5 GeV", nBins, min_mass, max_mass), "DiLepMass", "Weight_MuonSFUp"))
            hists.append(d.Histo1D(RDF.TH1DModel(region+"_DiLepMass_Weight_ElectronSFUp", ";Dilepton Mass [GeV];Events / 5 GeV", nBins, min_mass, max_mass), "DiLepMass", "Weight_ElectronSFUp"))
            hists.append(d.Histo2D(RDF.TH2DModel(region+"_DiLepMass_bjetpt", ";Dilepton Mass [GeV];b jet p_T [GeV]", nBins, min_mass, max_mass, 100, 20, 400), "DiLepMass","leading_b_jet_pt", "Weight"))
            hists.append(d.Histo2D(RDF.TH2DModel(region+"_DiLepMass_jetpt", ";Dilepton Mass [GeV]; jet p_T [GeV]", nBins, min_mass, max_mass, 100, 20, 400), "DiLepMass","leading_jet_pt", "Weight"))
           
            if not is_inclusive:
              hists.append(d.Histo1D(RDF.TH1DModel(region+"_DiLepMass_Weight_BTagDown", ";Dilepton Mass [GeV];Events / 5 GeV", nBins, min_mass, max_mass), "DiLepMass", "Weight_BTagDown"))
              hists.append(d.Histo1D(RDF.TH1DModel(region+"_DiLepMass_Weight_PUIDDown", ";Dilepton Mass [GeV];Events / 5 GeV", nBins, min_mass, max_mass), "DiLepMass", "Weight_PUIDDown"))

              hists.append(d.Histo1D(RDF.TH1DModel(region+"_DiLepMass_Weight_BTagUp", ";Dilepton Mass [GeV];Events / 5 GeV", nBins, min_mass, max_mass), "DiLepMass", "Weight_BTagUp"))
              hists.append(d.Histo1D(RDF.TH1DModel(region+"_DiLepMass_Weight_PUIDUp", ";Dilepton Mass [GeV];Events / 5 GeV", nBins, min_mass, max_mass), "DiLepMass", "Weight_PUIDUp"))
              
            if is_quick: continue

            if is_inclusive:
              hists.append(d.Histo1D(RDF.TH1DModel(region+"_DiLepPt", ";Dilepton p_{T} [GeV];Events / 5 GeV", 100, 0, 500), "DiLepPt", "Weight"))
              underz = d.Filter('DiLepMass>70 && DiLepMass<105')
              hists.append(underz.Histo1D(RDF.TH1DModel(region+"_DiLepPtUnderZ", ";Dilepton p_{T} [GeV];Events / 5 GeV", 100, 0, 500), "DiLepPt", "Weight"))
            hists.append(d.Histo1D(RDF.TH1DModel(region+"_nPVs", ";# Pileup Vertices;Events", 101, -0.5, 100.5), "PV_npvs", "Weight"))
            hists.append(d.Histo1D(RDF.TH1DModel(region+"_nPVs_PuUp", ";# Pileup Vertices;Events", 101, -0.5, 100.5), "PV_npvs", "Weight_PuUp"))
            hists.append(d.Histo1D(RDF.TH1DModel(region+"_nPVs_PuDown", ";# Pileup Vertices;Events", 101, -0.5, 100.5), "PV_npvs", "Weight_PuDown"))
            hists.append(d.Histo1D(RDF.TH1DModel(region+"_Weight", ";Event Weight;Events", 100, -10, 10), "Weight"))
            hists.append(d.Histo1D(RDF.TH1DModel(region+"_MuonPt", ";Muon p_{T} [GeV];Events / 5 GeV", 100, 0, 500), "GoodMuonPtLow", "Weight"))
            hists.append(d.Histo1D(RDF.TH1DModel(region+"_MuonEta", ";Muon #eta;Events / 5 GeV", 80, -4, 4), "GoodMuonEtaLow", "Weight"))
            hists.append(d.Histo1D(RDF.TH1DModel(region+"_MuonPhi", ";Muon #phi;Events / 5 GeV", 128, -3.2, 3.2), "GoodMuonPhi", "Weight"))
            hists.append(d.Histo1D(RDF.TH1DModel(region+"_nMuon", ";# Muon;Events", 15, -0.5, 14.5), "nMuons", "Weight"))
            hists.append(d.Histo1D(RDF.TH1DModel(region+"_ElectronPt", ";Electron p_{T} [GeV];Events / 5 GeV", 100, 0, 500), "GoodElePtLow", "Weight"))
            hists.append(d.Histo1D(RDF.TH1DModel(region+"_ElectronEta", ";Electron #eta;Events / 5 GeV", 80, -4, 4), "GoodEleEtaLow", "Weight"))
            hists.append(d.Histo1D(RDF.TH1DModel(region+"_nElectron", ";# Electron;Events", 15, -0.5, 14.5), "nEle", "Weight"))
            hists.append(d.Histo1D(RDF.TH1DModel(region+"_JetPt", ";Jet p_{T} [GeV];Events / 5 GeV", 100, 0, 500), "GoodJetPt", "Weight"))
            hists.append(d.Histo1D(RDF.TH1DModel(region+"_JetEta", ";Jet #eta;Events / 5 GeV", 80, -4, 4), "GoodJetEta", "Weight"))


            if not is_inclusive:
              hists.append(d.Histo1D(RDF.TH1DModel(region+"_JetEta_bTagUp", ";Jet #eta;Events / 5 GeV", 80, -4, 4), "GoodJetEta", "Weight_BTagUp"))
              hists.append(d.Histo1D(RDF.TH1DModel(region+"_JetEta_bTagDown", ";Jet #eta;Events / 5 GeV", 80, -4, 4), "GoodJetEta", "Weight_BTagDown"))
              hists.append(d.Histo1D(RDF.TH1DModel(region+"_JetPt_bTagUp", ";Jet p_{T} [GeV];Events / 5 GeV", 100, 0, 500), "GoodJetPt", "Weight_BTagUp"))
              hists.append(d.Histo1D(RDF.TH1DModel(region+"_JetPt_bTagDown", ";Jet p_{T} [GeV];Events / 5 GeV", 100, 0, 500), "GoodJetPt", "Weight_BTagDown"))
            hists.append(d.Histo1D(RDF.TH1DModel(region+"_nJet", ";# Jet;Events", 15, -0.5, 14.5), "nJets", "Weight"))
            #systematic weight distributions for error-finding purposes

        outdir.cd()
        if not is_quick:
          total_Weight_h = next(h for h in hists if h.GetName() == "total_Weight")
          if is_inclusive:
            IncMuMu_DiLepMass_h = next(h for h in hists if h.GetName() == "IncMuMu_DiLepMass")
            IncEe_DiLepMass_h = next(h for h in hists if h.GetName() == "IncEe_DiLepMass")
            IncEmu_DiLepMass_h = next(h for h in hists if h.GetName() == "IncEmu_DiLepMass")
            print("{}\t ee: {}\t mumu: {}\t emu: {}".format(
                total_Weight_h.Integral(),
                IncEe_DiLepMass_h.Integral(),
                IncMuMu_DiLepMass_h.Integral(),
                IncEmu_DiLepMass_h.Integral(),
              ))
          else:
            print("total events: {}".format(total_Weight_h.Integral()))
          
        for h in hists:
          h.Write()
    out.Close()


if __name__ == '__main__':
    import argparse
    import sys
    parser = argparse.ArgumentParser(description="Produce histograms from BFF Z' reduced nano samples")
    parser.add_argument('samples', type=str, help='samples.txt file containing info on samples to run. See README for format')
    parser.add_argument('output', type=str, help='filename for the output ROOT file')
    parser.add_argument('--inclusive', '-i', default=False, action='store_true', help='Run the inclusive selections (default is running control+signal region selections)')
    parser.add_argument('--spark', '-s', default=False, action='store_true', help='Run on a spark cluster (needs PyRDF)')
    parser.add_argument('--npartitions', '-n', type=int, default=32, help='Number of data partitions (for spark)')
    parser.add_argument('--quick', '-q', default=False, action='store_true', help='Run reduced selection for faster protyping)')
    if len(sys.argv) < 3:
        parser.print_help()
        exit(123)
    args = parser.parse_args()
    process_samples(args.samples, args.output, args.inclusive, args.spark, args.npartitions, args.quick)

