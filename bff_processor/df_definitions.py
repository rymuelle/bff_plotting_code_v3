from ROOT import RDataFrame, TFile
import ROOT
from glob import glob

def setup_btag_puid(ismc, era, bTagEff):
    '''read only one btag eff file'''
    bTagEff_files = glob(bTagEff)
    print(bTagEff_files[0])

    if ismc:
        bTagFile = TFile(bTagEff_files[0], 'read')
        ROOT.bTagEff = bTagFile.Get("bTagEff")
        PUIDSFfile = TFile("data/PUID_80XTraining_EffSFandUncties.root","read")
        ROOT.PUIDEff_true  = PUIDSFfile.Get("h2_eff_mc{}_T".format(era))
        ROOT.PUIDEff_false = PUIDSFfile.Get("h2_mistag_mc{}_T".format(era))
        ROOT.PUIDSF_true   = PUIDSFfile.Get("h2_eff_sf{}_T".format(era))
        ROOT.PUIDSF_false  = PUIDSFfile.Get("h2_mistag_sf{}_T".format(era))
        ROOT.PUIDUnc_true  = PUIDSFfile.Get("h2_eff_sf{}_T_Systuncty".format(era))
        ROOT.PUIDUnc_false = PUIDSFfile.Get("h2_mistag_sf{}_T_Systuncty".format(era))
    else: bTagFile, PUIDSFfile = None, None
    return bTagFile, PUIDSFfile

def def_good_jet(df, ismc, bDiscValue):
    df = df.Define("GoodJet", "Jet_pt>20 && abs(Jet_eta)<2.4 && (Jet_jetId > 3) && ((Jet_puId & 1) || Jet_pt>50.) && !(Jet_btagDeepB<={} && Jet_pt<=30.)".format(bDiscValue))\
           .Define("GoodJetPt", "Jet_pt[GoodJet]")\
           .Define("GoodJetEta", "Jet_eta[GoodJet]")\
           .Define("nJets", "GoodJetPt.size()")\
           .Define("BJet", "Jet_btagDeepB>{}".format(bDiscValue))\
           .Define("GoodBJet", "BJet[GoodJet]")
    df = df.Define("leading_b_jet_pt", "(GoodBJet.size()>0) ? Jet_pt[GoodBJet[0]] : 0")
    df = df.Define("leading_jet_pt", "(GoodJet.size()>0) ? Jet_pt[GoodJet[0]] : 0")
    df = df.Define("NoPUID_GoodJet", "Jet_pt>20 && abs(Jet_eta)<2.4 && (Jet_jetId & (1 << 1)) && !(Jet_btagDeepB<={} && Jet_pt<=30.) && Jet_pt<50.".format(bDiscValue))\
           .Define("NoPUID_GoodJetPt", "Jet_pt[NoPUID_GoodJet]")\
           .Define("NoPUID_GoodJetEta", "Jet_eta[NoPUID_GoodJet]")\
           .Define("NoPUID_PUID", "Jet_puId[NoPUID_GoodJet]")
    
    df = df.Define("minJetDRmuon", "min_delta_r(Jet_eta[GoodJet], Muon_eta[GoodMuonLowPt], Jet_phi[GoodJet], Muon_phi[GoodMuonLowPt])")
    df = df.Define("minJetDRelectron", "min_delta_r(Jet_eta[GoodJet], Electron_eta[GoodElectronLowPt], Jet_phi[GoodJet], Electron_phi[GoodElectronLowPt])")
    df = df.Define("minGoodJetMuDR", "min_vec(minJetDRmuon)")
    df = df.Define("minGoodJetElDR", "min_vec(minJetDRelectron)")
    if int(ismc):
        df = df.Define("NoPUID_GoodJetGenJetIdx", "Jet_genJetIdx[NoPUID_GoodJet]")
    return df

def def_good_leptons(df, ismc):
    df = df.Define("GoodMuon", "Muon_corrected_pt > 53 && abs(Muon_eta) < 2.4 && Muon_tightId > 0 && Muon_pfRelIso04_all < 0.25")\
           .Define("GoodMuonPt", "Muon_corrected_pt[GoodMuon]")\
           .Define("GoodMuonEta", "Muon_eta[GoodMuon]")\
           .Define("GoodMuonPhi", "Muon_phi[GoodMuon]")\
           .Define("GoodMuonCharge", "Muon_charge[GoodMuon]")\
           .Define("nMuons", "GoodMuonPt.size()")
    df = df.Define("GoodMuonLowPt", "Muon_corrected_pt > 24 && abs(Muon_eta) < 2.4 && Muon_tightId > 0 && Muon_pfRelIso04_all < 0.25")\
           .Define("GoodMuonPtLow", "Muon_corrected_pt[GoodMuonLowPt]")\
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
    return df
def def_HLT(df, ismc, era):
    if int(ismc):
        if era == "2016":
            df = df.Define("TriggerRegion1", "HLT_Mu50>0 or HLT_DoubleEle33_CaloIdL_MW>0")
            df = df.Define("TriggerRegion2", "HLT_Mu50>0 or HLT_TkMu50>0 or HLT_DoubleEle33_CaloIdL_MW>0")
            df = df.Define("TriggerRegion3", "HLT_Mu50>0 or HLT_TkMu50>0 or HLT_DoubleEle33_CaloIdL_MW>0 or HLT_DoubleEle33_CaloIdL_GsfTrkIdVL>0")
            df = df.Define("TriggerWeight", "TriggerRegion1*0.073+TriggerRegion2*(1.0-0.073-0.3087428721)+TriggerRegion3*0.3087428721")
        elif era == "2017":
            df = df.Define("TriggerRegion1", "HLT_Mu50>0 or HLT_OldMu100>0 or HLT_TkMu100>0 or HLT_DoubleEle33_CaloIdL_MW>0")
            df = df.Define("TriggerRegion2", "HLT_Mu50>0 or HLT_OldMu100>0 or HLT_TkMu100>0 or HLT_DoubleEle33_CaloIdL_MW>0 or HLT_DoubleEle25_CaloIdL_MW>0")
            df = df.Define("TriggerWeight", "TriggerRegion1*(1-0.6528268793)+TriggerRegion2*0.6528268793")
        elif era == "2018":
            df = df.Define("TriggerRegion1", "HLT_Mu50>0 or HLT_OldMu100>0 or HLT_TkMu100>0 or HLT_DoubleEle25_CaloIdL_MW>0")
            df = df.Define("TriggerWeight", "TriggerRegion1")
    else:
        df = df.Define("TriggerWeight", "1.0")
    return df
    
def def_sf_and_weight(df,ismc, is_inclusive, name,sample_weight):
    if int(ismc):
        if name not in ['mc_zz', 'mc_wz','mc_ww']:
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
            df = df.Define("k_factor","k_factor(2016, DiLepMass)")
        else: df = df.Define("k_factor","1")
        if not is_inclusive:
            df = df.Define("PUIDWeights", "GetPUIDweight(NoPUID_GoodJetPt, NoPUID_GoodJetEta, NoPUID_GoodJetGenJetIdx, NoPUID_PUID)")
            df = df.Define("PUIDWeight", "map_zero_to_one(PUIDWeights[0])")
            df = df.Define("PUIDWeightUp", "map_zero_to_one(PUIDWeights[1])")
            df = df.Define("PUIDWeightDown", "map_zero_to_one(PUIDWeights[2])")
            df = df.Define("GoodJetBTagSF", "(Jet_btagSF_deepcsv_M[GoodJet])")
            df = df.Define("GoodJetBTagSFUp", "(Jet_btagSF_deepcsv_M_up[GoodJet])")
            df = df.Define("GoodJetBTagSFDown", "(Jet_btagSF_deepcsv_M_down[GoodJet])")
            df = df.Define("GoodJetHadronFlav","Jet_hadronFlavour[GoodJet]")
            df = df.Define("BTagWeight", "map_zero_to_one(GetBTagWeight(GoodBJet, GoodJetHadronFlav, GoodJetPt, GoodJetBTagSF))")
            df = df.Define("BTagWeightUp", "map_zero_to_one(GetBTagWeight(GoodBJet, GoodJetHadronFlav, GoodJetPt, GoodJetBTagSFUp))")
            df = df.Define("BTagWeightDown", "map_zero_to_one(GetBTagWeight(GoodBJet, GoodJetHadronFlav, GoodJetPt, GoodJetBTagSFDown))")
            #event weights
            df = df.Define("sample_weight", str(sample_weight))
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
            df = df.Define("sample_weight", "1")
            df = df.Define("Weight_PuUp", "Weight/puWeight*puWeightUp")
            df = df.Define("Weight_PuDown", "Weight/puWeight*puWeightDown")
            df = df.Define("Weight_PDF_ISRFSR_Up", "Weight*PDF_ISRFSR_uncertaintyUp")
            df = df.Define("Weight_PDF_ISRFSR_Down", "Weight*PDF_ISRFSR_uncertaintyDown")
            df = df.Define("Weight_MuonSFUp", "Weight/MuonSFweight*MuonSFweightUp")
            df = df.Define("Weight_MuonSFDown", "Weight/MuonSFweight*MuonSFweightDown")
            df = df.Define("Weight_ElectronSFUp", "Weight/ElectronSFweight*ElectronSFweightUp")
            df = df.Define("Weight_ElectronSFDown", "Weight/ElectronSFweight*ElectronSFweightDown")
    else:
        df = df.Define("Weight", "1.")
        df = df.Define("sample_weight", "1")
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
    return df

def def_lep_selections(df):
    df = df.Define("oppositeSignMuon", "oppositeSign(GoodMuonCharge)")
    df = df.Define("oppositeSignElectron", "oppositeSign(GoodEleCharge)")
    df = df.Define("oppositeSignEmu", "oppositeSign(GoodMuonCharge,GoodEleCharge)")
    df = df.Define("IncMumuLocal", "nEle==0 && nMuons==2 && oppositeSignMuon==1")
    df = df.Define("IncMumuLocalLow", "IncMumuLocal && nMuonsLowPt==2 && nEleLowPt==0")
    df = df.Define("IncEeLocal", "nEle==2 && nMuons==0 && oppositeSignElectron==1")
    df = df.Define("IncEeLocalLow", "IncEeLocal && nEleLowPt==2 && nMuonsLowPt==0")
    df = df.Define("IncEmuLocal", "nEle==1 && nMuons==1 && oppositeSignEmu==1")
    df = df.Define("IncEmuLocalLow", "IncEmuLocal && nEleLowPt==1 && nMuonsLowPt==1")
    return df
