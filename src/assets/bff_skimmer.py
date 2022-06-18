tmb_min = [
'TMBMin_nom',
'TMBMin_jesTotalUp',
'TMBMin_jesTotalDown',
'TMBMin_jerUp',
'TMBMin_jerDown']
tmb_max = [
'TMBMax_nom',
'TMBMax_jesTotalUp',
'TMBMax_jesTotalDown',
'TMBMax_jerUp',
'TMBMax_jerDown']

columns_data = [
#'HEM_event_veto',
'EE_L1_prefire_test',
'minGoodJetElDR',
'minGoodJetMuDR',
'Weight',
'sample_weight',
'DiLepMass',
'TriggerWeight',
'Flag_goodVertices',
'Flag_globalSuperTightHalo2016Filter',
'Flag_HBHENoiseFilter',
'Flag_HBHENoiseIsoFilter',
'Flag_EcalDeadCellTriggerPrimitiveFilter',
'Flag_BadPFMuonFilter',
'Flag_eeBadScFilter',
'Flag_METFilters',
'HTLT_nom',
'RelMET_nom',
'TMB_nom',
'SR2_nom',
'SR1_nom',
'CR10_nom',
'CR11_nom',
'CR12_nom',
'CR13_nom',
'CR14_nom',
'SR2_nom',
'CR20_nom',
'CR21_nom',
'CR22_nom',
'CR23_nom',
'CR24_nom',
'nLep',
'nLowPtLep',
]

columns_mc = [
"Weight_PuUp",
"Weight_PuDown",
"Weight_BTagUp",
"Weight_BTagDown",
"Weight_PUIDUp",
"Weight_PUIDDown",
"Weight_PDF_ISRFSR_Up",
"Weight_PDF_ISRFSR_Down",
"Weight_MuonSFUp",
"Weight_MuonSFDown",
"Weight_ElectronSFUp",
"Weight_ElectronSFDown",
"HTLT_jerDown",
"RelMET_jerDown",
"TMB_jerDown",
"SR2_jerDown",
"SR1_jerDown",
"CR10_jerDown",
"CR11_jerDown",
"CR12_jerDown",
"CR13_jerDown",
"CR14_jerDown",
"SR2_jerDown",
"CR20_jerDown",
"CR21_jerDown",
"CR22_jerDown",
"CR23_jerDown",
"CR24_jerDown",
"HTLT_jerUp",
"RelMET_jerUp",
"TMB_jerUp",
"SR2_jerUp",
"SR1_jerUp",
"CR10_jerUp",
"CR11_jerUp",
"CR12_jerUp",
"CR13_jerUp",
"CR14_jerUp",
"SR2_jerUp",
"CR20_jerUp",
"CR21_jerUp",
"CR22_jerUp",
"CR23_jerUp",
"CR24_jerUp",
"HTLT_jesTotalDown",
"RelMET_jesTotalDown",
"TMB_jesTotalDown",
"SR2_jesTotalDown",
"SR1_jesTotalDown",
"CR10_jesTotalDown",
"CR11_jesTotalDown",
"CR12_jesTotalDown",
"CR13_jesTotalDown",
"CR14_jesTotalDown",
"SR2_jesTotalDown",
"CR20_jesTotalDown",
"CR21_jesTotalDown",
"CR22_jesTotalDown",
"CR23_jesTotalDown",
"CR24_jesTotalDown",
"HTLT_jesTotalUp",
"RelMET_jesTotalUp",
"TMB_jesTotalUp",
"SR2_jesTotalUp",
"SR1_jesTotalUp",
"CR10_jesTotalUp",
"CR11_jesTotalUp",
"CR12_jesTotalUp",
"CR13_jesTotalUp",
"CR14_jesTotalUp",
"SR2_jesTotalUp",
"CR20_jesTotalUp",
"CR21_jesTotalUp",
"CR22_jesTotalUp",
"CR23_jesTotalUp",
"CR24_jesTotalUp",    
]
columns = columns_data + columns_mc

#define columns to save for hem root file
hem_columns = ['GoodJetEta','GoodJetPhi','GoodJetPt','GoodBJet','GoodMuonPt','GoodMuonEta','GoodMuonPhi','GoodElePt','GoodEleEta','GoodElePhi','DiLepMass','HTLT_nom','RelMET_nom','TMB_nom','SR2_nom','SR1_nom','CR10_nom','CR11_nom','CR12_nom','CR13_nom','CR14_nom','SR2_nom','CR20_nom','CR21_nom','CR22_nom','CR23_nom','CR24_nom','run',
              'MET_phi','MET_pt',]


columns+=tmb_min
columns+=tmb_max
columns_data+=['TMBMin_nom','TMBMax_nom']

def make_columns(era, columns):
    '''put in era as string to add on additional columns to save'''
    if era=="2016" or era=="2017":
        columns += ["Weight_L1Down","Weight_L1Up"]
    if era=="2018":
        columns += [
            "Weight_jesHEMDown"
            ,"Weight_jesHEMUp"
            ,"HTLT_jesHEMIssueUp"
            ,"HTLT_jesHEMIssueDown"
            ,"RelMET_jesHEMIssueUp"
            ,"RelMET_jesHEMIssueDown"
            ,"TMB_jesHEMIssueUp"
            ,"TMBMin_jesHEMIssueUp"
            ,"TMBMax_jesHEMIssueUp"
            ,"TMB_jesHEMIssueDown"
            ,"TMBMin_jesHEMIssueDown"
            ,"TMBMax_jesHEMIssueDown"
            ,"SR1_jesHEMIssueUp"
            ,"SR2_jesHEMIssueUp"
            ,"SR1_jesHEMIssueDown"
            ,"SR2_jesHEMIssueDown"
            ,"CR10_jesHEMIssueUp"
            ,"CR11_jesHEMIssueUp"
            ,"CR12_jesHEMIssueUp"
            ,"CR13_jesHEMIssueUp"
            ,"CR14_jesHEMIssueUp"
            ,"CR20_jesHEMIssueUp"
            ,"CR21_jesHEMIssueUp"
            ,"CR22_jesHEMIssueUp"
            ,"CR23_jesHEMIssueUp"
            ,"CR24_jesHEMIssueUp"
            ,"CR10_jesHEMIssueDown"
            ,"CR11_jesHEMIssueDown"
            ,"CR12_jesHEMIssueDown"
            ,"CR13_jesHEMIssueDown"
            ,"CR14_jesHEMIssueDown"
            ,"CR20_jesHEMIssueDown"
            ,"CR21_jesHEMIssueDown"
            ,"CR22_jesHEMIssueDown"
            ,"CR23_jesHEMIssueDown"
            ,"CR24_jesHEMIssueDown"
        ]
    return columns
