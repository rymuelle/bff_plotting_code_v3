columns_data = [
#'HEM_event_veto',
#'EE_L1_prefire_test',
'Weight',
'sample_weight',
'TriggerWeight',
'Flag_goodVertices',
'Flag_globalSuperTightHalo2016Filter',
'Flag_HBHENoiseFilter',
'Flag_HBHENoiseIsoFilter',
'Flag_EcalDeadCellTriggerPrimitiveFilter',
'Flag_BadPFMuonFilter',
'Flag_eeBadScFilter',
'Flag_METFilters',
]

variation_data = [
'minGoodJetElDR{}',
'minGoodJetMuDR{}',
'DiLepMass{}',
'HTLT{}',
'RelMET{}',
'TMB{}',
'TMBMin{}',
'TMBMax{}',
'SR2{}',
'SR1{}',
'CR10{}',
'CR11{}',
'CR12{}',
'CR13{}',
'CR14{}',
'SR2{}',
'CR20{}',
'CR21{}',
'CR22{}',
'CR23{}',
'CR24{}',
'nLep{}',
'nLowPtLep{}',
]

columns_data += [v.format('_jet_nom_muon_corrected_pt_ele_pt') for v in variation_data]


variations = [
"_jet_jesTotal{}_muon_corrected_pt_ele_pt",
"_jet_nom_muon_corrected{}_pt_ele_pt",
"_jet_jer{}_muon_corrected_pt_ele_pt",
"_jet_jesTotal{}_muon_corrected_pt_ele_pt",
]

def make_columns(era, columns_data):
    '''put in era as string to add on additional columns to save'''
    columns_mc = [
        "Weight_PuUp",
        "Weight_PuDown",
        "Weight_BTagCorrUp",
        "Weight_BTagUncorrUp",
        "Weight_BTagCorrDown",
        "Weight_BTagUncorrDown",
        "Weight_PUIDUp",
        "Weight_PUIDDown",
        "Weight_PDF_Up",
        "Weight_PDF_Down",
        "Weight_ISRFSR_Up",
        "Weight_ISRFSR_Down",
        "Weight_MuonSFUp",
        "Weight_MuonSFDown",
        "Weight_ElectronSFUp",
        "Weight_ElectronSFDown",    
        ]
    columns_mc += [c for c in columns_data]
    if era=="2016":
        columns_mc += ["Weight_L1Down","Weight_L1Up"]
        columns_mc += ["HLT_Mu50","HLT_TkMu50", "HLT_DoubleEle33_CaloIdL_MW"]
        columns_data += ["HLT_Mu50","HLT_TkMu50", "HLT_DoubleEle33_CaloIdL_MW"]
    if era=="2017":
        columns_mc += ["Weight_L1Down","Weight_L1Up"]
        columns_mc += ["HLT_Mu50","HLT_OldMu100", "HLT_TkMu100",
                       "HLT_DoubleEle33_CaloIdL_MW", "HLT_DoubleEle25_CaloIdL_MW"]
        columns_data += ["HLT_Mu50","HLT_OldMu100", "HLT_TkMu100",
                       "HLT_DoubleEle33_CaloIdL_MW", "HLT_DoubleEle25_CaloIdL_MW"]
    if era=="2018":
        variations.append("_jet_jesHEMIssue{}_muon_corrected_pt_ele_pt")
        columns_mc += ["HLT_Mu50","HLT_OldMu100", "HLT_TkMu100",
                       "HLT_DoubleEle25_CaloIdL_MW"]
        columns_data += ["HLT_Mu50","HLT_OldMu100", "HLT_TkMu100",
                       "HLT_DoubleEle25_CaloIdL_MW"]
    var_postfix = []
    for var_template in variations:
        for direction in ('Up', 'Down'):
            var = var_template.format(direction)
            var_postfix.append(var)
            columns_mc +=  [v.format(var) for v in variation_data]
    return columns_data, columns_mc, var_postfix
