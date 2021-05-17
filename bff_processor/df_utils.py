import pandas as pd

def make_df(file_list):
    import pandas as pd
    df = pd.DataFrame()
    for file in file_list:
        df = pd.concat([df,pd.read_csv(file)])
    return df

def preselection(df):
    mass_cut = (df.DiLepMass > 105)*(df.DiLepMass < 800)
    flag = (df.Flag_goodVertices > 0)*(df.Flag_globalSuperTightHalo2016Filter > 0)*(df.Flag_HBHENoiseFilter > 0)*(df.Flag_HBHENoiseIsoFilter > 0)*(df.Flag_EcalDeadCellTriggerPrimitiveFilter > 0)*(df.Flag_BadPFMuonFilter > 0)*(df.Flag_eeBadScFilter > 0)*(df.Flag_METFilters > 0)
    return df[mass_cut & flag]

def isin(column):
    def wrapper(df):
        return df[df[column]==1]
    return wrapper