combine_template = '''#higgs combine tool shape analysis card for z'to mumu 1 jet

-------------------------

imax 1  number of channels                                      #1 Jet
jmax 1  number of backgrounds -1                                    #following AN2015_207_v5, not sure why the -1 is there?
kmax *  number of nuisance parameters (sources of systematic uncertainties)

-------------------------

shapes data_obs * era_{era}_zp_ws.root zp:SR{nJets}_background  #dummy data for now

shapes abcd SR1     era_{era}_zp_ws.root zp:SR{nJets}_background          #parameterized ABCD pdf for modeling background

shapes BFFZp * era_{era}_zp_ws.root zp:{signame} zp:{signame}_$SYSTEMATIC

-------------------------

bin       SR1
observation   -1

-------------------------

bin       SR1       SR1
process     abcd    BFFZp
process     1     -1
rate      -1   -1

-------------------------
lumi lnN -      {lumi}

#delatB1 lnN       1.03648698929    - 
#deltaS1 lnN       -       1.00771086929


Weight_BTag    shape  -           1
Weight_PUID    shape  -           1
Weight_PDF_ISRFSR_ shape  -           1
Weight_MuonSF    shape  -           1
Weight_Pu      shape  -           1
Weight_ElectronSF  shape  -           1
jer  shape  -           1
jesTotal  shape  -           1'''