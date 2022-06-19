bck_dict = {
    "DB": ['mc_zz', 'mc_wz', 'mc_ww'],
    "ST": ['mc_santitop', 'mc_stop'],
    "TT": ['mc_ttbar'],
    "DY": ['ZToEE_M_120_200', 'ZToEE_M_200_400', 'ZToEE_M_400_800',
       'ZToEE_M_50_120', 'ZToEE_M_800_1400', 'ZToMuMu_M_120_200',
       'ZToMuMu_M_200_400', 'ZToMuMu_M_400_800', 'ZToMuMu_M_50_120',
       'ZToMuMu_M_800_1400'],
}

bck_list = [name for cat, background in bck_dict.items() for name in background]

bck_colors=['#ccffcc','#ffccff','#ffff99','#99ffff']