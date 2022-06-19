regions_template = ['SR{}', 'CR{}0', 'CR{}3', 'CR{}4']
nJets = [1,2]

regions_dict = {nJet:[temp.format(nJet) for temp in regions_template] for nJet in nJets}
regions = [temp.format(nJet) for temp in regions_template for nJet in nJets]


# This dict lists regions and provieds a latex string
region_and_label = {'SR1': "#mu#mu_{b}", 
           'CR10': "#mu#mu_{j}", 
           'CR13': "ee_{b}", 
           'CR14': "ee_{j}", 
           'SR2': "#mu#mu_{1,2 b}", 
           'CR20': "#mu#mu_{2 j}", 
           'CR23': "ee_{1,2 b}", 
           'CR24': "ee_{2 j}", 
          }
# This dict lists regions and provieds a latex string for AN
region_and_label_AN = {'SR1': "\SR", 
           'CR10': "\CRmmj", 
           'CR13': "\CReeb", 
           'CR14': "\CReej", 
           'SR2': "\SRTwo", 
           'CR20': "\CRmmjTwo", 
           'CR23': "\CReebTwo", 
           'CR24': "\CReejTwo", 
          }

