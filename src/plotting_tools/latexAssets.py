abseta	 = "$\abs{\eta}$"
absetacut	 = "$\abseta < 2.4$"
ptcut	 = "$\pt > 53\GeV$"
dbse	 = "\delta_{\cPqb\cPqs}"
gbe	 = "g_\cPqb"
gmue	 = "g_{\mu}"
dbs	 = "$\dbse$"
gb	 = "$\gbe$"
gmu	 = "$\gmue$"
sigmaZP	 = "$\sigma_{\Z'}$"
ptmiss = 'p_{\\textrm{T}}^{miss}'
METmm	 = "${}/m_{{\mu\mu}}$".format(ptmiss)
Lagr	 = "\mathcal{L}"
SR	 = "$SR_b^{\mu\mu}$"
CRmmj	 = "$CR_j^{\mu\mu}$"
CReeb	 = "$CR_b^{ee}$"
CReej	 = "$CR_j^{ee}$"
SRTwo	 = "$SR_{b+j/b}^{\mu\mu}$"
CRmmjTwo	 = "$CR_{2j}^{\mu\mu}$"
CReebTwo	 = "$CR_{b+j/b}^{ee}$"
CReejTwo	 = "$CR_{2j}^{ee}$"
HTLT	 = "$H_{\\textrm{T}}-L_{\\textrm{T}}$"
RelMET	 = METmm
TMB	 = "$TMB$"
dr	 = "$\delta_{r}$"
zpm	 = "$m_{\PZpr}$"
mll = "$m_{\ell\ell}$"

region_str_to_latex = {"CR10": CRmmj,
"CR14": CReej,
"CR24": CReejTwo,
"CR20": CRmmjTwo,
"CR13": CReeb,
"CR23": CReebTwo,
"SR1": SR,
"SR2": SRTwo}


signal_type_dict = {
    "DB": "VV",
    "ST": "single t",
    "TT": '$t\\bar{t}$',
    "DY": "DY+Jets"
}