! Note: this file is generated automatically by build/mk_functionals_list.pl
!
!%Variable XCFunctional
!%Type integer
!%Section Hamiltonian::XC
!%Description
!% Defines the exchange and correlation functional to be used;
!% they should be specified as a sum of a correlation term and an
!% exchange term. Defaults:
!% <br> 1D: lda_x_1d + lda_c_1d_csc
!% <br> 2D: lda_x_2d + lda_c_2d_amgb
!% <br> 3D: lda_x + lda_c_pz_mod
!%Option lda_x               1
!% Exchange
!%Option lda_c_wigner               2000
!% Wigner parametrization
!%Option lda_c_rpa               3000
!% Random Phase Approximation
!%Option lda_c_hl               4000
!% Hedin & Lundqvist
!%Option lda_c_gl               5000
!% Gunnarson & Lundqvist
!%Option lda_c_xalpha               6000
!% Slater Xalpha
!%Option lda_c_vwn               7000
!% Vosko, Wilk, & Nussair (5)
!%Option lda_c_vwn_rpa               8000
!% Vosko, Wilk, & Nussair (RPA)
!%Option lda_c_pz               9000
!% Perdew & Zunger
!%Option lda_c_pz_mod               10000
!% Perdew & Zunger (Modified)
!%Option lda_c_ob_pz               11000
!% Ortiz & Ballone (PZ)
!%Option lda_c_pw               12000
!% Perdew & Wang
!%Option lda_c_pw_mod               13000
!% Perdew & Wang (Modified)
!%Option lda_c_ob_pw               14000
!% Ortiz & Ballone (PW)
!%Option lda_c_2d_amgb               15000
!% Attacalite et al
!%Option lda_c_2d_prm               16000
!% Pittalis, Rasanen & Marques correlation in 2D
!%Option lda_c_vbh               17000
!% von Barth & Hedin
!%Option lda_c_1d_csc               18000
!% Casula, Sorella, and Senatore 1D correlation
!%Option lda_x_2d               19
!% Exchange in 2D
!%Option lda_xc_teter93               20000
!% Teter 93 parametrization
!%Option lda_x_1d               21
!% Exchange in 1D
!%Option lda_c_ml1               22000
!% Modified LSD (version 1) of Proynov and Salahub
!%Option lda_c_ml2               23000
!% Modified LSD (version 2) of Proynov and Salahub
!%Option lda_c_gombas               24000
!% Gombas parametrization
!%Option lda_c_pw_rpa               25000
!% Perdew & Wang fit of the RPA
!%Option lda_c_1d_loos               26000
!% P-F Loos correlation LDA
!%Option lda_c_rc04               27000
!% Ragot-Cortona
!%Option lda_c_vwn_1               28000
!% Vosko, Wilk, & Nussair (1)
!%Option lda_c_vwn_2               29000
!% Vosko, Wilk, & Nussair (2)
!%Option lda_c_vwn_3               30000
!% Vosko, Wilk, & Nussair (3)
!%Option lda_c_vwn_4               31000
!% Vosko, Wilk, & Nussair (4)
!%Option gga_c_op_xalpha               84000
!% one-parameter progressive functional (G96 version)
!%Option gga_c_op_g96               85000
!% one-parameter progressive functional (G96 version)
!%Option gga_c_op_pbe               86000
!% one-parameter progressive functional (PBE version)
!%Option gga_c_op_b88               87000
!% one-parameter progressive functional (B88 version)
!%Option gga_c_ft97               88000
!% Filatov & Thiel correlation
!%Option gga_c_spbe               89000
!% PBE correlation to be used with the SSB exchange
!%Option gga_x_ssb_sw               90
!% Swarta, Sola and Bickelhaupt correction to PBE
!%Option gga_x_ssb               91
!% Swarta, Sola and Bickelhaupt
!%Option gga_x_ssb_d               92
!% Swarta, Sola and Bickelhaupt dispersion
!%Option gga_xc_hcth_407p               93000
!% HCTH/407+
!%Option gga_xc_hcth_p76               94000
!% HCTH p=7/6
!%Option gga_xc_hcth_p14               95000
!% HCTH p=1/4
!%Option gga_xc_b97_gga1               96000
!% Becke 97 GGA-1
!%Option gga_xc_hcth_a               97000
!% HCTH-A
!%Option gga_x_bpccac               98
!% BPCCAC (GRAC for the energy)
!%Option gga_c_revtca               99000
!% Tognetti, Cortona, Adamo (revised)
!%Option gga_c_tca               100000
!% Tognetti, Cortona, Adamo
!%Option gga_x_pbe               101
!% Perdew, Burke & Ernzerhof exchange
!%Option gga_x_pbe_r               102
!% Perdew, Burke & Ernzerhof exchange (revised)
!%Option gga_x_b86               103
!% Becke 86 Xalfa,beta,gamma
!%Option gga_x_herman               104
!% Herman et al original GGA
!%Option gga_x_b86_mgc               105
!% Becke 86 Xalfa,beta,gamma (with mod. grad. correction)
!%Option gga_x_b88               106
!% Becke 88
!%Option gga_x_g96               107
!% Gill 96
!%Option gga_x_pw86               108
!% Perdew & Wang 86
!%Option gga_x_pw91               109
!% Perdew & Wang 91
!%Option gga_x_optx               110
!% Handy & Cohen OPTX 01
!%Option gga_x_dk87_r1               111
!% dePristo & Kress 87 (version R1)
!%Option gga_x_dk87_r2               112
!% dePristo & Kress 87 (version R2)
!%Option gga_x_lg93               113
!% Lacks & Gordon 93
!%Option gga_x_ft97_a               114
!% Filatov & Thiel 97 (version A)
!%Option gga_x_ft97_b               115
!% Filatov & Thiel 97 (version B)
!%Option gga_x_pbe_sol               116
!% Perdew, Burke & Ernzerhof exchange (solids)
!%Option gga_x_rpbe               117
!% Hammer, Hansen & Norskov (PBE-like)
!%Option gga_x_wc               118
!% Wu & Cohen
!%Option gga_x_mpw91               119
!% Modified form of PW91 by Adamo & Barone
!%Option gga_x_am05               120
!% Armiento & Mattsson 05 exchange
!%Option gga_x_pbea               121
!% Madsen (PBE-like)
!%Option gga_x_mpbe               122
!% Adamo & Barone modification to PBE
!%Option gga_x_xpbe               123
!% xPBE reparametrization by Xu & Goddard
!%Option gga_x_2d_b86_mgc               124
!% Becke 86 MGC for 2D systems
!%Option gga_x_bayesian               125
!% Bayesian best fit for the enhancement factor
!%Option gga_x_pbe_jsjr               126
!% JSJR reparametrization by Pedroza, Silva & Capelle
!%Option gga_x_2d_b88               127
!% Becke 88 in 2D
!%Option gga_x_2d_b86               128
!% Becke 86 Xalfa,beta,gamma
!%Option gga_x_2d_pbe               129
!% Perdew, Burke & Ernzerhof exchange in 2D
!%Option gga_c_pbe               130000
!% Perdew, Burke & Ernzerhof correlation
!%Option gga_c_lyp               131000
!% Lee, Yang & Parr
!%Option gga_c_p86               132000
!% Perdew 86
!%Option gga_c_pbe_sol               133000
!% Perdew, Burke & Ernzerhof correlation SOL
!%Option gga_c_pw91               134000
!% Perdew & Wang 91
!%Option gga_c_am05               135000
!% Armiento & Mattsson 05 correlation
!%Option gga_c_xpbe               136000
!% xPBE reparametrization by Xu & Goddard
!%Option gga_c_lm               137000
!% Langreth and Mehl correlation
!%Option gga_c_pbe_jrgx               138000
!% JRGX reparametrization by Pedroza, Silva & Capelle
!%Option gga_x_optb88_vdw               139
!% Becke 88 reoptimized to be used with vdW functional of Dion et al
!%Option gga_x_pbek1_vdw               140
!% PBE reparametrization for vdW
!%Option gga_x_optpbe_vdw               141
!% PBE reparametrization for vdW
!%Option gga_x_rge2               142
!% Regularized PBE
!%Option gga_c_rge2               143000
!% Regularized PBE
!%Option gga_x_rpw86               144
!% refitted Perdew & Wang 86
!%Option gga_x_kt1               145
!% Keal and Tozer version 1
!%Option gga_xc_kt2               146000
!% Keal and Tozer version 2
!%Option gga_c_wl               147000
!% Wilson & Levy
!%Option gga_c_wi               148000
!% Wilson & Ivanov
!%Option gga_x_mb88               149
!% Modified Becke 88 for proton transfer
!%Option gga_x_sogga               150
!% Second-order generalized gradient approximation
!%Option gga_x_sogga11               151
!% Second-order generalized gradient approximation 2011
!%Option gga_c_sogga11               152000
!% Second-order generalized gradient approximation 2011
!%Option gga_c_wi0               153000
!% Wilson & Ivanov initial version
!%Option gga_xc_th1               154000
!% Tozer and Handy v. 1
!%Option gga_xc_th2               155000
!% Tozer and Handy v. 2
!%Option gga_xc_th3               156000
!% Tozer and Handy v. 3
!%Option gga_xc_th4               157000
!% Tozer and Handy v. 4
!%Option gga_x_c09x               158
!% C09x to be used with the VdW of Rutgers-Chalmers
!%Option gga_c_sogga11_x               159000
!% To be used with hyb_gga_x_SOGGA11-X
!%Option gga_x_lb               160
!% van Leeuwen & Baerends
!%Option gga_xc_hcth_93               161000
!% HCTH functional fitted to  93 molecules
!%Option gga_xc_hcth_120               162000
!% HCTH functional fitted to 120 molecules
!%Option gga_xc_hcth_147               163000
!% HCTH functional fitted to 147 molecules
!%Option gga_xc_hcth_407               164000
!% HCTH functional fitted to 407 molecules
!%Option gga_xc_edf1               165000
!% Empirical functionals from Adamson, Gill, and Pople
!%Option gga_xc_xlyp               166000
!% XLYP functional
!%Option gga_xc_b97               167000
!% Becke 97
!%Option gga_xc_b97_1               168000
!% Becke 97-1
!%Option gga_xc_b97_2               169000
!% Becke 97-2
!%Option gga_xc_b97_d               170000
!% Grimme functional to be used with C6 vdW term
!%Option gga_xc_b97_k               171000
!% Boese-Martin for Kinetics
!%Option gga_xc_b97_3               172000
!% Becke 97-3
!%Option gga_xc_pbe1w               173000
!% Functionals fitted for water
!%Option gga_xc_mpwlyp1w               174000
!% Functionals fitted for water
!%Option gga_xc_pbelyp1w               175000
!% Functionals fitted for water
!%Option gga_xc_sb98_1a               176000
!% Schmider-Becke 98 parameterization 1a
!%Option gga_xc_sb98_1b               177000
!% Schmider-Becke 98 parameterization 1b
!%Option gga_xc_sb98_1c               178000
!% Schmider-Becke 98 parameterization 1c
!%Option gga_xc_sb98_2a               179000
!% Schmider-Becke 98 parameterization 2a
!%Option gga_xc_sb98_2b               180000
!% Schmider-Becke 98 parameterization 2b
!%Option gga_xc_sb98_2c               181000
!% Schmider-Becke 98 parameterization 2c
!%Option gga_x_lbm               182
!% van Leeuwen & Baerends modified
!%Option gga_x_ol2               183
!% Exchange form based on Ou-Yang and Levy v.2
!%Option gga_x_apbe               184
!% mu fixed from the semiclassical neutral atom
!%Option gga_c_apbe               186000
!% mu fixed from the semiclassical neutral atom
!%Option gga_x_htbs               191
!% Haas, Tran, Blaha, and Schwarz
!%Option gga_x_airy               192
!% Constantin et al based on the Airy gas
!%Option gga_x_lag               193
!% Local Airy Gas
!%Option gga_xc_mohlyp               194000
!% Functional for organometallic chemistry
!%Option gga_xc_mohlyp2               195000
!% Functional for barrier heights
!%Option gga_xc_th_fl               196000
!% Tozer and Handy v. FL
!%Option gga_xc_th_fc               197000
!% Tozer and Handy v. FC
!%Option gga_xc_th_fcfo               198000
!% Tozer and Handy v. FCFO
!%Option gga_xc_th_fco               199000
!% Tozer and Handy v. FCO
!%Option gga_c_optc               200000
!% Optimized correlation functional of Cohen and Handy
!%Option gga_x_wpbeh               524
!% short-range version of the PBE
!%Option gga_x_hjs_pbe               525
!% HJS screened exchange PBE version
!%Option gga_x_hjs_pbe_sol               526
!% HJS screened exchange PBE_SOL version
!%Option gga_x_hjs_b88               527
!% HJS screened exchange B88 version
!%Option gga_x_hjs_b97x               528
!% HJS screened exchange B97x version
!%Option gga_x_ityh               529
!% short-range recipe for exchange GGA functionals
!%Option hyb_gga_xc_b3pw91               401000
!% The original (ACM) hybrid of Becke
!%Option hyb_gga_xc_b3lyp               402000
!% The (in)famous B3LYP
!%Option hyb_gga_xc_b3p86               403000
!% Perdew 86 hybrid similar to B3PW91
!%Option hyb_gga_xc_o3lyp               404000
!% hybrid using the optx functional
!%Option hyb_gga_xc_mpw1k               405000
!% mixture of mPW91 and PW91 optimized for kinetics
!%Option hyb_gga_xc_pbeh               406000
!% aka PBE0 or PBE1PBE
!%Option hyb_gga_xc_b97               407000
!% Becke 97
!%Option hyb_gga_xc_b97_1               408000
!% Becke 97-1
!%Option hyb_gga_xc_b97_2               410000
!% Becke 97-2
!%Option hyb_gga_xc_x3lyp               411000
!% maybe the best hybrid
!%Option hyb_gga_xc_b1wc               412000
!% Becke 1-parameter mixture of WC and PBE
!%Option hyb_gga_xc_b97_k               413000
!% Boese-Martin for Kinetics
!%Option hyb_gga_xc_b97_3               414000
!% Becke 97-3
!%Option hyb_gga_xc_mpw3pw               415000
!% mixture with the mPW functional
!%Option hyb_gga_xc_b1lyp               416000
!% Becke 1-parameter mixture of B88 and LYP
!%Option hyb_gga_xc_b1pw91               417000
!% Becke 1-parameter mixture of B88 and PW91
!%Option hyb_gga_xc_mpw1pw               418000
!% Becke 1-parameter mixture of mPW91 and PW91
!%Option hyb_gga_xc_mpw3lyp               419000
!% mixture of mPW and LYP
!%Option hyb_gga_xc_sb98_1a               420000
!% Schmider-Becke 98 parameterization 1a
!%Option hyb_gga_xc_sb98_1b               421000
!% Schmider-Becke 98 parameterization 1b
!%Option hyb_gga_xc_sb98_1c               422000
!% Schmider-Becke 98 parameterization 1c
!%Option hyb_gga_xc_sb98_2a               423000
!% Schmider-Becke 98 parameterization 2a
!%Option hyb_gga_xc_sb98_2b               424000
!% Schmider-Becke 98 parameterization 2b
!%Option hyb_gga_xc_sb98_2c               425000
!% Schmider-Becke 98 parameterization 2c
!%Option hyb_gga_x_sogga11_x               426
!% Hybrid based on SOGGA11 form
!%Option hyb_gga_xc_bhandh               435000
!% Becke half-and-half
!%Option hyb_gga_xc_bhandhlyp               436000
!% Becke half-and-half with B88 exchange
!%Option hyb_gga_xc_mb3lyp_rc04               437000
!% B3LYP with RC04 LDA
!%Option mgga_x_lta               201
!% Local tau approximation of Ernzerhof & Scuseria
!%Option mgga_x_tpss               202
!% Perdew, Tao, Staroverov & Scuseria exchange
!%Option mgga_x_m06_l               203
!% M06-Local functional of Minnesota
!%Option mgga_x_gvt4               204
!% GVT4 from Van Voorhis and Scuseria
!%Option mgga_x_tau_hcth               205
!% tau-HCTH from Boese and Handy
!%Option mgga_x_br89               206
!% Becke-Roussel 89
!%Option mgga_x_bj06               207
!% Becke & Johnson correction to Becke-Roussel 89
!%Option mgga_x_tb09               208
!% Tran & Blaha correction to Becke & Johnson
!%Option mgga_x_rpp09               209
!% Rasanen, Pittalis, and Proetto correction to Becke & Johnson
!%Option mgga_x_2d_prhg07               210
!% Pittalis, Rasanen, Helbig, Gross Exchange Functional
!%Option mgga_x_2d_prhg07_prp10               211
!% PRGH07 with PRP10 correction
!%Option mgga_x_revtpss               212
!% revised Perdew, Tao, Staroverov & Scuseria exchange
!%Option mgga_x_pkzb               213
!% Perdew, Kurth, Zupan, and Blaha
!%Option mgga_x_m05               214
!% M05 functional of Minnesota
!%Option mgga_x_m05_2x               215
!% M05-2X functional of Minnesota
!%Option mgga_x_m06_hf               216
!% M06-HF functional of Minnesota
!%Option mgga_x_m06               217
!% M06 functional of Minnesota
!%Option mgga_x_m06_2x               218
!% M06-2X functional of Minnesota
!%Option mgga_x_m08_hx               219
!% M08-HX functional of Minnesota
!%Option mgga_x_m08_so               220
!% M08-SO functional of Minnesota
!%Option mgga_c_tpss               231000
!% Perdew, Tao, Staroverov & Scuseria correlation
!%Option mgga_c_vsxc               232000
!% VSxc from Van Voorhis and Scuseria (correlation part)
!%Option mgga_c_m06_l               233000
!% M06-Local functional of Minnesota
!%Option mgga_c_m06_hf               234000
!% M06-HF functional of Minnesota
!%Option mgga_c_m06               235000
!% M06 functional of Minnesota
!%Option mgga_c_m06_2x               236000
!% M06-2X functional of Minnesota
!%Option mgga_c_m05               237000
!% M05 functional of Minnesota
!%Option mgga_c_m05_2x               238000
!% M05-2X functional of Minnesota
!%Option mgga_c_pkzb               239000
!% Perdew, Kurth, Zupan, and Blaha
!%Option mgga_c_bc95               240000
!% Becke correlation 95
!%Option hyb_mgga_xc_m05               438000
!% M05 functional of Minnesota
!%Option hyb_mgga_xc_m05_2x               439000
!% M05-2X functional of Minnesota
!%Option hyb_mgga_xc_b88b95               440000
!% Mixture of B88 with BC95 (B1B95)
!%Option hyb_mgga_xc_b86b95               441000
!% Mixture of B86 with BC95
!%Option hyb_mgga_xc_pw86b95               442000
!% Mixture of PW86 with BC95
!%Option hyb_mgga_xc_bb1k               443000
!% Mixture of B88 with BC95 from Zhao and Truhlar
!%Option hyb_mgga_xc_m06_hf               444000
!% M06-HF functional of Minnesota
!%Option hyb_mgga_xc_mpw1b95               445000
!% Mixture of mPW91 with BC95 from Zhao and Truhlar
!%Option hyb_mgga_xc_mpwb1k               446000
!% Mixture of mPW91 with BC95 for kinetics
!%Option hyb_mgga_xc_x1b95               447000
!% Mixture of X with BC95
!%Option hyb_mgga_xc_xb1k               448000
!% Mixture of X with BC95 for kinetics
!%Option hyb_mgga_xc_m06               449000
!% M06 functional of Minnesota
!%Option hyb_mgga_xc_m06_2x               450000
!% M06-2X functional of Minnesota
!%Option oep_x                    901
!% OEP: Exact exchange
!%Option lda_xc_cmplx             701
!% LDA complex scaled exchange-correlation.
!%Option rdmft_xc_m               601
!% RDMFT Mueller functional
!%Option none                       0
!% Exchange and correlation set to zero.
!%End
