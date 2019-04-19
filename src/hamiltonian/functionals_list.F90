! Note: this file is generated automatically by build/mk_functionals_list.pl
!
!%Variable XCFunctional
!%Type integer
!%Section Hamiltonian::XC
!%Description
!% Defines the exchange and correlation functionals to be used,
!% specified as a sum of an exchange functional and a
!% correlation functional, or a single exchange-correlation functional
!% (<i>e.g.</i> <tt>hyb_gga_xc_pbeh</tt>). For more information on the functionals, see
!% <a href=https://gitlab.com/libxc/libxc/wikis/Functionals-list-4.0.4>
!% Libxc documentation</a>. The list provided here is from libxc 4; if you have
!% linked against a different libxc version, you may have a somewhat different set
!% of available functionals. Note that kinetic-energy functionals are not supported.
!%
!% The default functional will be selected by Octopus to be consistent
!% with the pseudopotentials you are using. If you are not using
!% pseudopotentials, Octopus cannot determine the functional used to
!% generate the pseudopotential, or the pseudopotential functionals
!% are inconsistent, Octopus will use the following defaults:
!%
!% <br>1D: <tt>lda_x_1d + lda_c_1d_csc</tt>
!% <br>2D: <tt>lda_x_2d + lda_c_2d_amgb</tt>
!% <br>3D: <tt>lda_x + lda_c_pz_mod</tt>
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
!% Vosko, Wilk, & Nusair (5)
!%Option lda_c_vwn_rpa               8000
!% Vosko, Wilk, & Nusair (RPA)
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
!% Attaccalite et al
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
!% Vosko, Wilk, & Nusair (1)
!%Option lda_c_vwn_2               29000
!% Vosko, Wilk, & Nusair (2)
!%Option lda_c_vwn_3               30000
!% Vosko, Wilk, & Nusair (3)
!%Option lda_c_vwn_4               31000
!% Vosko, Wilk, & Nusair (4)
!%Option lda_xc_zlp               43000
!% Zhao, Levy & Parr, Eq. (20)
!%Option lda_xc_ksdt               259000
!% Karasiev et al. parametrization
!%Option lda_c_chachiyo               287000
!% Chachiyo simple 2 parameter correlation
!%Option lda_c_lp96               289000
!% Liu-Parr correlation
!%Option lda_x_rel               532
!% Relativistic exchange
!%Option lda_xc_1d_ehwlrg_1               536000
!% LDA constructed from slab-like systems of 1 electron
!%Option lda_xc_1d_ehwlrg_2               537000
!% LDA constructed from slab-like systems of 2 electrons
!%Option lda_xc_1d_ehwlrg_3               538000
!% LDA constructed from slab-like systems of 3 electrons
!%Option lda_x_erf               546
!% Attenuated exchange LDA (erf)
!%Option lda_xc_lp_a               547000
!% Lee-Parr reparametrization B
!%Option lda_xc_lp_b               548000
!% Lee-Parr reparametrization B
!%Option lda_x_rae               549
!% Rae self-energy corrected exchange
!%Option lda_c_mcweeny               551000
!% McWeeny 76
!%Option lda_c_br78               552000
!% Brual & Rothstein 78
!%Option lda_c_pk09               554000
!% Proynov and Kong 2009
!%Option lda_c_ow_lyp               573000
!% Wigner with corresponding LYP parameters
!%Option lda_c_ow               574000
!% Optimized Wigner
!%Option lda_xc_gdsmfb               577000
!% Groth et al. parametrization
!%Option lda_c_gk72               578000
!% Gordon and Kim 1972
!%Option lda_c_karasiev               579000
!% Karasiev reparameterization of Chachiyo
!%Option gga_x_gam               32
!% GAM functional from Minnesota
!%Option gga_c_gam               33000
!% GAM functional from Minnesota
!%Option gga_x_hcth_a               34
!% HCTH-A
!%Option gga_x_ev93               35
!% Engel and Vosko
!%Option gga_x_bcgp               38
!% Burke, Cancio, Gould, and Pittalis
!%Option gga_c_bcgp               39000
!% Burke, Cancio, Gould, and Pittalis
!%Option gga_x_lambda_oc2_n               40
!% lambda_OC2(N) version of PBE
!%Option gga_x_b86_r               41
!% Revised Becke 86 Xalpha,beta,gamma (with mod. grad. correction)
!%Option gga_x_lambda_ch_n               44
!% lambda_CH(N) version of PBE
!%Option gga_x_lambda_lo_n               45
!% lambda_LO(N) version of PBE
!%Option gga_x_hjs_b88_v2               46
!% HJS screened exchange corrected B88 version
!%Option gga_c_q2d               47000
!% Chiodo et al
!%Option gga_x_q2d               48
!% Chiodo et al
!%Option gga_x_pbe_mol               49
!% Del Campo, Gazquez, Trickey and Vela (PBE-like)
!%Option gga_x_ak13               56
!% Armiento & Kuemmel 2013
!%Option gga_x_lv_rpw86               58
!% Berland and Hyldgaard
!%Option gga_x_pbe_tca               59
!% PBE revised by Tognetti et al
!%Option gga_x_pbeint               60
!% PBE for hybrid interfaces
!%Option gga_c_zpbeint               61000
!% spin-dependent gradient correction to PBEint
!%Option gga_c_pbeint               62000
!% PBE for hybrid interfaces
!%Option gga_c_zpbesol               63000
!% spin-dependent gradient correction to PBEsol
!%Option gga_xc_opbe_d               65000
!% oPBE_D functional of Goerigk and Grimme
!%Option gga_xc_opwlyp_d               66000
!% oPWLYP-D functional of Goerigk and Grimme
!%Option gga_xc_oblyp_d               67000
!% oBLYP-D functional of Goerigk and Grimme
!%Option gga_x_vmt84_ge               68
!% VMT{8,4} with constraint satisfaction with mu = mu_GE
!%Option gga_x_vmt84_pbe               69
!% VMT{8,4} with constraint satisfaction with mu = mu_PBE
!%Option gga_x_vmt_ge               70
!% Vela, Medel, and Trickey with mu = mu_GE
!%Option gga_x_vmt_pbe               71
!% Vela, Medel, and Trickey with mu = mu_PBE
!%Option gga_c_n12_sx               79000
!% N12-SX functional from Minnesota
!%Option gga_c_n12               80000
!% N12 functional from Minnesota
!%Option gga_x_n12               82
!% N12 functional from Minnesota
!%Option gga_c_regtpss               83000
!% Regularized TPSS correlation (ex-VPBE)
!%Option gga_c_op_xalpha               84000
!% one-parameter progressive functional (XALPHA version)
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
!% Swart, Sola and Bickelhaupt correction to PBE
!%Option gga_x_ssb               91
!% Swart, Sola and Bickelhaupt
!%Option gga_x_ssb_d               92
!% Swart, Sola and Bickelhaupt dispersion
!%Option gga_xc_hcth_407p               93000
!% HCTH/407+
!%Option gga_xc_hcth_p76               94000
!% HCTH p=7/6
!%Option gga_xc_hcth_p14               95000
!% HCTH p=1/4
!%Option gga_xc_b97_gga1               96000
!% Becke 97 GGA-1
!%Option gga_c_hcth_a               97000
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
!% Becke 86 Xalpha,beta,gamma
!%Option gga_x_herman               104
!% Herman et al original GGA
!%Option gga_x_b86_mgc               105
!% Becke 86 Xalpha,beta,gamma (with mod. grad. correction)
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
!% Becke 86 Xalpha,beta,gamma
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
!% Exchange part of Keal and Tozer version 1
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
!% To be used with HYB_GGA_X_SOGGA11_X
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
!%Option gga_xc_kt1               167000
!% Keal and Tozer version 1
!%Option gga_xc_b97_d               170000
!% Grimme functional to be used with C6 vdW term
!%Option gga_xc_pbe1w               173000
!% Functionals fitted for water
!%Option gga_xc_mpwlyp1w               174000
!% Functionals fitted for water
!%Option gga_xc_pbelyp1w               175000
!% Functionals fitted for water
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
!%Option gga_c_pbeloc               246000
!% Semilocal dynamical correlation
!%Option gga_xc_vv10               255000
!% Vydrov and Van Voorhis
!%Option gga_c_pbefe               258000
!% PBE for formation energies
!%Option gga_c_op_pw91               262000
!% one-parameter progressive functional (PW91 version)
!%Option gga_x_pbefe               265
!% PBE for formation energies
!%Option gga_x_cap               270
!% Correct Asymptotic Potential
!%Option gga_x_eb88               271
!% Non-empirical (excogitated) B88 functional of Becke and Elliott
!%Option gga_c_pbe_mol               272000
!% Del Campo, Gazquez, Trickey and Vela (PBE-like)
!%Option gga_c_bmk               280000
!% Boese-Martin for kinetics
!%Option gga_c_tau_hcth               281000
!% correlation part of tau-hcth
!%Option gga_c_hyb_tau_hcth               283000
!% correlation part of hyb_tau-hcth
!%Option gga_x_beefvdw               285
!% BEEF-vdW exchange
!%Option gga_xc_beefvdw               286000
!% BEEF-vdW exchange-correlation
!%Option gga_x_pbetrans               291
!% Gradient-based interpolation between PBE and revPBE
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
!%Option gga_x_sfat               530
!% short-range recipe for exchange GGA functionals
!%Option gga_x_sg4               533
!% Semiclassical GGA at fourth order
!%Option gga_c_sg4               534000
!% Semiclassical GGA at fourth order
!%Option gga_x_gg99               535
!% Gilbert and Gill 1999
!%Option gga_x_pbepow               539
!% PBE power
!%Option gga_x_kgg99               544
!% Gilbert and Gill 1999 (mixed)
!%Option gga_xc_hle16               545000
!% high local exchange 2016
!%Option gga_c_scan_e0               553000
!% GGA component of SCAN
!%Option gga_c_gapc               555000
!% GapC
!%Option gga_c_gaploc               556000
!% Gaploc
!%Option gga_c_zvpbeint               557000
!% another spin-dependent correction to PBEint
!%Option gga_c_zvpbesol               558000
!% another spin-dependent correction to PBEsol
!%Option gga_c_tm_lyp               559000
!% Takkar and McCarthy reparametrization
!%Option gga_c_tm_pbe               560000
!% Thakkar and McCarthy reparametrization
!%Option gga_c_w94               561000
!% Wilson 94 (Eq. 25)
!%Option gga_c_cs1               565000
!% A dynamical correlation functional
!%Option gga_x_b88m               570
!% Becke 88 reoptimized to be used with mgga_c_tau1
!%Option hyb_gga_x_n12_sx               81
!% N12-SX functional from Minnesota
!%Option hyb_gga_xc_b97_1p               266000
!% version of B97 by Cohen and Handy
!%Option hyb_gga_xc_pbe_mol0               273000
!% PBEmol0
!%Option hyb_gga_xc_pbe_sol0               274000
!% PBEsol0
!%Option hyb_gga_xc_pbeb0               275000
!% PBEbeta0
!%Option hyb_gga_xc_pbe_molb0               276000
!% PBEmolbeta0
!%Option hyb_gga_xc_pbe50               290000
!% PBE0 with 50% exx
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
!% hybrid by Xu and Goddard
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
!%Option hyb_gga_xc_hse03               427000
!% the 2003 version of the screened hybrid HSE
!%Option hyb_gga_xc_hse06               428000
!% the 2006 version of the screened hybrid HSE
!%Option hyb_gga_xc_hjs_pbe               429000
!% HJS hybrid screened exchange PBE version
!%Option hyb_gga_xc_hjs_pbe_sol               430000
!% HJS hybrid screened exchange PBE_SOL version
!%Option hyb_gga_xc_hjs_b88               431000
!% HJS hybrid screened exchange B88 version
!%Option hyb_gga_xc_hjs_b97x               432000
!% HJS hybrid screened exchange B97x version
!%Option hyb_gga_xc_cam_b3lyp               433000
!% CAM version of B3LYP
!%Option hyb_gga_xc_tuned_cam_b3lyp               434000
!% CAM version of B3LYP tuned for excitations
!%Option hyb_gga_xc_bhandh               435000
!% Becke half-and-half
!%Option hyb_gga_xc_bhandhlyp               436000
!% Becke half-and-half with B88 exchange
!%Option hyb_gga_xc_mb3lyp_rc04               437000
!% B3LYP with RC04 LDA
!%Option hyb_gga_xc_mpwlyp1m               453000
!% MPW with 1 par. for metals/LYP
!%Option hyb_gga_xc_revb3lyp               454000
!% Revised B3LYP
!%Option hyb_gga_xc_camy_blyp               455000
!% BLYP with yukawa screening
!%Option hyb_gga_xc_pbe0_13               456000
!% PBE0-1/3
!%Option hyb_gga_xc_b3lyps               459000
!% B3LYP* functional
!%Option hyb_gga_xc_wb97               463000
!% Chai and Head-Gordon
!%Option hyb_gga_xc_wb97x               464000
!% Chai and Head-Gordon
!%Option hyb_gga_xc_lrc_wpbeh               465000
!% Long-range corrected functional by Rorhdanz et al
!%Option hyb_gga_xc_wb97x_v               466000
!% Mardirossian and Head-Gordon
!%Option hyb_gga_xc_lcy_pbe               467000
!% PBE with yukawa screening
!%Option hyb_gga_xc_lcy_blyp               468000
!% BLYP with yukawa screening
!%Option hyb_gga_xc_lc_vv10               469000
!% Vydrov and Van Voorhis
!%Option hyb_gga_xc_camy_b3lyp               470000
!% B3LYP with Yukawa screening
!%Option hyb_gga_xc_wb97x_d               471000
!% Chai and Head-Gordon
!%Option hyb_gga_xc_hpbeint               472000
!% hPBEint
!%Option hyb_gga_xc_lrc_wpbe               473000
!% Long-range corrected functional by Rorhdanz et al
!%Option hyb_gga_xc_b3lyp5               475000
!% B3LYP with VWN functional 5 instead of RPA
!%Option hyb_gga_xc_edf2               476000
!% Empirical functional from Lin, George and Gill
!%Option hyb_gga_xc_cap0               477000
!% Correct Asymptotic Potential hybrid
!%Option hyb_gga_xc_lc_wpbe               478000
!% Long-range corrected functional by Vydrov and Scuseria
!%Option hyb_gga_xc_hse12               479000
!% HSE12 by Moussa, Schultz and Chelikowsky
!%Option hyb_gga_xc_hse12s               480000
!% Short-range HSE12 by Moussa, Schultz, and Chelikowsky
!%Option hyb_gga_xc_hse_sol               481000
!% HSEsol functional by Schimka, Harl, and Kresse
!%Option hyb_gga_xc_cam_qtp_01               482000
!% CAM-QTP(01): CAM-B3LYP retuned using ionization potentials of water
!%Option hyb_gga_xc_mpw1lyp               483000
!% Becke 1-parameter mixture of mPW91 and LYP
!%Option hyb_gga_xc_mpw1pbe               484000
!% Becke 1-parameter mixture of mPW91 and PBE
!%Option hyb_gga_xc_kmlyp               485000
!% Kang-Musgrave hybrid
!%Option hyb_gga_xc_b5050lyp               572000
!% Like B3LYP but more exact exchange
!%Option mgga_c_dldf               37000
!% Dispersionless Density Functional
!%Option mgga_xc_zlp               42000
!% Zhao, Levy & Parr, Eq. (21)
!%Option mgga_xc_otpss_d               64000
!% oTPSS_D functional of Goerigk and Grimme
!%Option mgga_c_cs               72000
!% Colle and Salvetti
!%Option mgga_c_mn12_sx               73000
!% MN12-SX correlation functional from Minnesota
!%Option mgga_c_mn12_l               74000
!% MN12-L correlation functional from Minnesota
!%Option mgga_c_m11_l               75000
!% M11-L correlation functional from Minnesota
!%Option mgga_c_m11               76000
!% M11 correlation functional from Minnesota
!%Option mgga_c_m08_so               77000
!% M08-SO correlation functional from Minnesota
!%Option mgga_c_m08_hx               78000
!% M08-HX correlation functional from Minnesota
!%Option mgga_x_lta               201
!% Local tau approximation of Ernzerhof & Scuseria
!%Option mgga_x_tpss               202
!% Tao, Perdew, Staroverov & Scuseria exchange
!%Option mgga_x_m06_l               203
!% M06-L exchange functional from Minnesota
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
!% revised Tao, Perdew, Staroverov & Scuseria exchange
!%Option mgga_x_pkzb               213
!% Perdew, Kurth, Zupan, and Blaha
!%Option mgga_x_ms0               221
!% MS exchange of Sun, Xiao, and Ruzsinszky
!%Option mgga_x_ms1               222
!% MS1 exchange of Sun, et al
!%Option mgga_x_ms2               223
!% MS2 exchange of Sun, et al
!%Option mgga_x_m11_l               226
!% M11-L exchange functional from Minnesota
!%Option mgga_x_mn12_l               227
!% MN12-L exchange functional from Minnesota
!%Option mgga_xc_cc06               229000
!% Cancio and Chou 2006
!%Option mgga_x_mk00               230
!% Exchange for accurate virtual orbital energies
!%Option mgga_c_tpss               231000
!% Tao, Perdew, Staroverov & Scuseria correlation
!%Option mgga_c_vsxc               232000
!% VSxc from Van Voorhis and Scuseria (correlation part)
!%Option mgga_c_m06_l               233000
!% M06-L correlation functional from Minnesota
!%Option mgga_c_m06_hf               234000
!% M06-HF correlation functional from Minnesota
!%Option mgga_c_m06               235000
!% M06 correlation functional from Minnesota
!%Option mgga_c_m06_2x               236000
!% M06-2X correlation functional from Minnesota
!%Option mgga_c_m05               237000
!% M05 correlation functional from Minnesota
!%Option mgga_c_m05_2x               238000
!% M05-2X correlation functional from Minnesota
!%Option mgga_c_pkzb               239000
!% Perdew, Kurth, Zupan, and Blaha
!%Option mgga_c_bc95               240000
!% Becke correlation 95
!%Option mgga_c_revtpss               241000
!% revised TPSS correlation
!%Option mgga_xc_tpsslyp1w               242000
!% Functionals fitted for water
!%Option mgga_x_mk00b               243
!% Exchange for accurate virtual orbital energies (v. B)
!%Option mgga_x_bloc               244
!% functional with balanced localization
!%Option mgga_x_modtpss               245
!% Modified Tao, Perdew, Staroverov & Scuseria exchange
!%Option mgga_c_tpssloc               247000
!% Semilocal dynamical correlation
!%Option mgga_x_mbeef               249
!% mBEEF exchange
!%Option mgga_x_mbeefvdw               250
!% mBEEF-vdW exchange
!%Option mgga_xc_b97m_v               254000
!% Mardirossian and Head-Gordon
!%Option mgga_x_mvs               257
!% MVS exchange of Sun, Perdew, and Ruzsinszky
!%Option mgga_x_mn15_l               260
!% MN15-L exhange functional from Minnesota
!%Option mgga_c_mn15_l               261000
!% MN15-L correlation functional from Minnesota
!%Option mgga_x_scan               263
!% SCAN exchange of Sun, Ruzsinszky, and Perdew
!%Option mgga_c_scan               267000
!% SCAN correlation
!%Option mgga_c_mn15               269000
!% MN15 correlation functional from Minnesota
!%Option mgga_x_b00               284
!% Becke 2000
!%Option mgga_xc_hle17               288000
!% high local exchange 2017
!%Option mgga_c_scan_rvv10               292000
!% SCAN correlation + rVV10 correlation
!%Option mgga_x_revm06_l               293
!% revised M06-L exchange functional from Minnesota
!%Option mgga_c_revm06_l               294000
!% Revised M06-L correlation functional from Minnesota
!%Option mgga_x_tm               540
!% Tao and Mo 2016
!%Option mgga_x_vt84               541
!% meta-GGA version of VT{8,4} GGA
!%Option mgga_x_sa_tpss               542
!% TPSS with correct surface asymptotics
!%Option mgga_c_kcis               562000
!% Krieger, Chen, Iafrate, and Savin
!%Option mgga_xc_lp90               564000
!% Lee & Parr, Eq. (56)
!%Option mgga_c_b88               571000
!% Meta-GGA correlation by Becke
!%Option mgga_x_gx               575
!% GX functional of Loos
!%Option mgga_x_pbe_gx               576
!% PBE-GX functional of Loos
!%Option hyb_mgga_x_dldf               36
!% Dispersionless Density Functional
!%Option hyb_mgga_x_ms2h               224
!% MS2 hybrid exchange of Sun, et al
!%Option hyb_mgga_x_mn12_sx               248
!% MN12-SX hybrid exchange functional from Minnesota
!%Option hyb_mgga_x_scan0               264
!% SCAN hybrid exchange
!%Option hyb_mgga_x_mn15               268
!% MN15 hybrid exchange functional from Minnesota
!%Option hyb_mgga_x_bmk               279
!% Boese-Martin for kinetics
!%Option hyb_mgga_x_tau_hcth               282
!% Hybrid version of tau-HCTH
!%Option hyb_mgga_x_m08_hx               295
!% M08-HX exchange functional from Minnesota
!%Option hyb_mgga_x_m08_so               296
!% M08-SO exchange functional from Minnesota
!%Option hyb_mgga_x_m11               297
!% M11 hybrid exchange functional from Minnesota
!%Option hyb_mgga_x_m05               438
!% M05 hybrid exchange functional from Minnesota
!%Option hyb_mgga_x_m05_2x               439
!% M05-2X hybrid exchange functional from Minnesota
!%Option hyb_mgga_xc_b88b95               440000
!% Mixture of B88 with BC95 (B1B95)
!%Option hyb_mgga_xc_b86b95               441000
!% Mixture of B86 with BC95
!%Option hyb_mgga_xc_pw86b95               442000
!% Mixture of PW86 with BC95
!%Option hyb_mgga_xc_bb1k               443000
!% Mixture of B88 with BC95 from Zhao and Truhlar
!%Option hyb_mgga_x_m06_hf               444
!% M06-HF hybrid exchange functional from Minnesota
!%Option hyb_mgga_xc_mpw1b95               445000
!% Mixture of mPW91 with BC95 from Zhao and Truhlar
!%Option hyb_mgga_xc_mpwb1k               446000
!% Mixture of mPW91 with BC95 for kinetics
!%Option hyb_mgga_xc_x1b95               447000
!% Mixture of X with BC95
!%Option hyb_mgga_xc_xb1k               448000
!% Mixture of X with BC95 for kinetics
!%Option hyb_mgga_x_m06               449
!% M06 hybrid exchange functional from Minnesota
!%Option hyb_mgga_x_m06_2x               450
!% M06-2X hybrid exchange functional from Minnesota
!%Option hyb_mgga_xc_pw6b95               451000
!% Mixture of PW91 with BC95 from Zhao and Truhlar
!%Option hyb_mgga_xc_pwb6k               452000
!% Mixture of PW91 with BC95 from Zhao and Truhlar for kinetics
!%Option hyb_mgga_xc_tpssh               457000
!% TPSS hybrid
!%Option hyb_mgga_xc_revtpssh               458000
!% revTPSS hybrid
!%Option hyb_mgga_x_mvsh               474
!% MVSh hybrid
!%Option hyb_mgga_xc_wb97m_v               531000
!% Mardirossian and Head-Gordon
!%Option hyb_mgga_xc_b0kcis               563000
!% Hybrid based on KCIS
!%Option hyb_mgga_xc_mpw1kcis               566000
!% Modified Perdew-Wang + KCIS hybrid
!%Option hyb_mgga_xc_mpwkcis1k               567000
!% Modified Perdew-Wang + KCIS hybrid with more exact exchange
!%Option hyb_mgga_xc_pbe1kcis               568000
!% Perdew-Burke-Ernzerhof + KCIS hybrid
!%Option hyb_mgga_xc_tpss1kcis               569000
!% TPSS hybrid with KCIS correlation
!%Option oep_x                    901
!% OEP: Exact exchange (not from libxc).
!%Option ks_inversion             801
!% Inversion of KS potential (not from libxc).
!%Option lda_xc_cmplx             701
!% Complex-scaled LDA exchange and correlation (not from libxc).
!%Option pbe_xc_cmplx             702
!% Complex-scaled PBE exchange and correlation (not from libxc).
!%Option lb94_xc_cmplx            703
!% Complex-scaled LB94 exchange and correlation (not from libxc).
!%Option rdmft_xc_m               601
!% RDMFT Mueller functional (not from libxc).
!%Option xc_half_hartree          917
!% Half-Hartree exchange for two electrons (supports complex scaling) (not from libxc).
!% Defined by <math>v_{xc}(r) = v_H(r) / 2</math>.
!%Option vdw_c_vdwdf      918000
!% van der Waals density functional vdW-DF correlation from libvdwxc (not from libxc).  Use with gga_x_pbe_r.
!%Option vdw_c_vdwdf2     919000
!% van der Waals density functional vdW-DF2 correlation from libvdwxc (not from libxc).  Use with gga_x_rpw86.
!%Option vdw_c_vdwdfcx    920000
!% van der Waals density functional vdW-DF-cx correlation from libvdwxc (not from libxc).  Use with gga_x_lv_rpw86.
!%Option none                       0
!% Exchange and correlation set to zero (not from libxc).
!%End
