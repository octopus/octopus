! Note: this file is generated automatically by build/mk_kfunctionals_list.pl
!
!%Variable TnaddFunctional
!%Type integer
!%Section Hamiltonian
!%Description
!% Defines the Kinetic Functional to be used in a Subsystem calculation,
!% For more information on the functionals, see
!% <a href=http://www.tddft.org/programs/octopus/wiki/index.php/Libxc:manual#Available_functionals>
!% Libxc documentation</a>. The list provided here is from libxc 2.2.1; if you have
!% linked against a different libxc version, you may have a somewhat different set
!% of available functionals.
!% <br>Default: <tt>none</tt>
!%Option lda_k_tf               50
!% Thomas-Fermi kinetic energy functional
!%Option lda_k_lp               51
!% Lee and Parr Gaussian ansatz
!%Option gga_k_tfvw               52
!% Thomas-Fermi plus von Weiszaecker correction
!%Option gga_k_revapbeint               53
!% interpolated version of REVAPBE
!%Option gga_k_apbeint               54
!% interpolated version of APBE
!%Option gga_k_revapbe               55
!% revised APBE
!%Option gga_k_meyer               57
!% Meyer,  Wang, and Young
!%Option gga_k_apbe               185
!% mu fixed from the semiclassical neutral atom
!%Option gga_k_tw1               187
!% Tran and Wesolowski set 1 (Table II)
!%Option gga_k_tw2               188
!% Tran and Wesolowski set 2 (Table II)
!%Option gga_k_tw3               189
!% Tran and Wesolowski set 3 (Table II)
!%Option gga_k_tw4               190
!% Tran and Wesolowski set 4 (Table II)
!%Option gga_k_vw               500
!% von Weiszaecker functional
!%Option gga_k_ge2               501
!% Second-order gradient expansion (l = 1/9)
!%Option gga_k_golden               502
!% TF-lambda-vW form by Golden (l = 13/45)
!%Option gga_k_yt65               503
!% TF-lambda-vW form by Yonei and Tomishima (l = 1/5)
!%Option gga_k_baltin               504
!% TF-lambda-vW form by Baltin (l = 5/9)
!%Option gga_k_lieb               505
!% TF-lambda-vW form by Lieb (l = 0.185909191)
!%Option gga_k_absp1               506
!% gamma-TFvW form by Acharya et al [g = 1 - 1.412/N^(1/3)]
!%Option gga_k_absp2               507
!% gamma-TFvW form by Acharya et al [g = 1 - 1.332/N^(1/3)]
!%Option gga_k_gr               508
!% gamma-TFvW form by Gázquez and Robles
!%Option gga_k_ludena               509
!% gamma-TFvW form by Ludeña
!%Option gga_k_gp85               510
!% gamma-TFvW form by Ghosh and Parr
!%Option gga_k_pearson               511
!% Pearson
!%Option gga_k_ol1               512
!% Ou-Yang and Levy v.1
!%Option gga_k_ol2               513
!% Ou-Yang and Levy v.2
!%Option gga_k_fr_b88               514
!% Fuentealba & Reyes (B88 version)
!%Option gga_k_fr_pw86               515
!% Fuentealba & Reyes (PW86 version)
!%Option gga_k_dk               516
!% DePristo and Kress
!%Option gga_k_perdew               517
!% Perdew
!%Option gga_k_vsk               518
!% Vitos, Skriver, and Kollar
!%Option gga_k_vjks               519
!% Vitos, Johansson, Kollar, and Skriver
!%Option gga_k_ernzerhof               520
!% Ernzerhof
!%Option gga_k_lc94               521
!% Lembarki & Chermette
!%Option gga_k_llp               522
!% Lee, Lee & Parr
!%Option gga_k_thakkar               523
!% Thakkar 1992
!%Option none                       0
!% Exchange and correlation set to zero (not from libxc).
!%End
