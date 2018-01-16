!> @file
!!   Filters for finite difference first derivative of a one dimensional function (order 16)
!! @author
!!    Copyright (C) 2014 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


  ! Setting of 'nord' order accurate first derivative coefficient from 'Matematica'.
  !Only nord=2,4,6,8,16

   ! Coefficients from Matematica program.
select case(nord)
case(2)

!In[292]:= f[x_]:=Evaluate[InterpolatingPolynomial[{{-1,ym1},{0,y0},{1,yp1}},x]]
!In[293]:= ff=FortranForm[Expand[Simplify[ReplaceAll[D[f[x],x],{x->0}]]]]
     c1D(-m:m,0) = (/&
     -1.d0/2.d0,0.d0,1.d0/2.d0&
     /)

!In[294]:= ff=FortranForm[Expand[Simplify[ReplaceAll[D[f[x],x],{x->-1}]]]]
     c1D(-m:m,-1) = (/&
     -3.d0/2.d0, 2.d0,-1.d0/2.d0&
     /)
!In[295]:= ff=FortranForm[Expand[Simplify[ReplaceAll[D[f[x],x],{x->1}]]]]
     c1D(-m:m,1) = (/&
     1.d0/2.d0,-2.d0,3.d0/2.d0&
     /)

   case(4)

!In[296]:= f[x_]:=Evaluate[InterpolatingPolynomial[{{-2,ym2},{-1,ym1},{0,y0},{1,yp1},{2,yp2}},x]]
!In[297]:= ff=FortranForm[Expand[Simplify[ReplaceAll[D[f[x],x],{x->0}]]]]


     c1D(-m:m,0) = (/&
     1.d0/12.d0,-2.d0/3.d0,0.d0,2.d0/3.d0,-1.d0/12.d0&
     /)
!In[298]:= ff=FortranForm[Expand[Simplify[ReplaceAll[D[f[x],x],{x->-1}]]]]
     c1D(-m:m,-1) = (/&
     -1.d0/4.d0,-5.d0/6.d0,3.d0/2.d0,-1.d0/2.d0,1.d0/12.d0&
     /)
!In[299]:= ff=FortranForm[Expand[Simplify[ReplaceAll[D[f[x],x],{x->-2}]]]]
     c1D(-m:m,-2) = (/&
     -25.d0/12.0,4.d0,-3.d0,4.d0/3.d0,-1.d0/4.d0&
     /)
!In[300]:= ff=FortranForm[Expand[Simplify[ReplaceAll[D[f[x],x],{x->1}]]]]
     c1D(-m:m,1) = (/&
     -1.d0/12.d0,+1.d0/2.d0,-3.d0/2.d0,+5.d0/6.d0,+1.d0/4.d0&
     /)
!In[301]:= ff=FortranForm[Expand[Simplify[ReplaceAll[D[f[x],x],{x->2}]]]]
     c1D(-m:m,2) = (/&
     1.d0/4.d0,-4.d0/3.d0,3.d0,-4.d0,25.d0/12.d0&
     /)

!In[302]:= f[x_]:=Evaluate[InterpolatingPolynomial[{{-3,ym3},{-2,ym2},{-1,ym1},{0,y0},{1,yp1},{2,yp2},{3,yp3}},x]]
!In[303]:= ff=FortranForm[Expand[Simplify[ReplaceAll[D[f[x],x],{x->0}]]]]

   case(6)

     c1D(-m:m,0) = (/&
     -1.d0/60.d0,3.d0/20.d0,-3.d0/4.d0,0.d0,3.d0/4.d0,-3.d0/20.d0,1.d0/60.d0&
     /)
!In[304]:= ff=FortranForm[Expand[Simplify[ReplaceAll[D[f[x],x],{x->-1}]]]]
     c1D(-m:m,-1) = (/&
     1.d0/30.d0,-2/5.d0,-7.d0/12.d0,4.d0/3.d0,-1.d0/2.d0,2.d0/15.d0,-1.d0/60.d0&
     /)
!In[305]:= ff=FortranForm[Expand[Simplify[ReplaceAll[D[f[x],x],{x->1}]]]]
     c1D(-m:m,1) = (/&
     1.d0/60.d0,-2.d0/15.d0,1.d0/2.d0,-4.d0/3.d0,7.d0/12.d0,2.d0/5.d0,-1.d0/30.d0&
     /)
!In[306]:= ff=FortranForm[Expand[Simplify[ReplaceAll[D[f[x],x],{x->-2}]]]]
     c1D(-m:m,-2) = (/&
     -1.d0/6.d0,-77.d0/60.d0,5.d0/2.d0,-5.d0/3.d0,5.d0/6.d0,-1.d0/4.d0,1.d0/30.d0&
     /)
!In[307]:= ff=FortranForm[Expand[Simplify[ReplaceAll[D[f[x],x],{x->+2}]]]]
     c1D(-m:m,2) = (/&
     -1.d0/30.d0,1.d0/4.d0,-5/6.d0,5.d0/3.d0,-5.d0/2.d0,77.d0/60.d0,1.d0/6.d0&
     /)
!In[308]:= ff=FortranForm[Expand[Simplify[ReplaceAll[D[f[x],x],{x->-3}]]]]
     c1D(-m:m,-3) = (/&
     -49.d0/20.d0,6.d0,-15.d0/2.d0,20.d0/3.d0,-15.d0/4.d0,6.d0/5.d0,-1.d0/6.d0&
     /)
!In[309]:= ff=FortranForm[Expand[Simplify[ReplaceAll[D[f[x],x],{x->+3}]]]]
     c1D(-m:m,3) = (/&
     1.d0/6.d0,-6.d0/5.d0,15.d0/4.d0,-20.d0/3.d0,15.d0/2.d0,-6.d0,49.d0/20.d0&
     /)
!In[310]:= f[x_]:=Evaluate[InterpolatingPolynomial[{{-4,ym4},{-3,ym3},{-2,ym2},{-1,ym1},{0,y0},{1,yp1},{2,yp2},{3,yp3},{4,yp4}},x]]
!In[311]:= ff=FortranForm[Expand[Simplify[ReplaceAll[D[f[x],x],{x->0}]]]]

   case(8)

     c1D(-m:m,0) = (/&
     1.d0/280.d0,-4.d0/105.d0,1.d0/5.d0,-4.d0/5.d0,0.d0,4.d0/5.d0,-1.d0/5.d0,4.d0/105.d0,-1.d0/280.d0&
     /)
!In[312]:= ff=FortranForm[Expand[Simplify[ReplaceAll[D[f[x],x],{x->-1}]]]]
     c1D(-m:m,-1) = (/&
     -1.d0/168.d0,1.d0/14.d0,-1.d0/2.d0,-9.d0/20.d0,5.d0/4.d0,-1.d0/2.d0,1.d0/6.d0,-1.d0/28.d0,1.d0/280.d0&
     /)
!In[313]:= ff=FortranForm[Expand[Simplify[ReplaceAll[D[f[x],x],{x->1}]]]]
     c1D(-m:m,1) = (/&
     -1.d0/280.d0,1.d0/28.d0,-1.d0/6.d0,1.d0/2.d0,-5.d0/4.d0,9.d0/20.d0,1.d0/2.d0,-1.d0/14.d0,1.d0/168.d0&
     /)
!In[314]:= ff=FortranForm[Expand[Simplify[ReplaceAll[D[f[x],x],{x->-2}]]]]
     c1D(-m:m,-2) = (/&
     1.d0/56.d0,-2.d0/7.d0,-19.d0/20.d0,2.d0,-5.d0/4.d0,2.d0/3.d0,-1.d0/4.d0,2.d0/35.d0,-1.d0/168.d0&
     /)
!In[315]:= ff=FortranForm[Expand[Simplify[ReplaceAll[D[f[x],x],{x->+2}]]]]
     c1D(-m:m,2) = (/&
      1.d0/168.d0,-2.d0/35.d0,1.d0/4.d0,-2.d0/3.d0,5.d0/4.d0,-2.d0,19.d0/20.d0,2.d0/7.d0,-1.d0/56.d0&
     /)
!In[316]:= ff=FortranForm[Expand[Simplify[ReplaceAll[D[f[x],x],{x->-3}]]]]
     c1D(-m:m,-3) = (/&
     -1.d0/8.d0,-223.d0/140.d0,7.d0/2.d0,-7.d0/2.d0,35.d0/12.d0,-7.d0/4.d0,7.d0/10.d0,& 
     -1.d0/6.d0,1.d0/56.d0&
     /)
!In[317]:= ff=FortranForm[Expand[Simplify[ReplaceAll[D[f[x],x],{x->+3}]]]]
     c1D(-m:m,3) = (/&
     -1.d0/56.d0,1.d0/6.d0,-7.d0/10.d0,7.d0/4.d0,-35.d0/12.d0,7.d0/2.d0,-7.d0/2.d0,&
     223.d0/140.d0,1.d0/8.d0&
     /)
!In[318]:= ff=FortranForm[Expand[Simplify[ReplaceAll[D[f[x],x],{x->-4}]]]]
     c1D(-m:m,-4) = (/&
     -761.d0/280.d0,8.d0,-14.d0,56.d0/3.d0,-35.d0/2.d0,56.d0/5.d0,-14.d0/3.d0,& 
     8.d0/7.d0,-1.d0/8.d0&
     /)
!In[319]:= ff=FortranForm[Expand[Simplify[ReplaceAll[D[f[x],x],{x->+4}]]]]
     c1D(-m:m,4) = (/&
     1.d0/8.d0,-8.d0/7.d0,14.d0/3.d0,-56.d0/5.d0,35.d0/2.d0,-56.d0/3.d0,14.d0,-8.d0,& 
     761.d0/280.d0&
     /)

   case(10)

!In[29]:= f[x_]:=Evaluate[InterpolatingPolynomial[{{-5,ym5},{-4,ym4},{-3,ym3},{-2,ym2},{-1,ym1},{0,y0},{1,yp1},{2,yp2},{3,yp3},{4,yp4},{5,yp5}},x]]
!In[30]:= 
!In[31]:= ff=FortranForm[Expand[Simplify[ReplaceAll[D[f[x],x],{x->0}]]]]
!Out[31]//FortranForm=
!        (-5*ym1)/6. + (5*ym2)/21. - (5*ym3)/84. + (5*ym4)/504. - ym5/1260. + (5*yp1)/6. - (5*yp2)/21. + 
!     -  (5*yp3)/84. - (5*yp4)/504. + yp5/1260.
!In[32]:= ff=FortranForm[Expand[Simplify[ReplaceAll[D[f[x],x],{x->-1}]]]]
!Out[32]//FortranForm=
!        (6*y0)/5. - (11*ym1)/30. - (4*ym2)/7. + (3*ym3)/28. - ym4/63. + ym5/840. - yp1/2. + 
!     -  (4*yp2)/21. - (3*yp3)/56. + yp4/105. - yp5/1260.
!In[33]:= ff=FortranForm[Expand[Simplify[ReplaceAll[D[f[x],x],{x->1}]]]]
!Out[33]//FortranForm=
!        (-6*y0)/5. + ym1/2. - (4*ym2)/21. + (3*ym3)/56. - ym4/105. + ym5/1260. + (11*yp1)/30. + 
!     -  (4*yp2)/7. - (3*yp3)/28. + yp4/63. - yp5/840.
!In[34]:= ff=FortranForm[Expand[Simplify[ReplaceAll[D[f[x],x],{x->-2}]]]]
!Out[34]//FortranForm=
!        (-21*y0)/20. + (7*ym1)/4. - (319*ym2)/420. - (3*ym3)/8. + ym4/24. - ym5/360. + (7*yp1)/12. - 
!     -  yp2/4. + (3*yp3)/40. - yp4/72. + yp5/840.
!In[35]:= ff=FortranForm[Expand[Simplify[ReplaceAll[D[f[x],x],{x->+2}]]]]
!Out[35]//FortranForm=
!        (21*y0)/20. - (7*ym1)/12. + ym2/4. - (3*ym3)/40. + ym4/72. - ym5/840. - (7*yp1)/4. + 
!     -  (319*yp2)/420. + (3*yp3)/8. - yp4/24. + yp5/360.
!In[36]:= ff=FortranForm[Expand[Simplify[ReplaceAll[D[f[x],x],{x->-3}]]]]
!Out[36]//FortranForm=
!        (28*y0)/15. - (7*ym1)/3. + (8*ym2)/3. - (341*ym3)/280. - (2*ym4)/9. + ym5/90. - (7*yp1)/6. + 
!     -  (8*yp2)/15. - yp3/6. + (2*yp4)/63. - yp5/360.
!In[37]:= ff=FortranForm[Expand[Simplify[ReplaceAll[D[f[x],x],{x->+3}]]]]
!Out[37]//FortranForm=
!        (-28*y0)/15. + (7*ym1)/6. - (8*ym2)/15. + ym3/6. - (2*ym4)/63. + ym5/360. + (7*yp1)/3. - 
!     -  (8*yp2)/3. + (341*yp3)/280. + (2*yp4)/9. - yp5/90.
!In[38]:= ff=FortranForm[Expand[Simplify[ReplaceAll[D[f[x],x],{x->-4}]]]]
!Out[38]//FortranForm=
!        (-63*y0)/10. + 7*ym1 - 6*ym2 + (9*ym3)/2. - (4609*ym4)/2520. - ym5/10. + (21*yp1)/5. - 2*yp2 + 
!     -  (9*yp3)/14. - yp4/8. + yp5/90.
!In[39]:= ff=FortranForm[Expand[Simplify[ReplaceAll[D[f[x],x],{x->+4}]]]]
!Out[39]//FortranForm=
!        (63*y0)/10. - (21*ym1)/5. + 2*ym2 - (9*ym3)/14. + ym4/8. - ym5/90. - 7*yp1 + 6*yp2 - 
!     -  (9*yp3)/2. + (4609*yp4)/2520. + yp5/10.
!In[40]:= ff=FortranForm[Expand[Simplify[ReplaceAll[D[f[x],x],{x->-5}]]]]
!Out[40]//FortranForm=
!        (252*y0)/5. - (105*ym1)/2. + 40*ym2 - (45*ym3)/2. + 10*ym4 - (7381*ym5)/2520. - 35*yp1 + 
!     -  (120*yp2)/7. - (45*yp3)/8. + (10*yp4)/9. - yp5/10.
!In[41]:= ff=FortranForm[Expand[Simplify[ReplaceAll[D[f[x],x],{x->+5}]]]]
!
!Out[41]//FortranForm=
!        (-252*y0)/5. + 35*ym1 - (120*ym2)/7. + (45*ym3)/8. - (10*ym4)/9. + ym5/10. + (105*yp1)/2. - 
!     -  40*yp2 + (45*yp3)/2. - 10*yp4 + (7381*yp5)/2520.

   case(12)

!In[42]:= f[x_]:=Evaluate[InterpolatingPolynomial[{{-6,ym6},{-5,ym5},{-4,ym4},{-3,ym3},{-2,ym2},{-1,ym1},{0,y0},{1,yp1},{2,yp2},{3,yp3},{4,yp4},{5,yp5},{6,yp6}},x]]
!In[43]:= ff=FortranForm[Expand[Simplify[ReplaceAll[D[f[x],x],{x->0}]]]]
!Out[43]//FortranForm=
!        (-6*ym1)/7. + (15*ym2)/56. - (5*ym3)/63. + ym4/56. - ym5/385. + ym6/5544. + (6*yp1)/7. - 
!     -  (15*yp2)/56. + (5*yp3)/63. - yp4/56. + yp5/385. - yp6/5544.
!In[44]:= ff=FortranForm[Expand[Simplify[ReplaceAll[D[f[x],x],{x->-1}]]]]
!Out[44]//FortranForm=
!        (7*y0)/6. - (13*ym1)/42. - (5*ym2)/8. + (5*ym3)/36. - ym4/36. + ym5/264. - ym6/3960. - yp1/2. + 
!     -  (5*yp2)/24. - (5*yp3)/72. + yp4/60. - yp5/396. + yp6/5544.
!In[45]:= ff=FortranForm[Expand[Simplify[ReplaceAll[D[f[x],x],{x->1}]]]]
!Out[45]//FortranForm=
!        (-7*y0)/6. + ym1/2. - (5*ym2)/24. + (5*ym3)/72. - ym4/60. + ym5/396. - ym6/5544. + 
!     -  (13*yp1)/42. + (5*yp2)/8. - (5*yp3)/36. + yp4/36. - yp5/264. + yp6/3960.
!In[46]:= ff=FortranForm[Expand[Simplify[ReplaceAll[D[f[x],x],{x->-2}]]]]
!Out[46]//FortranForm=
!        (-14*y0)/15. + (8*ym1)/5. - (533*ym2)/840. - (4*ym3)/9. + ym4/15. - (4*ym5)/495. + ym6/1980. + 
!     -  (8*yp1)/15. - yp2/4. + (4*yp3)/45. - yp4/45. + (4*yp5)/1155. - yp6/3960.
!In[47]:= ff=FortranForm[Expand[Simplify[ReplaceAll[D[f[x],x],{x->+2}]]]]
!Out[47]//FortranForm=
!        (14*y0)/15. - (8*ym1)/15. + ym2/4. - (4*ym3)/45. + ym4/45. - (4*ym5)/1155. + ym6/3960. - 
!     -  (8*yp1)/5. + (533*yp2)/840. + (4*yp3)/9. - yp4/15. + (4*yp5)/495. - yp6/1980.
!In[48]:= 
!In[49]:= ff=FortranForm[Expand[Simplify[ReplaceAll[D[f[x],x],{x->-3}]]]]
!Out[49]//FortranForm=
!        (7*y0)/5. - (9*ym1)/5. + (9*ym2)/4. - (2509*ym3)/2520. - (3*ym4)/10. + (3*ym5)/110. - 
!     -  ym6/660. - (9*yp1)/10. + (9*yp2)/20. - yp3/6. + (3*yp4)/70. - (3*yp5)/440. + yp6/1980.
!In[50]:= ff=FortranForm[Expand[Simplify[ReplaceAll[D[f[x],x],{x->+3}]]]]
!Out[50]//FortranForm=
!        (-7*y0)/5. + (9*ym1)/10. - (9*ym2)/20. + ym3/6. - (3*ym4)/70. + (3*ym5)/440. - ym6/1980. + 
!     -  (9*yp1)/5. - (9*yp2)/4. + (2509*yp3)/2520. + (3*yp4)/10. - (3*yp5)/110. + yp6/660.
!In[51]:= ff=FortranForm[Expand[Simplify[ReplaceAll[D[f[x],x],{x->-4}]]]]
!Out[51]//FortranForm=
!        (-7*y0)/2. + 4*ym1 - (15*ym2)/4. + (10*ym3)/3. - (3601*ym4)/2520. - (2*ym5)/11. + ym6/132. + 
!     -  (12*yp1)/5. - (5*yp2)/4. + (10*yp3)/21. - yp4/8. + (2*yp5)/99. - yp6/660.
!In[52]:= 
!In[53]:= ff=FortranForm[Expand[Simplify[ReplaceAll[D[f[x],x],{x->+4}]]]]
!Out[53]//FortranForm=
!        (7*y0)/2. - (12*ym1)/5. + (5*ym2)/4. - (10*ym3)/21. + ym4/8. - (2*ym5)/99. + ym6/660. - 4*yp1 + 
!     -  (15*yp2)/4. - (10*yp3)/3. + (3601*yp4)/2520. + (2*yp5)/11. - yp6/132.
!In[54]:= ff=FortranForm[Expand[Simplify[ReplaceAll[D[f[x],x],{x->-5}]]]]
!Out[54]//FortranForm=
!        (77*y0)/5. - (33*ym1)/2. + (55*ym2)/4. - (55*ym3)/6. + (11*ym4)/2. - (55991*ym5)/27720. - 
!     -  ym6/12. - 11*yp1 + (165*yp2)/28. - (55*yp3)/24. + (11*yp4)/18. - yp5/10. + yp6/132.
!In[55]:= 
!In[56]:= ff=FortranForm[Expand[Simplify[ReplaceAll[D[f[x],x],{x->+5}]]]]
!Out[56]//FortranForm=
!        (-77*y0)/5. + 11*ym1 - (165*ym2)/28. + (55*ym3)/24. - (11*ym4)/18. + ym5/10. - ym6/132. + 
!     -  (33*yp1)/2. - (55*yp2)/4. + (55*yp3)/6. - (11*yp4)/2. + (55991*yp5)/27720. + yp6/12.
!In[57]:= ff=FortranForm[Expand[Simplify[ReplaceAll[D[f[x],x],{x->-6}]]]]
!Out[57]//FortranForm=
!        -154*y0 + (792*ym1)/5. - (495*ym2)/4. + (220*ym3)/3. - 33*ym4 + 12*ym5 - (86021*ym6)/27720. + 
!     -  (792*yp1)/7. - (495*yp2)/8. + (220*yp3)/9. - (33*yp4)/5. + (12*yp5)/11. - yp6/12.
!In[58]:= ff=FortranForm[Expand[Simplify[ReplaceAll[D[f[x],x],{x->+6}]]]]
!Out[58]//FortranForm=
!        154*y0 - (792*ym1)/7. + (495*ym2)/8. - (220*ym3)/9. + (33*ym4)/5. - (12*ym5)/11. + ym6/12. - 
!     -  (792*yp1)/5. + (495*yp2)/4. - (220*yp3)/3. + 33*yp4 - 12*yp5 + (86021*yp6)/27720.

   case(14)

!In[59]:= f[x_]:=Evaluate[InterpolatingPolynomial[{{-7,ym7},{-6,ym6},{-5,ym5},{-4,ym4},{-3,ym3},{-2,ym2},{-1,ym1},{0,y0},{1,yp1},{2,yp2},{3,yp3},{4,yp4},{5,yp5},{6,yp6},{7,yp7}},x]]
!In[60]:= ff=FortranForm[Expand[Simplify[ReplaceAll[D[f[x],x],{x->0}]]]]
!Out[60]//FortranForm=
!        (-7*ym1)/8. + (7*ym2)/24. - (7*ym3)/72. + (7*ym4)/264. - (7*ym5)/1320. + (7*ym6)/10296. - 
!     -  ym7/24024. + (7*yp1)/8. - (7*yp2)/24. + (7*yp3)/72. - (7*yp4)/264. + (7*yp5)/1320. - 
!     -  (7*yp6)/10296. + yp7/24024.
!In[61]:= ff=FortranForm[Expand[Simplify[ReplaceAll[D[f[x],x],{x->-1}]]]]
!Out[61]//FortranForm=
!        (8*y0)/7. - (15*ym1)/56. - (2*ym2)/3. + ym3/6. - (4*ym4)/99. + ym5/132. - (2*ym6)/2145. + 
!     -  ym7/18018. - yp1/2. + (2*yp2)/9. - yp3/12. + (4*yp4)/165. - yp5/198. + (2*yp6)/3003. - 
!     -  yp7/24024.
!In[62]:= ff=FortranForm[Expand[Simplify[ReplaceAll[D[f[x],x],{x->1}]]]]
!Out[62]//FortranForm=
!        (-8*y0)/7. + ym1/2. - (2*ym2)/9. + ym3/12. - (4*ym4)/165. + ym5/198. - (2*ym6)/3003. + 
!     -  ym7/24024. + (15*yp1)/56. + (2*yp2)/3. - yp3/6. + (4*yp4)/99. - yp5/132. + (2*yp6)/2145. - 
!     -  yp7/18018.
!In[63]:= ff=FortranForm[Expand[Simplify[ReplaceAll[D[f[x],x],{x->-2}]]]]
!Out[63]//FortranForm=
!        (-6*y0)/7. + (3*ym1)/2. - (275*ym2)/504. - ym3/2. + ym4/11. - ym5/66. + ym6/572. - ym7/10010. + 
!     -  yp1/2. - yp2/4. + yp3/10. - yp4/33. + yp5/154. - yp6/1144. + yp7/18018.
!In[64]:= ff=FortranForm[Expand[Simplify[ReplaceAll[D[f[x],x],{x->+2}]]]]
!Out[64]//FortranForm=
!       (6*y0)/7. - ym1/2. + ym2/4. - ym3/10. + ym4/33. - ym5/154. + ym6/1144. - ym7/18018. - 
!     -  (3*yp1)/2. + (275*yp2)/504. + yp3/2. - yp4/11. + yp5/66. - yp6/572. + yp7/10010.
!In[65]:= ff=FortranForm[Expand[Simplify[ReplaceAll[D[f[x],x],{x->-3}]]]]
!Out[65]//FortranForm=
!        (8*y0)/7. - (3*ym1)/2. + 2*ym2 - (2131*ym3)/2520. - (4*ym4)/11. + ym5/22. - (2*ym6)/429. + 
!     -  ym7/4004. - (3*yp1)/4. + (2*yp2)/5. - yp3/6. + (4*yp4)/77. - yp5/88. + (2*yp6)/1287. - 
!     -  yp7/10010.
!In[66]:= ff=FortranForm[Expand[Simplify[ReplaceAll[D[f[x],x],{x->+3}]]]]
!Out[66]//FortranForm=
!        (-8*y0)/7. + (3*ym1)/4. - (2*ym2)/5. + ym3/6. - (4*ym4)/77. + ym5/88. - (2*ym6)/1287. + 
!     -  ym7/10010. + (3*yp1)/2. - 2*yp2 + (2131*yp3)/2520. + (4*yp4)/11. - yp5/22. + (2*yp6)/429. - 
!     -  yp7/4004.
!In[67]:= ff=FortranForm[Expand[Simplify[ReplaceAll[D[f[x],x],{x->-4}]]]]
!Out[67]//FortranForm=
!        (-33*y0)/14. + (11*ym1)/4. - (11*ym2)/4. + (11*ym3)/4. - (32891*ym4)/27720. - ym5/4. + 
!     -  ym6/52. - ym7/1092. + (33*yp1)/20. - (11*yp2)/12. + (11*yp3)/28. - yp4/8. + yp5/36. - 
!    -  yp6/260. + yp7/4004.
!In[68]:= ff=FortranForm[Expand[Simplify[ReplaceAll[D[f[x],x],{x->+4}]]]]
!Out[68]//FortranForm=
!        (33*y0)/14. - (33*ym1)/20. + (11*ym2)/12. - (11*ym3)/28. + ym4/8. - ym5/36. + ym6/260. - 
!     -  ym7/4004. - (11*yp1)/4. + (11*yp2)/4. - (11*yp3)/4. + (32891*yp4)/27720. + yp5/4. - yp6/52. + 
!     -  yp7/1092.
!In[69]:= ff=FortranForm[Expand[Simplify[ReplaceAll[D[f[x],x],{x->-5}]]]]
!Out[69]//FortranForm=
!        (264*y0)/35. - (33*ym1)/4. + (22*ym2)/3. - (11*ym3)/2. + 4*ym4 - (44441*ym5)/27720. - 
!     -  (2*ym6)/13. + ym7/182. - (11*yp1)/2. + (22*yp2)/7. - (11*yp3)/8. + (4*yp4)/9. - yp5/10. + 
!     -  (2*yp6)/143. - yp7/1092.
!In[70]:= ff=FortranForm[Expand[Simplify[ReplaceAll[D[f[x],x],{x->+5}]]]]
!Out[70]//FortranForm=
!        (-264*y0)/35. + (11*ym1)/2. - (22*ym2)/7. + (11*ym3)/8. - (4*ym4)/9. + ym5/10. - (2*ym6)/143. + 
!     -  ym7/1092. + (33*yp1)/4. - (22*yp2)/3. + (11*yp3)/2. - 4*yp4 + (44441*yp5)/27720. + 
!     -  (2*yp6)/13. - yp7/182.
!In[71]:= ff=FortranForm[Expand[Simplify[ReplaceAll[D[f[x],x],{x->-6}]]]]
!Out[71]//FortranForm=
!        (-286*y0)/7. + (429*ym1)/10. - (143*ym2)/4. + (143*ym3)/6. - 13*ym4 + (13*ym5)/2. - 
!     -  (785633*ym6)/360360. - ym7/14. + (429*yp1)/14. - (143*yp2)/8. + (143*yp3)/18. - (13*yp4)/5. + 
!     -  (13*yp5)/22. - yp6/12. + yp7/182.
!In[72]:= ff=FortranForm[Expand[Simplify[ReplaceAll[D[f[x],x],{x->+6}]]]]
!Out[72]//FortranForm=
!        (286*y0)/7. - (429*ym1)/14. + (143*ym2)/8. - (143*ym3)/18. + (13*ym4)/5. - (13*ym5)/22. + 
!     -  ym6/12. - ym7/182. - (429*yp1)/10. + (143*yp2)/4. - (143*yp3)/6. + 13*yp4 - (13*yp5)/2. + 
!     -  (785633*yp6)/360360. + yp7/14.
!In[73]:= ff=FortranForm[Expand[Simplify[ReplaceAll[D[f[x],x],{x->-7}]]]]
!Out[73]//FortranForm=
!        (3432*y0)/7. - (1001*ym1)/2. + (2002*ym2)/5. - (1001*ym3)/4. + (364*ym4)/3. - (91*ym5)/2. + 
!     -  14*ym6 - (1171733*ym7)/360360. - (3003*yp1)/8. + (2002*yp2)/9. - (1001*yp3)/10. + 
!     -  (364*yp4)/11. - (91*yp5)/12. + (14*yp6)/13. - yp7/14.
!In[74]:= ff=FortranForm[Expand[Simplify[ReplaceAll[D[f[x],x],{x->+7}]]]]
!Out[74]//FortranForm=
!        (-3432*y0)/7. + (3003*ym1)/8. - (2002*ym2)/9. + (1001*ym3)/10. - (364*ym4)/11. + (91*ym5)/12. - 
!     -  (14*ym6)/13. + ym7/14. + (1001*yp1)/2. - (2002*yp2)/5. + (1001*yp3)/4. - (364*yp4)/3. + 
!     -  (91*yp5)/2. - 14*yp6 + (1171733*yp7)/360360.

   case(16)

!In[75]:= f[x_]:=Evaluate[InterpolatingPolynomial[{{-8,ym8},{-7,ym7},{-6,ym6},{-5,ym5},{-4,ym4},{-3,ym3},{-2,ym2},{-1,ym1},{0,y0},{1,yp1},{2,yp2},{3,yp3},{4,yp4},{5,yp5},{6,yp6},{7,yp7},{8,yp8}},x]]
!In[76]:= ff=FortranForm[Expand[Simplify[ReplaceAll[D[f[x],x],{x->0}]]]]
     c1D(-m:m,0) = (/&
     1.d0/102960.d0,-8.d0/45045.d0,2.d0/1287.d0,-56.d0/6435.d0,7.d0/198.d0,-56.d0/495.d0,14.d0/45.d0,&
     -8.d0/9.d0,0.d0,8.d0/9.d0,-14.d0/45.d0,56.d0/495.d0,-7.d0/198.d0,&
     56.d0/6435.d0,-2.d0/1287.d0,8.d0/45045.d0,-1.d0/102960.d0&
     /)
!In[77]:= ff=FortranForm[Expand[Simplify[ReplaceAll[D[f[x],x],{x->-1}]]]]
     c1D(-m:m,-1) = (/&
     -1.d0/80080.d0,1.d0/4290.d0,-3.d0/1430.d0,7.d0/572.d0,-7.d0/132.d0,21.d0/110.d0,-7.d0/10.d0,-17.d0/72.d0,&
     9.d0/8.d0,-1.d0/2.d0,7.d0/30.d0,-21.d0/220.d0,7.d0/220.d0,&
     -7.d0/858.d0,3.d0/2002.d0,-1.d0/5720.d0,1.d0/102960.d0&
     /)
!In[78]:= ff=FortranForm[Expand[Simplify[ReplaceAll[D[f[x],x],{x->1}]]]]
     c1D(-m:m,1) = (/&
     -1.d0/102960.d0,1.d0/5720.d0,-3.d0/2002.d0,7.d0/858.d0,-7.d0/220.d0,21.d0/220.d0,-7.d0/30.d0,1.d0/2.d0,&
     -9.d0/8.d0,17.d0/72.d0,7.d0/10.d0,-21.d0/110.d0,&
     7.d0/132.d0,-7.d0/572.d0,3.d0/1430.d0,-1.d0/4290.d0,1.d0/80080.d0&
     /)
!In[79]:= ff=FortranForm[Expand[Simplify[ReplaceAll[D[f[x],x],{x->-2}]]]]
!Out[79]//FortranForm=
     c1D(-m:m,-2) = (/&
      1.d0/48048.d0,-2.d0/5005.d0,15.d0/4004.d0,-10.d0/429.d0,5.d0/44.d0,-6.d0/11.d0,-1207.d0/2520.d0,10.d0/7.d0,&
      -45.d0/56.d0,10.d0/21.d0,-1.d0/4.d0,6.d0/55.d0,-5.d0/132.d0,10.d0/1001.d0,-15.d0/8008.d0,2.d0/9009.d0,-1.d0/80080.d0&
      /)
!In[80]:= ff=FortranForm[Expand[Simplify[ReplaceAll[D[f[x],x],{x->+2}]]]]
!Out[80]//FortranForm=
     c1D(-m:m,2) = (/&
     1.d0/80080.d0,-2.d0/9009.d0,15.d0/8008.d0,-10.d0/1001.d0,5.d0/132.d0,-6.d0/55.d0,1.d0/4.d0,-10.d0/21.d0,45.d0/56.d0,&      
     -10.d0/7.d0,1207.d0/2520.d0,6.d0/11.d0,-5.d0/44.d0,10.d0/429.d0,-15.d0/4004.d0,2.d0/5005.d0,-1.d0/48048.d0&
     /)
!In[81]:= ff=FortranForm[Expand[Simplify[ReplaceAll[D[f[x],x],{x->-3}]]]]
!Out[81]//FortranForm=
     c1D(-m:m,-3) = (/&
     -1.d0/21840.d0,1.d0/1092.d0,-5.d0/546.d0,5.d0/78.d0,-5.d0/12.d0,-20417.d0/27720.d0,11.d0/6.d0,-55.d0/42.d0,55.d0/56.d0,&
     -55.d0/84.d0,11.d0/30.d0,-1.d0/6.d0,5.d0/84.d0,-5.d0/312.d0,5.d0/1638.d0,-1.d0/2730.d0,1.d0/48048.d0&
     /)
!In[82]:= ff=FortranForm[Expand[Simplify[ReplaceAll[D[f[x],x],{x->+3}]]]]
!Out[82]//FortranForm=
     c1D(-m:m,3) = (/&
     -1.d0/48048.d0,1.d0/2730.d0,-5.d0/1638.d0,5.d0/312.d0,-5.d0/84.d0,1.d0/6.d0,-11.d0/30.d0,55.d0/84.d0,&
     -55.d0/56.d0,55.d0/42.d0,-11.d0/6.d0,20417.d0/27720.d0,5.d0/12.d0,-5.d0/78.d0,5.d0/546.d0,-1.d0/1092.d0,1.d0/21840.d0&
     /)
!In[83]:= ff=FortranForm[Expand[Simplify[ReplaceAll[D[f[x],x],{x->-4}]]]]
!Out[83]//FortranForm=
     c1D(-m:m,-4) = (/&
     1.d0/7280.d0,-4.d0/1365.d0,3.d0/91.d0,-4.d0/13.d0,-28271.d0/27720.d0,12.d0/5.d0,-11.d0/5.d0,+44.d0/21.d0,&
     -99.d0/56.d0,44.d0/35.d0,-11.d0/15.d0,12.d0/35.d0,-1.d0/8.d0,&
     4.d0/117.d0,-3.d0/455.d0,4.d0/5005.d0,-1.d0/21840.d0&
     /)
!In[84]:= ff=FortranForm[Expand[Simplify[ReplaceAll[D[f[x],x],{x->+4}]]]]
!Out[84]//FortranForm=
     c1D(-m:m,4) = (/&
     1.d0/21840.d0,-4.d0/5005.d0,3.d0/455.d0,-4.d0/117.d0,1.d0/8.d0,-12.d0/35.d0,11.d0/15.d0,-44.d0/35.d0,&
     99.d0/56.d0,-44.d0/21.d0,11.d0/5.d0,-12.d0/5.d0,&
     28271.d0/27720.d0,4.d0/13.d0,-3.d0/91.d0,4.d0/1365.d0,-1.d0/7280.d0&
     /)
!In[85]:= ff=FortranForm[Expand[Simplify[ReplaceAll[D[f[x],x],{x->-5}]]]]
!Out[85]//FortranForm=
     c1D(-m:m,-5) = (/&
     -1.d0/1680.d0,1.d0/70.d0,-3.d0/14.d0,-485333.d0/360360.d0,13.d0/4.d0,-39.d0/10.d0,143.d0/30.d0,-143.d0/28.d0,&
     1287.d0/280.d0,-143.d0/42.d0,143.d0/70.d0,-39.d0/40.d0,13.d0/36.d0,-1.d0/10.d0,3.d0/154.d0,-1.d0/420.d0,1.d0/7280.d0&
     /)
!In[86]:= ff=FortranForm[Expand[Simplify[ReplaceAll[D[f[x],x],{x->+5}]]]]
!Out[86]//FortranForm=
    c1D(-m:m,5) = (/&
    -1.d0/7280.d0,1.d0/420.d0,-3.d0/154.d0,1.d0/10.d0,-13.d0/36.d0,39.d0/40.d0,-143.d0/70.d0,143.d0/42.d0,&
    -1287.d0/280.d0,143.d0/28.d0,-143.d0/30.d0,39.d0/10.d0,-13.d0/4.d0,485333.d0/360360.d0,3.d0/14.d0,-1.d0/70.d0,1.d0/1680.d0&
    /)
!In[87]:= ff=FortranForm[Expand[Simplify[ReplaceAll[D[f[x],x],{x->-6}]]]]
!Out[87]//FortranForm=
    c1D(-m:m,-6) = (/&
    1.d0/240.d0,-2.d0/15.d0,-631193.d0/360360.d0,14.d0/3.d0,-91.d0/12.d0,182.d0/15.d0,-1001.d0/60.d0,286.d0/15.d0,&
    -143.d0/8.d0,286.d0/21.d0,-1001.d0/120.d0,182.d0/45.d0,-91.d0/60.d0,14.d0/33.d0,-1.d0/12.d0,2.d0/195.d0,-1.d0/1680.d0&
    /)
!In[88]:= ff=FortranForm[Expand[Simplify[ReplaceAll[D[f[x],x],{x->+6}]]]]
!Out[88]//FortranForm=
    c1D(-m:m,6) = (/&
    1.d0/1680.d0,-2.d0/195.d0,1.d0/12.d0,-14.d0/33.d0,91.d0/60.d0,-182.d0/45.d0,1001.d0/120.d0,-286.d0/21.d0,&
    143.d0/8.d0,-286.d0/15.d0,1001.d0/60.d0,-182.d0/15.d0,91.d0/12.d0,-14.d0/3.d0,631193.d0/360360.d0,2.d0/15.d0,-1.d0/240.d0&
    /)
!In[89]:= ff=FortranForm[Expand[Simplify[ReplaceAll[D[f[x],x],{x->-7}]]]]
!Out[89]//FortranForm=
    c1D(-m:m,-7) = (/&
    -1.d0/16.d0,-835397.d0/360360.d0,15.d0/2.d0,-35.d0/2.d0,455.d0/12.d0,-273.d0/4.d0,1001.d0/10.d0,-715.d0/6.d0,&
    6435.d0/56.d0,-715.d0/8.d0,1001.d0/18.d0,-273.d0/10.d0,455.d0/44.d0,-35.d0/12.d0,15.d0/26.d0,-1.d0/14.d0,1.d0/240.d0&
    /)
!In[90]:= ff=FortranForm[Expand[Simplify[ReplaceAll[D[f[x],x],{x->+7}]]]]
!Out[90]//FortranForm=
    c1D(-m:m,7) = (/&
    -1.d0/240.d0,1.d0/14.d0,-15.d0/26.d0,35.d0/12.d0,-455.d0/44.d0,273.d0/10.d0,-1001.d0/18.d0,715.d0/8.d0,&
    -6435.d0/56.d0,715.d0/6.d0,-1001.d0/10.d0,273.d0/4.d0,-455.d0/12.d0,35.d0/2.d0,-15.d0/2.d0,835397.d0/360360.d0,1.d0/16.d0&
    /)
!In[91]:= ff=FortranForm[Expand[Simplify[ReplaceAll[D[f[x],x],{x->-8}]]]]
!Out[91]//FortranForm=
    c1D(-m:m,-8) = (/&
    -2436559.d0/720720.d0,16.d0,-60.d0,560.d0/3.d0,-455.d0,4368.d0/5.d0,-4004.d0/3.d0,11440.d0/7.d0,&
    -6435.d0/4.d0,11440.d0/9.d0,-4004.d0/5.d0,4368.d0/11.d0,-455.d0/3.d0,560.d0/13.d0,-60.d0/7.d0,16.d0/15.d0,-1.d0/16.d0&
    /)
!In[92]:= ff=FortranForm[Expand[Simplify[ReplaceAll[D[f[x],x],{x->+8}]]]]
!Out[92]//FortranForm=
    c1D(-m:m,8) = (/&
    1.d0/16.d0,-16.d0/15.d0,60.d0/7.d0,-560.d0/13.d0,455.d0/3.d0,-4368.d0/11.d0,4004.d0/5.d0,-11440.d0/9.d0,&
    6435.d0/4.d0,-11440.d0/7.d0,4004.d0/3.d0,-4368.d0/5.d0,455.d0,-560.d0/3.d0,60.d0,-16.d0,2436559.d0/720720.d0&
    /)

   end select
