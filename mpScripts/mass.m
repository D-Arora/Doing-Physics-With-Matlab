%   *******************************************************************
function m = mass(Z,A)
% NUCLEAR MASSES (not atom) [amu u] 
%  Element: atomic number Z  / Isotope: mass number A
%   New values for the mass of an isotope can be easily added to the data
%   base. Mass values can be found at
%         https://wwwndc.jaea.go.jp/NuC/


% mass electron [amu]
me = 5.485799093287202e-04;
  
if Z == 1     % hydrogen  
   if A == 1; m = 1.00782503223; end
   if A == 2; m = 2.01410177812; end
   if A == 3; m = 3.01604927790; end
end   

if Z == 2     % heliumn  
   if A == 3; m = 3.01602932007; end
   if A == 4; m = 4.00260325413; end
end
   
if Z == 3     % lithium 
   if A == 6;  m = 6.01512288741;   end 
   if A == 7;  m = 7.01600342665;   end
   if A == 8;  m = 8.022486238;     end
   if A == 9;  m = 9.026790199;     end
   if A == 11; m = 11.043723754;    end
end   

if Z == 4     % beryllium 
   if A == 7;  m = 7.016928707;   end 
   if A == 8;  m = 8.005305102;   end
   if A == 9;  m = 9.012183050;   end
   if A == 10; m = 10.013534679;  end
   if A == 11; m = 11.021661078;  end
   if A == 12; m = 12.026920848 ; end
   if A == 14; m = 14.042892920;  end
end

if Z == 6     % carbon
   if A == 11; m = 11.011433611;   end
   if A == 12; m = 12.0000000000;  end
   if A == 13; m = 13.00335483507; end
   if A == 14; m = 14.00324198842; end
end


if Z == 7     % nitrogenn  
   if A == 13; m = 13.005738609; end
   if A == 14; m = 14.00307400443; end
   if A == 15; m = 15.00010889888; end
end
  
if Z == 8     % oxygen  
   if A == 15; m = 15.003065618; end  
   if A == 16; m = 15.99491461957; end  
   if A == 17; m = 16.99913175650; end  
end   
   
if Z == 9     % fluorine 
   if A == 18; m = 18.000937325; end  
   if A == 19; m = 18.99840316273; end  
   if A == 20; m = 19.999981252 ; end  
end      
   
if Z == 10     % neon  
   if A == 18; m = 18.005708704; end  
   if A == 19; m = 19.001880907; end  
   if A == 20; m = 19.99244017617; end   
end 

if Z == 25     % iron  
   if A == 55; m = 54.938043937; end  
end 

if Z == 26     % iron  
   if A == 55; m = 54.938292023; end  
   if A == 56; m = 55.934936355; end 
end 

if Z == 36     % krypton  
   if A == 83; m = 82.914127143;   end 
   if A == 84; m = 83.91149772821; end
   if A == 85; m = 84.912527262;   end
   if A == 86; m = 85.910610627;   end
   if A == 87; m = 86.913354760;   end
   if A == 88; m = 87.914447881;   end
   if A == 89; m = 88.917835451;   end
   if A == 90; m = 90.923806311;   end
   if A == 91; m = 90.923806311;   end
   if A == 92; m = 91.926173095;   end
   if A == 93; m = 92.931147175;   end
   if A == 94; m = 93.934140455;   end
   if A == 95; m = 94.939710924;   end
   if A == 96; m = 95.943016618;   end
   if A == 97; m = 96.949088785;   end
end 

if Z == 37     % rubidium  
   if A == 90; m = 89.914798243; end 
end

if Z == 38     % strontium
   if A == 88; m = 87.905612509; end 
   if A == 89; m = 88.907451063; end                
end


if Z == 40     % zirconium 
   if A == 99;  m = 98.916664513; end 
   if A == 100; m =  99.917997918  ; end                
end

if Z == 52     % tellurium 
   if A == 134; m = 133.911368729; end 
   if A == 135; m = 134.916367067; end 
end

if Z == 54     % Xenon
   if A == 136; m = 135.907214484; end 
   if A == 137; m = 136.911557815; end 
end

if Z == 55     % cesium  
   if A == 143; m = 142.927347969; end 
end

if Z == 56     % barium  
   if A == 138; m = 137.905246899; end 
   if A == 139 ; m = 138.908841004; end
   if A == 140; m = 139.910605633; end
   if A == 141; m = 140.914402957; end
   if A == 142; m = 141.916429935; end
   if A == 143; m = 142.920625195; end
   if A == 144; m = 143.922954824; end
   if A == 145; m = 144.927518400; end
   if A == 146; m = 145.930283286; end
   if A == 147; m = 146.935303900; end
   if A == 148; m = 147.938170578; end
   if A == 149; m = 148.942920; end
   if A == 150; m = 149.945950; end
   if A == 151; m = 150.951080; end
end 

if Z == 86     % radon 
   if A == 222; m = 222.017577269; end  
   if A == 226; m = 226.030852862; end  
end 

if Z == 88     % radium  
   if A == 226; m = 226.025409353; end  
   if A == 228; m = 228.031074994; end 
end

if Z == 90     % thorium  
   if A == 226; m = 226.024891; end   
   if A == 228; m = 228.028734823; end
   if A == 232; m = 232.038060026; end
   if A == 234; m = 234.043602450; end
end

if Z == 91     % protactinium 
   if A == 231; m = 231.035884277; end 
   if A == 237; m = 237.051147110; end  
end

if Z == 92     % uranium  
   if A == 217; m = 217.024382652; end  
   if A == 218; m = 218.023522502; end
   if A == 219; m = 219.025016432; end
   if A == 222; m = 222.026081   ; end
   if A == 223; m = 223.027737909; end
   if A == 224; m = 224.027603882; end
   if A == 225; m = 225.029389831; end
   if A == 226; m = 226.029337817; end
   if A == 227; m = 227.031155482; end
   if A == 228; m = 228.031373109; end
   if A == 229; m = 229.033505084; end
   if A == 230; m = 230.033938943; end
   if A == 231; m = 231.036293977; end
   if A == 232; m = 232.037149836; end
   if A == 233; m = 233.039636574; end
   if A == 234; m = 234.040953616; end
   if A == 235; m = 235.043931368; end
   if A == 236; m = 236.045569468; end
   if A == 237; m = 237.048731636; end
   if A == 238; m = 238.050789466; end
   if A == 239; m = 239.054294518; end
   if A == 240; m = 240.056593384; end
   if A == 242; m = 242.062933   ; end
end 
   
  m = m - Z*me;
end


