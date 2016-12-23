$Title: Biodiesel TES Model - 12/23/2016
$ontext
Modified: 12/23/2016
This contains Technological and Ecological systems and choice is between the best option in Summer to meet air quality demand
$offtext

$OFFSYMXREF
$OFFSYMLIST

Options LimRow = 70, LimCol = 0;


Set
*Global settings
*        Define units: Src=source, Snk=sink, Mix=mixer, Spl=splitter
*                      Str=storage, HX = Heat exchangers

         unit    units
         /Reactor3, Reactor4, Sep2*Sep3, Mix2*Mix3, Col3*Col5,
          Snk3*Snk7, Src3, Src5, Src6, Spl3, Src, HX10*HX23, HX45, Src7, CHP, SrcWa, Wetland/

*        Define subsets of units
         HX(unit) heat exchangers
        /HX10*HX23, HX45 /

         Mix(unit) mixers
         /Mix2*Mix3/

         Spl(unit) splitter
         /Spl3/

         Src(unit) sources
         /Src3, Src5, Src6, Src, Src7/

         Snk(unit) sinks
         /Snk3*Snk7/

         Sep(unit) separators
         /Sep2*Sep3/

         Column(unit) distillation columns
         /Col3*Col5/


*        Define Components
         J       components
         /Wa, MeOH, Glycerol, KOH, Oil, FFA, FAME, PhosA, PotPhos/

*        Define Reactions: KOH_Neut - neutralization reaction of KOH with phosphoric acid
         react
         /KOH_neut/

*        running variable for vapor pressure correlation
         l /1*3/

*         Pollutant emissions
         K       Pollutants
         /NO2, SO2, PM10, CO2, O3, GWater/

*        running variable for flux correlation
         p  /1*3/

Alias(unit, unit1)
Alias (J,h)


Parameters
*        individual liquid heat capacity of a component (average in a range 20 C - 100 C)
*        in kJ/(kg*C), assume: constant heat capacities and c_p_ind('CellM')=c_p_ind('Wa')

         c_p_ind(J)
         /Wa             4.19
          MeOH           2.534
          Glycerol       2.409
          KOH            1.179
          Oil            2.000
          FFA            2.07
          FAME           2.01
          PhosA          1.758
          PotPhos        4.07/

*        Molecular weight in g/mol
         MW(J)
         /Wa             18.015
          MeOH           32.04
          Glycerol       92.093
          KOH            56.1056
          Oil            920
          FFA            280.5
          FAME           292.2
          PhosA          98
          PotPhos        174.18/


*        Standard heats of vapourization Kj/Kg
         dH_vap_0(J)
         /Wa             2254
          MeOH           1000
          Glycerol       974
          Oil            170.137
          FAME           460/;



Table
*        vapor pressure coefficients: mmHq, T= C
*        750 mmHq = 1bar, 760mmHq = 1atm
*        ln(p_k)j= coef_p(h,1)-coef_p(h,2)/(T+coef_p(h,3))
         coef_p(J,l)
                      1       2       3
         Wa         5.0768  1659.73  -45.854
         MeOH       5.20409 1581.341 -33.50
         Oil        11.4785 -708.72  -167.48
         FAME       8.2175  1450.62  88.03
         Glycerol   41.9806 6024.4   10.3272
         KOH        7.06223 6300.48  172.76
         FFA        3.82846 2066.995 -116.87;

*FAME in K, MeOH in K, Oil in K, Wa in K, KOH in C, FFA in K
*MeOH in bar, FAME in Pa, Oil in Pa, glycerol in Torr, Wa in bar, KOH in mmHg, FFA in bar;
*** NOTE ANTOINE EQUATION FOR OIL AND GLYCEROL ARE SLIGHTLY DIFF FROM THE ORIGINAL EQUATION **
*Vp for KOH is coming to zero now
* KOH - KOH        4.02783 5854.906 144.113

SET      Arc(unit,unit1) stream matrix;

*setting entries in stream matrix
***Arc is unknown symbol***
Arc(unit, unit1)=No;

*Define all existing streams Arc('1','2') is stream from unit 1 to unit 2


*Feed to reactors
Arc('Mix2','Mix3') = Yes;
Arc('Src','Mix3') = Yes;
Arc('Col3','Mix2') = Yes;
Arc('Src3','Mix2') = Yes;
Arc('Mix3','HX45') = Yes;
Arc('Src5','Mix3') = Yes;

*Reaction transesterification unit
Arc('HX45','Reactor3') = Yes;
Arc('Mix3','HX45') = Yes;
Arc('Reactor3','HX10') = Yes;

*MeOH distillation unit
Arc('HX10','Col3') = Yes;
Arc('Col3','HX13') = Yes;
Arc('Col3','Mix2') = Yes;

*Washing
Arc('HX13','Sep2') = Yes;
Arc('Src6','Sep2') = Yes;
Arc('Sep2','Reactor4') = Yes;
Arc('Sep2','HX19') = Yes;
Arc('SrcWa','Src6') = Yes;

*Neutralization alkali
Arc('Reactor4','HX14') = Yes;
Arc('Src7','Reactor4')  = Yes;
Arc('Reactor4','Snk3') = Yes;

*FAME Distillation
Arc('HX19','Col5') = Yes;
Arc('Col5','HX22') = Yes;
Arc('Col5','HX23') = Yes;

*Glycerol distillation
Arc('HX14','Col4') = Yes;
Arc('Col4','HX17') = Yes;
Arc('Col4','HX18') = Yes;

*Wetland
Arc('Wetland','Src6') = Yes;

*CHP
Arc('SrcWa','Src6') = Yes;
Arc('Src6','CHP') = Yes;



Positive Variables
*        streams and mass fractions: all in kg/s
         F(unit,unit1)           total streams in kg s^-1
         fc(J,unit,unit1)        individual components streams
         x(J,unit,unit1)         mass fraction of comp J in stream

*        heat of evaporation
         dH_v(J,unit,unit1)      individual heat of vap. (KJ per kg)

*        vapor pressure
         p_v(J,unit,unit1)       vapor pressure in bar

*        temperatures in C
         T(unit,unit1)          temperature of stream in C

*        Separation efficiencies
         Separation(J)                  Separation efficiencies

Variables
*        heat
         Q(Unit)         heat produced or consumed in unit (efficiency included)
*         Q_cond(column)  heat load of condenser of column
*         Q_reb(column)   heat load of reboiler of column

         Z               objective function value ;


*Defining global bounds and fixing some variables

*mass fractions
x.UP(J,unit,unit1)$Arc(unit,unit1)=1;

*Total streams. #Setting this at 1000 but look into Python code for changing value#
F.up(unit,unit1)$Arc(unit,unit1)=100;
F.lo(unit,unit1)$Arc(unit,unit1) = 0.001;

*Component streams #Setting this at 1000 but look into Python code for changing value#
fc.up(J,unit,unit1)$Arc(unit,unit1)=20;
fc.lo(J,unit,unit1)$Arc(unit,unit1)=0.001;

F.up(unit,unit1)$(not Arc(unit,unit1)) = 0;
fc.up(J,unit,unit1)$(not Arc(unit,unit1)) = 0;


*Heat consumption of certain units
Q.Fx('Sep2') = 0;
Q.Fx('Mix2') = 0;
Q.Fx('Mix3') = 0;
Q.Fx('Col3') = 0;
Q.Fx('Col4') = 0;
Q.Fx('Col5') = 0;
Q.Fx('Snk3') = 0;
Q.Fx('Snk4') = 0;
Q.Fx('Snk5') = 0;
Q.Fx('Snk6') = 0;
Q.Fx('Snk7') = 0;
Q.Fx('Src3') = 0;
Q.Fx('Src5') = 0;
Q.Fx('Src7') = 0;
Q.Fx('Spl3') = 0;
Q.Fx('Src') = 0;

*Global bounds for separation efficiencies of components
*MeOH
Separation.lo('MeOH')  = 0.5;
Separation.up('MeOH') = 0.9999;
Separation.l('MeOH') = (Separation.up('MeOH') + Separation.lo('MeOH'))/2;
*Water
Separation.lo('Wa')  = 0.5;
Separation.up('Wa') = 0.9999;
Separation.l('Wa') = (Separation.up('Wa') + Separation.lo('Wa'))/2;
*Separation of water is always 0.5, so giving it a higher value makes this equation infeasible.

*Glycerol
Separation.lo('Glycerol')  = 0.5;
Separation.up('Glycerol') = 0.9999;
Separation.l('Glycerol') = (Separation.up('Glycerol') + Separation.lo('Glycerol'))/2;
*Oil
Separation.lo('Oil')  = 0.01;
Separation.up('Oil') = 0.3;
Separation.l('Oil') = (Separation.up('Oil') + Separation.lo('Oil'))/2;
*FFA
Separation.lo('FFA')  = 0.01;
Separation.up('FFA') = 0.5;
Separation.l('FFA') = (Separation.up('FFA') + Separation.lo('FFA'))/2;
*KOH
Separation.lo('KOH')  = 0.5;
Separation.up('KOH') = 0.9999;
Separation.l('KOH') = (Separation.up('KOH') + Separation.lo('KOH'))/2;
*FAME
Separation.lo('FAME')  = 0.01;
Separation.up('FAME') = 0.5;
Separation.l('FAME') = (Separation.up('FAME') + Separation.lo('FAME'))/2;

Scalar
         dp                     pressure loss in the column /0.1/
         alpha_c4               relative volatility in column 4 /3.8/
         alpha_c3                relative volatility in column 3 /4.5/
         dhreactn                Neutralization heat Kj per mol /-200.18/
         dhreactneu              Neutralization heat in Kh per mol /-173/
         washing_wa              Water washing recovery        /0.01/
         Temp_cooldown           Cooldown temperature of FAME   /25/
         T_amb                   Ambient teperature              /25/;

*recover_oil             recovery of Oil that is fixed /0.997/
*recover_FAME            recovery of FAME that is fixed /0.997/

*global temperature bounds - bounds get redefined for specific streams
T.lo(unit,unit1)=0.01;
*Set the lowest temp to be a value greater than zero
T.up(unit,unit1)=300;
T.l(unit, unit1) = (T.lo(unit, unit1) + T.up(unit, unit1))/2;

*Specifying temperatures in C
T.Fx('Src3','Spl3') = 25;
T.Fx('Spl3','Mix2') = 25;
T.Fx('Src','Mix3') = 60;
T.Fx('Mix3','HX45') = 37.23;
T.Fx('HX11','Mix2') = 28.2;
T.Fx('Src5','Mix3') = 25;
T.Fx('HX45','Reactor3') = 60;
T.Fx('Reactor3','HX10') = 60;
T.Fx('HX10','Col3') = 80;
*T.Fx('Col3','HX13') = 145;
T.Fx('HX13','Sep2') = 60;
T.Fx('Src7','Reactor4') = 60;
T.Fx('Reactor4','Snk3') = 60;
*Temp value from HX14 to Col 4 value needs to be checked again
T.Fx('HX14','Col4') = 60;
T.Fx('HX17','Snk4') = 42.8;
T.Fx('HX16','HX18') = 112;
T.Fx('HX18','Snk5') = 148.6;
T.Fx('Sep2','Reactor4') = 60;
T.Fx('HX19','Col5') = 80;
T.Fx('HX22','Snk6') = 167.8;
T.Fx('Src6','Sep2') = 20;
*Inlet temperature of water



$ontext
This section  describes all src streams that have a pre-fixed value. In this model,
mass flow of methanol stream entering from src3 to mix2, oil from src to mix3, KOH from src5 to mix3
and Wa from src6 to Sep2 are decision variables
$offtext
*################################################*************
*fc.up('MeOH','Src3','Mix2') = 0.059137;
*fc.lo('MeOH','Src3','Mix2') = 0.038382;
*################################################*************
fc.up('MeOH','Src3','Mix2') = 10;
fc.lo('MeOH','Src3','Mix2') = 0.01;
fc.l('MeOH','Src3','Mix2') = (fc.lo('MeOH','Src3','Mix2') + fc.up('MeOH','Src3','Mix2'))/2;
*fc.L = inital point will be min max of this

fc.FX('Wa','Src3','Mix2') = 0;
fc.FX('Glycerol','Src3','Mix2') = 0;
fc.FX('KOH','Src3','Mix2') = 0;
fc.FX('Oil','Src3','Mix2') = 0;
fc.FX('FFA','Src3','Mix2') = 0;
fc.FX('FAME','Src3','Mix2') = 0;
fc.FX('PhosA','Src3','Mix2') = 0;
fc.FX('PotPhos','Src3','Mix2') = 0;

x.FX('MeOH','Src3','Mix2') = 1;

x.FX('Wa','Src3','Mix2') = 0;
x.FX('Glycerol','Src3','Mix2') = 0;
x.FX('KOH','Src3','Mix2') = 0;
x.FX('Oil','Src3','Mix2') = 0;
x.FX('FFA','Src3','Mix2') = 0;
x.FX('FAME','Src3','Mix2') = 0;
x.FX('PhosA','Src3','Mix2') = 0;
x.FX('PotPhos','Src3','Mix2') = 0;

****************


fc.FX('Wa','Src5','Mix3') = 0;
fc.FX('MeOH','Src5','Mix3') = 0;
fc.FX('Glycerol','Src5','Mix3') = 0;
fc.FX('Oil','Src5','Mix3') = 0;
fc.FX('FFA','Src5','Mix3') = 0;
fc.FX('FAME','Src5','Mix3') = 0;
fc.FX('PhosA','Src5','Mix3') = 0;
fc.FX('PotPhos','Src5','Mix3') = 0;

x.FX('KOH','Src5','Mix3') = 1;

x.FX('Wa','Src5','Mix3') = 0;
x.FX('MeOH','Src5','Mix3') = 0;
x.FX('Glycerol','Src5','Mix3') = 0;
x.FX('Oil','Src5','Mix3') = 0;
x.FX('FFA','Src5','Mix3') = 0;
x.FX('FAME','Src5','Mix3') = 0;
x.FX('PhosA','Src5','Mix3') = 0;
x.FX('PotPhos','Src5','Mix3') = 0;



*Src has only Oil as the input which is a variable
*################################################*************
*fc.up('Oil','Src','Mix3') = 0.59222947;
*fc.lo('Oil','Src','Mix3') = 0.3947620;
*################################################*************
fc.up('Oil','Src','Mix3') = 5;
fc.lo('Oil','Src','Mix3') = 0.01;

fc.l('Oil','Src','Mix3') = (fc.lo('Oil','Src','Mix3') + fc.up('Oil','Src','Mix3'))/2;
*fc.L


fc.FX('Wa','Src','Mix3') = 0;
fc.FX('MeOH','Src','Mix3') = 0;
fc.FX('Glycerol','Src','Mix3') = 0;
fc.FX('KOH','Src','Mix3') = 0;
fc.FX('FFA','Src','Mix3') = 0;
fc.FX('FAME','Src','Mix3') = 0;
fc.FX('PhosA','Src','Mix3') = 0;
fc.FX('PotPhos','Src','Mix3') = 0;

x.FX('Oil','Src','Mix3') = 1;

x.FX('Wa','Src','Mix3') = 0;
x.FX('MeOH','Src','Mix3') = 0;
x.FX('Glycerol','Src','Mix3') = 0;
x.FX('KOH','Src','Mix3') = 0;
x.FX('FFA','Src','Mix3') = 0;
x.FX('FAME','Src','Mix3') = 0;
x.FX('PhosA','Src','Mix3') = 0;
x.FX('PotPhos','Src','Mix3') = 0;


************* Water to washing stream ******

fc.FX('MeOH','Src6','Sep2') = 0;
fc.FX('Glycerol','Src6','Sep2') = 0;
fc.FX('KOH','Src6','Sep2') = 0;
fc.FX('Oil','Src6','Sep2') = 0;
fc.FX('FFA','Src6','Sep2') = 0;
fc.FX('FAME','Src6','Sep2') = 0;
fc.FX('PhosA','Src6','Sep2') = 0;
fc.FX('PotPhos','Src6','Sep2') = 0;

x.FX('Wa','Src6','Sep2') = 1;

x.FX('MeOH','Src6','Sep2') = 0;
x.FX('Glycerol','Src6','Sep2') = 0;
x.FX('KOH','Src6','Sep2') = 0;
x.FX('Oil','Src6','Sep2') = 0;
x.FX('FFA','Src6','Sep2') = 0;
x.FX('FAME','Src6','Sep2') = 0;
x.FX('PhosA','Src6','Sep2') = 0;
x.FX('PotPhos','Src6','Sep2') = 0;



*Src7 has PhosA as input which can be
*fc.up('PhosA','Src7','Reactor4') = 200;
*fc.lo('PhosA','Src7','Reactor4') = 0.001;
*fc.l('PhosA','Src7','Reactor4') = (fc.lo('PhosA','Src7','Reactor4') + fc.up('PhosA','Src7','Reactor4'))/2;

fc.Fx('Wa','Src7','Reactor4') = 0;
fc.Fx('MeOH','Src7','Reactor4') = 0;
fc.Fx('Glycerol','Src7','Reactor4') = 0;
fc.Fx('KOH','Src7','Reactor4') = 0;
fc.Fx('Oil','Src7','Reactor4') = 0;
fc.Fx('FFA','Src7','Reactor4') = 0;
fc.Fx('FAME','Src7','Reactor4') = 0;
fc.Fx('PotPhos','Src7','Reactor4') = 0;

x.Fx('PhosA','Src7','Reactor4') = 1;

x.Fx('Wa','Src7','Reactor4') = 0;
x.Fx('MeOH','Src7','Reactor4') = 0;
x.Fx('Glycerol','Src7','Reactor4') = 0;
x.Fx('KOH','Src7','Reactor4') = 0;
x.Fx('Oil','Src7','Reactor4') = 0;
x.Fx('FFA','Src7','Reactor4') = 0;
x.Fx('FAME','Src7','Reactor4') = 0;
x.Fx('PotPhos','Src7','Reactor4') = 0;


**------- Fixing flow streams for Phosphoric Acid ---- **
fc.Fx('PhosA','Reactor4','HX14') =0;
fc.Fx('PhosA','Sep2','Reactor4') =0;
fc.Fx('PhosA','Mix2','Mix3') =0;
fc.Fx('PhosA','Mix3','HX45') =0;
fc.Fx('PhosA','HX10','Col3') =0;
fc.Fx('PhosA','HX45','Reactor3') =0;
x.Fx('PhosA','Reactor4','HX14') =0;
x.Fx('PhosA','Sep2','Reactor4') =0;
x.Fx('PhosA','Mix2','Mix3') =0;
x.Fx('PhosA','Mix3','HX45') =0;
x.Fx('PhosA','HX10','Col3') =0;
x.Fx('PhosA','HX45','Reactor3') =0;



**------- Fixing flow streams for Potassium Phosphate ---- **
fc.Fx('PotPhos','Mix2','Mix3') = 0;
fc.Fx('PotPhos','Mix3','HX45') = 0;
fc.Fx('PotPhos','HX10','Col3') = 0;
fc.Fx('PotPhos','HX45','Reactor3') = 0;
x.Fx('PotPhos','Mix2','Mix3') = 0;
x.Fx('PotPhos','Mix3','HX45') = 0;
x.Fx('PotPhos','HX10','Col3') = 0;
x.Fx('PotPhos','HX45','Reactor3') = 0;


** ------ Fixing flow of FFA ---- **
fc.Fx('FFA','Col5','HX23') = 0;
x.FX('FFA','Col5','HX23') = 0;

** ----------- Fixing glycerol flow stream ---- **
fc.Fx('Glycerol','Mix2','Mix3') = 0;
x.Fx('Glycerol','Mix2','Mix3') = 0;
fc.Fx('Glycerol','Mix3','HX45') = 0;
x.Fx('Glycerol','Mix3','HX45') = 0;
fc.Fx('Glycerol','HX45','Reactor3') = 0;
x.Fx('Glycerol','HX45','Reactor3') = 0;
*fc.Fx('Glycerol','Col4','HX17') = 0;
*x.Fx('Glycerol','Col4','HX17') = 0;

** ------ Fixing FAME flow streams ---**
fc.Fx('FAME','Mix2','Mix3') = 0;
x.Fx('FAME','Mix2','Mix3') = 0;
fc.Fx('FAME','Mix3','HX45') = 0;
x.Fx('FAME','Mix3','HX45') = 0;
fc.Fx('FAME','HX45','Reactor3')= 0;
x.Fx('FAME','HX45','Reactor3') = 0;





Positive variables
         Tmp_reactor          Temperature of reactor 4
         cat                  Catalyst ratio in input
         ratio_met_alk        Ratio of methanol to alkali
         Op                   Output
         conver               conversion rate;



Tmp_reactor.lo = 45;
Tmp_reactor.up = 65;
Tmp_reactor.l = (Tmp_reactor.lo + Tmp_reactor.up)/2;
cat.lo = 0.005;
cat.up = 0.015;
cat.l = (cat.lo + cat.up)/2;
ratio_met_alk.lo = 4.5;
ratio_met_alk.up = 7.5;
ratio_met_alk.l = (ratio_met_alk.lo + ratio_met_alk.up)/2;

*Src5 has only KOH as the input which is a variable
*################################################*************
*fc.up('KOH','Src5','Mix3') = 0.005827;
*fc.lo('KOH','Src5','Mix3') = 0.001847;
*################################################*************
fc.up('KOH','Src5','Mix3') = 5;
fc.lo('KOH','Src5','Mix3') = 0.01;

fc.l('KOH','Src5','Mix3') = (fc.up('KOH','Src5','Mix3') + fc.lo('KOH','Src5','Mix3'))/2;
*fc.L will be 0.5(min + max)



** Defininte scalars for Cost of different materials **
Scalar

*         Cost_FAME       Cost of biodiesel /1.4/
*         Cost_Glycerol   Cost of glycerol /1.026/               (Old value is 0.6)
*         Cost_PotPhos    Cost of potassium phosphate /3.491/       (Old value is 1.76)
*         Cost_KOH        Cost of KOH /1.6/
*         Cost_MeOH       Cost of Methanol /3.85/                        (Old value is 0.28)
*         Cost_PhosA      Cost of Phosphoric acid /0.595/                 (Old value is 0.34)
*         Cost_Steam      Cost of steam /0.02/
*         Cost_Wa         Cost of process water /0.0004/                 (Old value if 0.0105)
*         Cost_Oil        Cost of soybean oil /0.8147/                          (Old value is 0.682)
*         Annual          Factor for converting Kg per s to kg per year /31622400/
*         Cost_tree       Restoration cost /196.34/
*         Cost_land       Land cost       /3300/ ;


         Cost_FAME       Cost of biodiesel /1.4/
         Cost_Glycerol   Cost of glycerol /1.026/
         Cost_PotPhos    Cost of potassium phosphate /3.491/
         Cost_KOH        Cost of KOH /1.6/
         Cost_MeOH       Cost of Methanol /3.85/
         Cost_PhosA      Cost of Phosphoric acid /0.595/
         Cost_Steam      Cost of steam /0.02/
         Cost_Wa         Cost of process water /0.002/
         Cost_Oil        Cost of soybean oil /0.682/
         Annual          Factor for converting Kg per s to kg per year /21168000/
         Cost_tree       Restoration cost /196.34/
         Cost_land       Land cost       /3300/ ;
*All costs are in $/kg. Check notes for references
*Replace cost of water with (0.0105,0.021) for extreme case




*         Rel_1, Rel_2;


*Obj..    Z =E= Q('Hx1') + Q('Jet1') + Q('Hx4') + Q('Hx6') + Q('Hx7') + Q('Hx8') + Q('Dry1') + Q_reb('BC1') + Q_reb('Rec1') + Q('Hx3');


*Rel_1(J,unit,unit1)$Arc(unit,unit1)..
*         fc(J,unit,unit1) =e= F(unit,unit1)*x(J,unit,unit1);

*Rel_2(unit,unit1)$Arc(unit,unit1)..
*         Sum(J,fc(J,unit,unit1)) =e= F(unit,unit1);


Equations
         Feed_1, Feed_2, Admixing_1;

Feed_1..
         fc('MeOH','Mix2','Mix3') =e= ratio_met_alk*(0.034826087)*fc('Oil','Src','Mix3');
*Calculating flow of methanol between the mix streams
*MW('MeOH')/MW('Oil') = 0.034826087

Feed_2..
         fc('MeOH','Col3','Mix2') =e= fc('MeOH','Mix2','Mix3') - fc('MeOH','Src3','Mix2');
*Calculating recycle stream flow
Admixing_1(J)..
         fc(J,'Mix3','HX45') =e= fc(J,'Mix2','Mix3') + fc(J,'Src','Mix3') + fc(J,'Src5','Mix3');


*Reaction transesterification
Equations
         Reactran_1, HX45_1, Yield, Convert, Reactran_2, Reactran_3, Reactran_4, Reactran_5, Reactran_6, Reactran_7, Reactran_8, Reactran_9, React3_1,Reactran_10,
         Reactran_11;


HX45_1..
         Q('HX45') =E= Sum(J,fc(J,'Mix3','HX45')*c_p_ind(J))*(T('HX45','Reactor3')-T('Mix3','HX45'));
Reactran_1(J)..
         fc(J,'HX45','Reactor3') =E= fc(J,'Mix3','HX45');

*Op refers to yield here
Yield..
         Op =E= (74.6301 + (0.4209*Tmp_reactor) + (15.1582*cat) + (3.1561*ratio_met_alk)
         -(0.0019*T('HX45','Reactor3')*T('HX45','Reactor3')) - (0.2022*T('HX45','Reactor3')*cat)
         -(0.01925*T('HX45','Reactor3')*ratio_met_alk) - (4.0143*(cat)*(cat))
         -(0.3400*cat*ratio_met_alk) - (0.1459*(ratio_met_alk)*(ratio_met_alk)));

Convert..
         conver =E= Op*0.01;


Reactran_2..
         fc('FAME','Reactor3','HX10') =E= fc('FAME','HX45','Reactor3')
         + 3*(0.317608696 )*conver*fc('Oil','HX45','Reactor3')
         + (1.04171123)*conver*fc('FFA','HX45','Reactor3');

*MW('FAME')/MW('Oil') = 0.317608696
*MW('FAME')/MW('FFA')= 1.04171123


Reactran_3..
         fc('Oil','Reactor3','HX10') =E=(1-conver)*fc('Oil','HX45','Reactor3');


Reactran_4..
         fc('FFA','Reactor3','HX10') =E=(1-conver)*fc('FFA','HX45','Reactor3');

Reactran_5..
         fc('MeOH','Reactor3','HX10') =E= fc('MeOH','HX45','Reactor3')
         -(conver*fc('Oil','HX45','Reactor3'))*3*(0.034826087)
         -(conver*fc('FFA','HX45','Reactor3'))*(0.114224599);
*MW('MeOH')/MW('Oil') = 0.034826087
*MW('MeOH')/MW('FFA') = 0.114224599




Reactran_6..
         fc('Glycerol','Reactor3','HX10') =E= conver*fc('Oil','HX45','Reactor3')*(0.100101087);
*MW('Glycerol')/MW('Oil') = 0.100101087


Reactran_7..
         fc('Wa','Reactor3','HX10') =E= fc('Wa','HX45','Reactor3')
         + conver*fc('FFA','HX45','Reactor3')*(0.064224599);
*MW('Wa')/MW('FFA') = 0.064224599

Reactran_8..
         fc('KOH','Reactor3','HX10') =E= fc('KOH','HX45','Reactor3');


Reactran_9..
         T('Reactor3','HX10') =E= T('HX45','Reactor3');

React3_1..
         Q('Reactor3') =E= 469*(fc('FAME','Reactor3','HX10') - fc('FAME','HX45','Reactor3'));

Reactran_10..
         fc('KOH','Mix2','Mix3') =L= 0.1;
*         fc('KOH','Src5','Mix3') =E= Cat*fc('Oil','Src','Mix3')*0.01;
Reactran_11..
        fc('Oil','Mix2','Mix3') =L= 0.4;
*Reactran_12..
*         fc('FFA','Mix2','Mix3') =L=0.1;

*****************************
*Methanol Distillation*
*****************************
Positive variables
         m_frac_c3(J,unit,unit1)    mol fraction of components in Column3
         recover_c3_MeOH            Recovery of methanol in column 3
         Reflux                     Intermediate: Reflux ratio in column 3
         R_Col3                     Reflux ratio in column 3
         P_Col3                     Pressure at feed tray in col3
         VpMeOH_1                   Vp of methanol in MD
         VpGlycerol_1               Vp of glycerol in MD
         VpKOH_1                    Vp of KOH in MD
         VpWa_1                     Vp of water in MD
         VpFAME_1                   Vp of FAME in MD
         P_Col3_Temp                Calculating pressure in col3
         Col3Mix2Temp               Temperature of stream Col3 to Mix2
         VpOil_1                    Vapor pressure of oil
         recover_c3_Wa              Recovery of water in col3
         recover_c3_FAME            Recovery of FAME in Col3
         recover_c3_Glycerol        Recovery of glycerol in Col3
         recover_c3_KOH             Recovery of KOH in col3
         recover_c3_Oil             Recovery of oil in col3
         recover_c3_FFA             Recovery of FFA in Col3
         R_Col3_Temp                Temp factor for R_Col3;
Equations
         MD_1,MD_2,MD_3,MD_4,MD_5, MD_6, MD_7, MD_8,MD_9, MD_10, MD_11, MD_12, MD_13, MD_14, MD_15, MD_16, MD_17, MD_31, MD_32,MD_33,MD_34, MD_35, MD_36
         MD_18, MD_19, MD_20, MD_21, MD_23, MD_24, MD_25, MD_26 , HX11_1, HX12_1, MD_27, MD_28, MD_29, MD_30,MD_37, MD_38;
*Scalar
*         dp              Pressure loss in column /0.1/
*         alpha_c3        Relative volatility in column /4.5/;


recover_c3_MeOH.lo = 0.94;
recover_c3_MeOH.up = 0.99;
recover_c3_MeOH.l = (recover_c3_MeOH.lo + recover_c3_MeOH.up)/2;
P_Col3_Temp.lo = 0.01;
recover_c3_Wa.lo = 0.75;
recover_c3_Wa.up = 0.99;
recover_c3_Wa.l = (recover_c3_Wa.lo + recover_c3_Wa.up)/2;
recover_c3_FAME.lo = 0.01;
recover_c3_FAME.up = 0.99;
recover_c3_FAME.l = (recover_c3_FAME.up + recover_c3_FAME.lo)/2;
recover_c3_Glycerol.lo = 0.75;
recover_c3_Glycerol.up = 0.99;
recover_c3_Glycerol.l = (recover_c3_Glycerol.lo + recover_c3_Glycerol.up)/2;
recover_c3_KOH.lo = 0.75;
recover_c3_KOH.up = 0.99;
recover_c3_KOH.L = (recover_c3_KOH.lo + recover_c3_KOH.up)/2;
recover_c3_Oil.lo = 0.75;
recover_c3_Oil.up = 0.99;
recover_c3_Oil.l = (recover_c3_Oil.lo + recover_c3_Oil.up)/2;
recover_c3_FFA.lo = 0.01;
recover_c3_FFA.up = 0.99;
recover_c3_FFA.l = (recover_c3_FFA.lo  + recover_c3_FFA.up)/2;


*for all components J in Wa,MeOH,Glycerol,KOH
MD_1..
         fc('Wa','HX10','Col3') =E= fc('Wa','Reactor3','HX10');
MD_2..
         fc('MeOH','HX10','Col3')=E= fc('MeOH','Reactor3','HX10');
MD_3..
         fc('Glycerol','HX10','Col3')=E= fc('Glycerol','Reactor3','HX10');
MD_4..
         fc('KOH','HX10','Col3')=E= fc('KOH','Reactor3','HX10');
MD_5..
         fc('Oil','HX10','Col3')=E= fc('Oil','Reactor3','HX10');

MD_6..
         fc('FAME','HX10','Col3')=E= fc('FAME','Reactor3','HX10');

MD_38..
         fc('FFA','HX10','Col3') =E= fc('FFA','Reactor3','HX10');


MD_7..
         F('HX10','Col3') =E=  fc('MeOH','HX10','Col3') + fc('Glycerol','HX10','Col3')
                               + fc('KOH','HX10','Col3') + fc('Wa','HX10','Col3') ;



MD_8(J)$((ord(J) ne 8) and (ord(J) ne 9) and (ord(J) ne 5) and (ord(J) ne 6) and (ord(J) ne 7))..
         m_frac_c3(J,'HX10','Col3') =E= ((fc(J,'HX10','Col3')/F('HX10','Col3')/MW(J)))
         /(((fc('MeOH','HX10','Col3')/F('HX10','Col3'))/MW('MeOH'))
         + ((fc('Glycerol','HX10','Col3')/F('HX10','Col3'))/MW('Glycerol'))
         + ((fc('KOH','HX10','Col3')/F('HX10','Col3'))/MW('KOH'))
         + ((fc('Wa','HX10','Col3')/F('HX10','Col3'))/MW('Wa')));

*+ ((fc('Oil','HX10','Col3')/F('HX10','Col3'))/MW('Oil')));



MD_9..
         VpFAME_1 =E=  (10**(coef_p('FAME','1')-(coef_p('FAME','2')/
         (coef_p('FAME','3')+T('HX10','Col3')+273))))*0.0075;

*Conversion from Pa to mmHG
MD_10..
         VpWa_1 =E= (10**(coef_p('Wa','1')-(coef_p('Wa','2')/
         (coef_p('Wa','3')+(T('HX10','Col3')+273)))))*750.06;

*Conversion from Bar to mmHg
MD_11..
         VpGlycerol_1 =E= 10**(coef_p('Glycerol','1')-(coef_p('Glycerol','2')/T('HX10','Col3')) -
           (coef_p('Glycerol','3')*log((T('HX10','Col3')+273))));
*Torr and mmhg are almost equal
*** NOTE ANTOINE EQUATION FOR OIL AND GLYCEROL ARE SLIGHTLY DIFF FROM THE ORIGINAL EQUATION **
MD_12..
         VpMeOH_1 =E= (10**(coef_p('MeOH','1')-(coef_p('MeOH','2')/
         (coef_p('MeOH','3')+(T('HX10','Col3')+273)))))*750.06;
*Conversion from Bar to mmHg


MD_13..
         VpKOH_1 =E= (10**(coef_p('KOH','1')-(coef_p('KOH','2')/
         (coef_p('KOH','3')+(T('HX10','Col3'))))));

*T in C and Vp in mmHG


MD_14..
         VpOil_1 =E= (exp(coef_p('Oil','1')+ (coef_p('Oil','2')/
         (coef_p('Oil','3')+(T('HX10','Col3')+273)))))*750.06;

*** NOTE ANTOINE EQUATION FOR OIL AND GLYCEROL ARE SLIGHTLY DIFF FROM THE ORIGINAL EQUATION **
*Oil is in kPA, Temp is in K
MD_15..
         P_Col3 =E= VpFAME_1
         + (m_frac_c3('Wa','HX10','Col3')*VpWa_1)
         + (m_frac_c3('MeOH','HX10','Col3')*VpMeOH_1)
         + (m_frac_c3('Glycerol','HX10','Col3')*VpGlycerol_1)
         + (m_frac_c3('KOH','HX10','Col3')*VpKOH_1);
*+ (m_frac_c3('Oil','HX10','Col3')*VpOil_1);

*         + (m_frac_c3('MeOH','HX10','Col3')*exp(coef_p('MeOH','1')-(coef_p('MeOH','2')/(coef_p('MeOH','3')+T('HX10','Col3')))))
*         + (m_frac_c3('Glycerol','HX10','Col3')*exp(coef_p('Glycerol','1')-(coef_p('Glycerol','2')/(coef_p('Glycerol','3')+T('HX10','Col3')))))
*         + (m_frac_c3('KOH','HX10','Col3')*exp(coef_p('KOH','1')-(coef_p('KOH','2')/(coef_p('KOH','3')+T('HX10','Col3')))));


*P_Col3 in mmhg


MD_16(J)..
        Q('HX10') =E= sum(h,fc(h,'HX10','Col3')*c_p_ind(h))*(T('HX10','Col3')-T('Reactor3','HX10'));
*This has highest energy/heat consumption. Does not make any sense!

$ontext
This part is slighlty modified than the python code. In python, using flow of each individual components through
reactros for the mass balance. In GAMS, using "J" of all components to do mass balance, same as what is there in the supporting information
Also, using recovery of MeOH alone since there is no recovery factor given for the rest of the components.
$offtext


MD_17..
         fc('Wa','Col3','HX13') =E= recover_c3_Wa*fc('Wa','HX10','Col3');

MD_37..
         fc('MeOH','Col3','HX13') =E= (1-recover_c3_MeOH)*fc('MeOH','HX10','Col3');

MD_32..
         fc('Glycerol','Col3','HX13') =E= recover_c3_Glycerol*fc('Glycerol','HX10','Col3');

MD_33..
         fc('KOH','Col3','HX13') =E= recover_c3_KOH*fc('KOH','HX10','Col3');

MD_34..
         fc('Oil','Col3','HX13') =E= recover_c3_Oil*fc('Oil','HX10','Col3');

MD_35..
         fc('FFA','Col3','HX13') =E= recover_c3_FFA*fc('FFA','HX10','Col3');

MD_36..
         fc('FAME','Col3','HX13') =E= recover_c3_FAME*fc('FAME','HX10','Col3');


MD_18..
         fc('MeOH','Col3','Mix2') =E= recover_c3_MeOH*fc('MeOH','HX10','Col3');


MD_19..
         (P_Col3_Temp) =E= (1-0.5*dp)*P_Col3*0.0013332239;

*P_Col3_Temp is mmhg. Convert to bar "x" b'Col3'y 0.0013332239 since coef of MeOH are in bar and K
MD_20..
                 (Col3Mix2Temp - 33.5)*(5.20409 - log10(P_Col3_Temp)) =E= 1581.341;
*Col3Mix2Temp =E= (1581.341/(5.20409 - log10(P_Col3_Temp)))  -  (-33.50);


MD_21..
         T('Col3','Mix2') =E= Col3Mix2Temp - 273;


*Temp is in C

*(coef_p('MeOH','2')/((coef_p('MeOH','1') -  log((1 - 0.5*dp)*P_Col3)))) - coef_p('MeOH','3');
*log((1 - 0.5*dp)*P_Col3) =E= coef_p('MeOH','1') - ((coef_p('MeOH','2')/(coef_p('MeOH','3') + T('Col3','Mix2'))));
*T('Col3','Mix2') =e= -(coef_p('MeOH','3')) + (-coef_p('MeOH','2')/(log((1-(0.5*dp))*P_Col3)) - coef_p('MeOH','1'));



MD_23..
         F('Col3','HX13') =E= fc('MeOH','Col3','HX13') + fc('Glycerol','Col3','HX13') + fc('KOH','Col3','HX13') + fc('Wa','Col3','HX13')  ;


MD_24(J)$((ord(J) ne 8) and (ord(J) ne 9) and (ord(J) ne 6) and (ord(J) ne 5) and (ord(J) ne 7))..
         m_frac_c3(J,'Col3','HX13') =E= ((fc(J,'Col3','HX13')/F('Col3','HX13')/MW(J)))
         /(((fc('MeOH','Col3','HX13')/F('Col3','HX13'))/MW('MeOH'))
         + ((fc('Glycerol','Col3','HX13')/F('Col3','HX13'))/MW('Glycerol'))
         + ((fc('KOH','Col3','HX13')/F('Col3','HX13'))/MW('KOH'))
         + ((fc('Wa','Col3','HX13')/F('Col3','HX13'))/MW('Wa')));

*+ ((fc('Oil','Col3','HX13')/F('Col3','HX13'))/MW('Oil'))
*+ ((fc('FAME','Col3','HX13')/F('Col3','HX13'))/MW('FAME')));



MD_25..
*(J)$((ord(J) ne 8) and (ord(J) ne 9) and (ord(J) ne 1) and (ord(J) ne 6))..
         P_Col3_Temp =E= (10**(coef_p('FAME','1')-coef_p('FAME','2')/
         (coef_p('FAME','3')+T('Col3','HX13'))))*0.0075
         + (m_frac_c3('MeOH','Col3','HX13')*(10**(coef_p('MeOH','1')-(coef_p('MeOH','2')/(coef_p('MeOH','3')+T('Col3','HX13')))))*750.06)
         + (m_frac_c3('KOH','Col3','HX13')*(10**(coef_p('KOH','1')-coef_p('KOH','2')/(coef_p('KOH','3')+(T('Col3','HX13') + 273)))))
         + (m_frac_c3('Wa','Col3','HX13')*(10**(coef_p('Wa','1')-(coef_p('Wa','2')/(coef_p('Wa','3')+(T('Col3','HX13'))))))*750.06)
         + (m_frac_c3('Glycerol','Col3','HX13')*(10**(coef_p('Glycerol','1')-(coef_p('Glycerol','2')/T('Col3','HX13')) - (coef_p('Glycerol','3')*log(T('Col3','HX13'))))));


*** NOTE ANTOINE EQUATION FOR OIL AND GLYCEROL ARE SLIGHTLY DIFF FROM THE ORIGINAL EQUATION **


**** Check this line. Probably calculating summation this way may not be right ***

MD_26..
         R_Col3 =E= 1.5*(1/(alpha_c3 - 1))*((0.9999/(m_frac_c3('MeOH','HX10','Col3')+0.001))- (alpha_c3*((1-0.999)/(1-m_frac_c3('MeOH','HX10','Col3')))));

R_Col3.lo = 1;
R_Col3.up = 3;
*R_Col3.l = (R_Col3.lo + R_Col3.up)/2;


HX11_1..
         Q('HX11') =E= -fc('MeOH','Col3','Mix2')*(R_Col3 + 1)*dH_vap_0('MeOH');
*Only methanol is recycled

MD_27..
          x('Wa','Col3','HX13')*F('Col3','HX13') =E= fc('Wa','Col3','HX13');
*x('Wa','Col3','HX13') =E= fc('Wa','Col3','HX13')/F('Col3','HX13');


MD_28..
           x('Glycerol','Col3','HX13')*F('Col3','HX13') =E= fc('Glycerol','Col3','HX13');
*x('Glycerol','Col3','HX13') =E= fc('Glycerol','Col3','HX13')/F('Col3','HX13');

MD_29..
           x('MeOH','Col3','HX13')*F('Col3','HX13') =E= fc('MeOH','Col3','HX13');
*x('MeOH','Col3','HX13') =E= fc('MeOH','Col3','HX13')/F('Col3','HX13');

MD_30..
          x('FAME','Col3','HX13')*F('Col3','HX13') =E= fc( 'FAME','Col3','HX13');
*x('FAME','Col3','HX13') =E= fc( 'FAME','Col3','HX13')/F('Col3','HX13');

MD_31..
          x('Oil','Col3','HX13')*F('Col3','HX13') =E= fc('Oil','Col3','HX13');
*x('Oil','Col3','HX13') =E= fc('Oil','Col3','HX13')/F('Col3','HX13');


HX12_1..
         Q('HX12') =E= F('Col3','HX13')*(R_Col3 + 1)*((x('Wa','Col3','HX13')*dH_vap_0('Wa')) +
               (x('Glycerol','Col3','HX13')*dH_vap_0('Glycerol')) +
               (x('MeOH','Col3','HX13')*dH_vap_0('MeOH')) +
               (x('FAME','Col3','HX13')*dH_vap_0('FAME')) +
               (x('Oil','Col3','HX13')*dH_vap_0('Oil')));

*****************************
********** WASHING **********
*****************************
*Src6 has water as input which is a variable
**##################################################****
*fc.lo('Wa','Src6','Sep2') = 0.01;
*fc.up('Wa','Src6','Sep2') = 10;
**##################################################****
*fc.lo('Wa','SrcWa','Src6') = 0.01;
*fc.up('Wa','SrcWa','Src6') = 10;

*fc.l('Wa','SrcWa','Src6') = (fc.lo('Wa','SrcWa','Src6') + fc.up('Wa','SrcWa','Src6'))/2;
*fc.L


Positive variable
         water_req      Water input required;


Equations
         Wash_1, Wash_2, HX13_1 ,Wash_3, Wash_4, Wash_5, Wash_6, Wash_7, Wash_8, Wash_9, Wash_10, Wash_11, Wash_12, Wash_14, Wash_15, Wash_16, Wash_17, Wash_18, Wash_19
         Wash_20, Wash_21, Wash_22, Wash_23, Wash_13,Wash_25, Wash_26 ;

Wash_1..
         fc('Wa','HX13','Sep2') =E= fc('Wa','Col3','HX13');
Wash_2..
         fc('MeOH','HX13','Sep2') =E= fc('MeOH','Col3','HX13');
Wash_3..
         fc('Glycerol','HX13','Sep2') =E= fc('Glycerol','Col3','HX13');
Wash_4..
         fc('KOH','HX13','Sep2') =E= fc('KOH','Col3','HX13');
Wash_5..
         fc('Oil','HX13','Sep2') =E= fc('Oil','Col3','HX13');
Wash_6..
         fc('FAME','HX13','Sep2') =E= fc('FAME','Col3','HX13');
Wash_7..
         fc('FFA','HX13','Sep2') =E= fc('FFA','Col3','HX13');
Wash_9..
         F('HX13','Sep2') =E= fc('Wa','Col3','HX13') + fc('MeOH','Col3','HX13')+ fc('Glycerol','Col3','HX13')+ fc('KOH','Col3','HX13')+ fc('Oil','Col3','HX13')+
         fc('FAME','Col3','HX13') + fc('FFA','Col3','HX13');

**** THIS LINE IS NOT WORKING. NOT ABLE TO ASSIGN A VARIABLE TP FC.UP. CHECK HOW THIS CAN BE FIXED *****
*Value of F('HX13','Sep2') is 21.15 so fix water input max as 10
*fc.UP('Wa','Src6','Sep2') = 0.01*F.l('HX13','Sep2');

Wash_8..
         water_req =E= 0.01*F('HX13','Sep2');
*Water input is 1% of total biodiesel phase.
Wash_26..
         fc('Wa','Src6','Sep2') =E= water_req;

HX13_1..
         Q('HX13') =E=(fc('Wa','Col3','HX13')*c_p_ind('Wa') + fc('MeOH','Col3','HX13')*c_p_ind('MeOH') +
                      fc('Glycerol','Col3','HX13')*c_p_ind('Glycerol') + fc('KOH','Col3','HX13')*c_p_ind('KOH') +
         fc('FAME','Col3','HX13')*c_p_ind('FAME') + fc('FFA','Col3','HX13')*c_p_ind('FFA') + fc('Oil','Col3','HX13')*c_p_ind('Oil'))*(T('HX13','Sep2')-(T('Col3','HX13')));

Wash_10..
         T('Sep2','HX19') =E= T('Sep2','Reactor4');

Wash_11..
         fc('Wa','Sep2','Reactor4') =E= Separation('Wa')*(fc('Wa','HX13','Sep2') +fc('Wa','Src6','Sep2'));
Wash_12..
         fc('MeOH','Sep2','Reactor4') =E= Separation('MeOH')*(fc('MeOH','HX13','Sep2') + fc('MeOH','Src6','Sep2'));
Wash_13..
         fc('Glycerol','Sep2','Reactor4') =E= Separation('Glycerol')*(fc('Glycerol','HX13','Sep2') + fc('Glycerol','Src6','Sep2'));
Wash_14..
         fc('KOH','Sep2','Reactor4') =E= Separation('KOH')*(fc('KOH','HX13','Sep2') + fc('KOH','Src6','Sep2'));
Wash_25..
         fc('Oil','Sep2','Reactor4') =E= Separation('Oil')*(fc('Oil','HX13','Sep2') + fc('Oil','Src6','Sep2'));
Wash_15..
         fc('FAME','Sep2','Reactor4') =E= Separation('FAME')*(fc('FAME','HX13','Sep2') + fc('FAME','Src6','Sep2'));
Wash_16..
         fc('FFA','Sep2','Reactor4') =E= Separation('FFA')*(fc('FFA','HX13','Sep2') + fc('FFA','Src6','Sep2'));
Wash_17..
         fc('Wa','Sep2','HX19') =E= (1-Separation('Wa'))*(fc('Wa','HX13','Sep2') +fc('Wa','Src6','Sep2'));
Wash_18..
         fc('MeOH','Sep2','HX19') =E= (1-Separation('MeOH'))*(fc('MeOH','HX13','Sep2') + fc('MeOH','Src6','Sep2'));
Wash_19..
         fc('Glycerol','Sep2','HX19') =E= (1-Separation('Glycerol'))*(fc('Glycerol','HX13','Sep2') + fc('Glycerol','Src6','Sep2'));
Wash_20..
         fc('KOH','Sep2','HX19') =E= (1-Separation('KOH'))*(fc('KOH','HX13','Sep2') + fc('KOH','Src6','Sep2'));
Wash_21..
         fc('Oil','Sep2','HX19') =E= (1-Separation('Oil'))*(fc('Oil','HX13','Sep2') + fc('Oil','Src6','Sep2'));
Wash_22..
         fc('FAME','Sep2','HX19') =E= (1-Separation('FAME'))*(fc('FAME','HX13','Sep2') + fc('FAME','Src6','Sep2'));
Wash_23..
         fc('FFA','Sep2','HX19') =E= (1-Separation('FFA'))*(fc('FFA','HX13','Sep2') + fc('FFA','Src6','Sep2'));

*******************************
** Neutralization Alkali reaction **
*******************************

Equations
         Neualk_1, Neualk_2, Neualk_3, Neualk_4, Neualk_5, Reactor4_1,  Neualk_6,  Neualk_7,  Neualk_8,  Neualk_9,  Neualk_10;

Neualk_1..
         fc('MeOH','Reactor4','HX14') =E= fc('MeOH','Sep2','Reactor4');
Neualk_6..
         fc('Glycerol','Reactor4','HX14') =E= fc('Glycerol','Sep2','Reactor4');
Neualk_7..
         fc('Oil','Reactor4','HX14') =E= fc('Oil','Sep2','Reactor4');
Neualk_8..
         fc('FFA','Reactor4','HX14') =E= fc('FFA','Sep2','Reactor4');
Neualk_9..
         fc('FAME','Reactor4','HX14') =E= fc('FAME','Sep2','Reactor4');
Neualk_10..
         fc('PhosA','Reactor4','HX14') =E= fc('PhosA','Sep2','Reactor4');

Neualk_2..
         fc('PhosA','Src7','Reactor4') =E= 0.333*fc('KOH','Sep2','Reactor4')*(1.746706211);
*MW('PhosA')/MW('KOH') = 1.746706211
Neualk_3..
         fc('PotPhos','Reactor4','Snk3') =E= 0.333*fc('KOH','Sep2','Reactor4')*(3.104502937);
*MW('PotPhos')/MW('KOH')  =3.104502937
Neualk_4..
         fc('Wa','Reactor4','HX14') =E= fc('Wa','Sep2','Reactor4') + fc('KOH','Sep2','Reactor4')*(0.321090943);
*MW('Wa')/MW('KOH') =     0.321090943

Reactor4_1..
         Q('Reactor4') =E= (fc('PhosA','Src7','Reactor4')/MW('PhosA'))*dhreactneu;

Neualk_5(J)$((ord(J) eq 1) and (ord(J)eq 2) and (ord(J)eq 3) and (ord(J)eq 4) and (ord(J)eq 8) and (ord(J)eq 9))..
         Q('Reactor4') =E= sum(h,fc(h,'Reactor4','HX14')*c_p_ind(h))*(T('Reactor4','HX14')-(T('Sep2','Reactor4')));
*Wa, MeOH, Gly, KOH, PhosA, PotPhoS
*Solve for T('Reactor4','HX14')


*******************************
*** FAME DISTILLATION ********
***************************
Positive variables
         m_frac_c5       Fraction of components in column 5
         VpOil_2         Vapor pressure of oil in Col 5
         VpFAME_2        Vapor pressure of FAME in Col 5
         VpFFA_2         Vapor pressure of FFA in Col 5
         P_Col5          Pressure in Col5
         Col5HX22Temp    Distillate temperature_temp var
         P_Col5_Temp_1   Distillate Pressure in Col5
         recover_oil     Oil recovery in FAME distaillation
         recover_FAME    FAME recovery in distillation column
         P_Col5_Temp_2   Reboiler pressure in Col5
         Col5HX23Temp    Bubble point temperature_temp var
         R_Col5          Reflux ratio in distillation;


recover_oil.lo = 0.5;
recover_oil.up = 0.997;
recover_oil.l = (recover_oil.lo + recover_oil.up)/2;
recover_FAME.lo = 0.5;
recover_FAME.up = 0.997;
recover_FAME.l = (recover_FAME.lo + recover_FAME.up)/2;
P_Col5_Temp_1.lo = 0.01;
P_Col5_Temp_2.lo = 0.01;

*** This is another variable in the design problem ***
R_Col5.lo = 2;
R_Col5.up = 3;
R_Col5.l = (R_Col5.lo + R_Col5.up)/2;

Equations
         FAME_1, FAME_2, FAME_3, HX19_1, FAME_4, FAME_5, FAME_6, FAME_7, FAME_8, FAME_9,FAME_10, FAME_11, FAME_12, FAME_13, FAME_14, FAME_15, FAME_16, FAME_17,
         FAME_18, FAME_19, FAME_21, FAME_22, HX20_1, HX21_1, HX23_1, HX22_1, FAME_23;


FAME_1..
         fc('Oil','HX19','Col5') =E= fc('Oil','Sep2','HX19');

FAME_2..
         fc('FAME','HX19','Col5') =E= fc('FAME','Sep2','HX19');

FAME_3..
         fc('FFA','HX19','Col5') =E= fc('FFA','Sep2','HX19');

HX19_1(J)$((ord(J) eq 5) and (ord(J)eq 6) and (ord(J)eq 7))..
         Q('HX19') =E= sum(h,fc(h,'HX19','Col5')*c_p_ind(h))*(T('HX19','Col5')-(T('Sep2','HX19')));


FAME_4..
         F('HX19','Col5') =E=  fc('Oil','HX19','Col5') + fc('FAME','HX19','Col5') +  fc('FFA','HX19','Col5');


FAME_5(J)$((ord(J) eq 5) and (ord(J) eq 6) and (ord(J) eq 7))..
         m_frac_c5(J,'HX19','Col5') =E= ((fc(J,'HX19','Col5')/F('HX19','Col5')/MW(J)))
         /(((fc('Oil','HX19','Col5')/F('HX19','Col5'))/MW('Oil'))
         + ((fc('FAME','HX19','Col5')/F('HX19','Col5'))/MW('FAME'))
         + ((fc('FFA','HX10','Col5')/F('HX19','Col5'))/MW('FFA')));

FAME_6..
         VpOil_2 =E= (exp(coef_p('Oil','1')+ (coef_p('Oil','2')/
         (coef_p('Oil','3')+(T('HX19','Col5')+273)))))*0.7506;


*Vp in kPa is converted to mmhg, Temp in C is converted to K

FAME_7..
         VpFAME_2 =E=  (10**(coef_p('FAME','1')-(coef_p('FAME','2')/
         (coef_p('FAME','3')+T('HX19','Col5')+273))))*0.0075;

*Conversion from Pa to mmHG

FAME_8..
         VpFFA_2 =E= (10**(coef_p('FFA','1')-(coef_p('FFA','2')/
         (coef_p('FFA','3')+(T('HX19','Col5')+273)))))*750.06;
*Conversion from Bar to mmhg

FAME_9..
         P_Col5 =E= (m_frac_c5('Oil','HX19','Col5')*VpOil_2) + (m_frac_c5('FAME','HX19','Col5')*VpFAME_2)
                 + (m_frac_c5('FFA','HX19','Col5')*VpFFA_2);

FAME_10..
         fc('FAME','Col5','HX22') =E= recover_FAME*fc('FAME','HX19','Col5');
FAME_11..
         fc('FFA','Col5','HX22') =E= fc('FFA','HX19','Col5');


FAME_12..
         fc('Oil','Col5','HX22') =E= (1-recover_oil)*fc('Oil','HX19','Col5');


FAME_13..
         fc('FAME','Col5','HX23') =E= (1-recover_FAME)*fc('FAME','HX19','Col5');


FAME_14..
         fc('Oil','Col5','HX23') =E= (recover_Oil)*fc('Oil','HX19','Col5');

FAME_15..

         P_Col5_Temp_1 =E= (1-0.5*dp)*P_Col5*133.322365;
*Pressure from mmHG to Pa since FAME is PA

FAME_16..
         Col5HX22Temp =E= (1450.62/(8.2175 - log10(P_Col5_Temp_1)))  -  (88.03);

*Temp is in K
FAME_17..
         T('Col5','HX22') =E= Col5HX22Temp - 273;


FAME_18..
         P_Col5_Temp_2 =E= (1+0.5*dp)*P_Col5*0.133322365;
FAME_23..
         fc('FAME','Col5','HX22') =E= 0.5318004825;

*Min production of 5 million gals per year
*FAME_24..
*         fc('PotPhos','Reactor4','Snk3') =L= 0.532;
*Constraint that PotPhos production should be less than 0.532.

FAME_19..
        T('Col5','HX23')  =E= (-(-167.48) - (-708.72/(log(P_Col5_Temp_2) - 11.4785)));

*Pressure is in mmHg converted to kPa since VP of oil is in KPA
*Converting Temp from K to C

FAME_21..
         F('Col5','HX22') =E= fc('Oil','Col5','HX22') + fc('FAME','Col5','HX22') + fc('FFA','Col5','HX22');

FAME_22..
         F('Col5','HX23') =E= fc('Oil','Col5','HX23') + fc('FAME','Col5','HX23') + fc('FFA','Col5','HX23');

HX20_1..
         Q('HX20') =E= -F('Col5','HX22')*(R_Col5 + 1)*dH_vap_0('FAME');

HX21_1..
         Q('HX21') =E= F('Col5','HX23')*(R_Col5 + 1)*dH_vap_0('Oil');



HX22_1(J)$((ord(J) eq 5) and (ord(J)eq 6) and (ord(J)eq 7))..
         Q('HX22') =E= (fc('Oil','Col5','HX22')*c_p_ind('Oil') + fc('FFA','Col5','HX22')*c_p_ind('FFA')+ fc('FAME','Col5','HX22')*c_p_ind('FAME'))*(T_amb-T('Col5','HX22'));

HX23_1(J)$((ord(J) eq 5) and (ord(J)eq 6) and (ord(J)eq 7))..
         Q('HX23') =E=(fc('Oil','Col5','HX22')*c_p_ind('Oil') + fc('FFA','Col5','HX22')*c_p_ind('FFA')+ fc('FAME','Col5','HX22')*c_p_ind('FAME'))*(T_amb-(T('Col5','HX23')));



*******************************
**** Glcyerol Distillation *****
*********************************

Equations
         Glydis_1, Glydis_2, Glydis_3, HX14_1, Glydis_4, Glydis_5, Glydis_6, Glydis_7, Glydis_8, Glydis_9, Glydis_10, Glydis_11, Glydis_12, Glydis_13,Glydis_14,
          Glydis_15, Glydis_16, Glydis_17, Glydis_18, Glydis_19 , Glydis_20, Glydis_21, Glydis_22, Glydis_23, Glydis_24, Glydis_25, Glydis_26, HX15_1, HX16_1
         HX17_1, HX18_1, Glydis_27, Glydis_28;

Positive variables
         m_frac_c4       Fraction of components in column 4
         P_Col4          Pressure in Dis col 4
         VpWa_3          Vapor pressure of water
         VpGly_3         Vapor pressure of glycerol
         VpMeOH_3        Vapor preesure of MeOH
         recover_c4_wa   Water recovery in col4
         x_gly           purity of glycerol
         VpWa_4          Vp of water in col4
         VpMeOH_4        Vp of methanol in col4
         VpGly_4         Vp of glyercol in col4
         Col4HX17Temp    Temp in Col4
         P_Col4_Temp_1   Distillate pressure in col 4
         P_Col4_Temp_2   Reboiler pressure in Col 4
         Col4HX18Temp    Temp in HX18Col4 stream
         R_Col4          Reflux ratio in Column 4;

**** This is another variable in the problem ****
recover_c4_wa.lo = 0.5;
recover_c4_wa.up = 0.9;
recover_c4_wa.l = (recover_c4_wa.lo + recover_c4_wa.up)/2;

*recover_c4_wa.l  = (recover_c4_wa.lo + recover_c4_wa.up)/2;
** This value is given **

*Giving some lower bounds to avoid division by zero *
P_Col4_Temp_1.lo = 0.01;
P_Col4_Temp_2.lo = 0.01;
P_Col4.lo = 0.01;
Col4HX17Temp.lo = 0.01;
Col4HX18Temp.lo = 0.01;




Glydis_1..
         fc('Wa','HX14','Col4') =E= fc('Wa','Reactor4','HX14');
Glydis_2..
         fc('MeOH','HX14','Col4') =E= fc('MeOH','Reactor4','HX14') ;
Glydis_3..
         fc('Glycerol','HX14','Col4') =E= fc('Glycerol','Reactor4','HX14') ;

HX14_1(J)$((ord(J) eq 1) and (ord(J) eq 2) and (ord(J) eq 3))..
         Q('HX14') =E=sum(h,fc(h,'HX14','Col4')*c_p_ind(h))*(T('HX14','Col4')-T('Reactor4','HX14'));


Glydis_4..
         F('HX14','Col4') =E= fc('Wa','HX14','Col4') + fc('MeOH','HX14','Col4') + fc('Glycerol','HX14','Col4');

Glydis_5(J)$((ord(J) eq 1) and (ord(J) eq 2) and (ord(J) eq 3))..
         m_frac_c4(J,'HX14','Col4') =E= ((fc(J,'HX14','Col4')/F('HX14','Col4')/MW(J)))
         /(((fc('Wa','HX14','Col4')/F('HX14','Col4'))/MW('Wa'))
         + ((fc('MeOH','HX14','Col4')/F('HX14','Col4'))/MW('MeOH'))
         + ((fc('Glycerol','HX14','Col4')/F('HX14','Col4'))/MW('Glycerol')));


Glydis_6..
         VpWa_3 =E= (10**(coef_p('Wa','1')-(coef_p('Wa','2')/
         (coef_p('Wa','3')+(T('HX14','Col4')+273)))))*750.06;

*Conversion from Bar to mmHg, converting temp from C to K

Glydis_7..
          VpGly_3 =E= 10**(coef_p('Glycerol','1')-(coef_p('Glycerol','2')/T('HX14','Col4')) -
           (coef_p('Glycerol','3')*log((T('HX14','Col4')+273))));
*Torr and mmhg are almost equal
*** NOTE ANTOINE EQUATION FOR OIL AND GLYCEROL ARE SLIGHTLY DIFF FROM THE ORIGINAL EQUATION **

Glydis_8..
          VpMeOH_3 =E= (10**(coef_p('MeOH','1')-(coef_p('MeOH','2')/
         (coef_p('MeOH','3')+(T('HX14','Col4')+273)))))*750.06;
*Conversion from Bar to mmHg


Glydis_9..
          P_Col4 =E= (m_frac_c4('Wa','HX14','Col4')*VpWa_3) + (m_frac_c4('MeOH','HX14','Col4')*VpMeOH_3) + (m_frac_c4('Gly','HX14','Col4')*VpGly_3);

*P in mmHg

Glydis_10..
         fc('MeOH','Col4','HX17') =E= fc('MeOH','HX14','Col4');

Glydis_11..
         fc('Glycerol','Col4','HX18') =E= fc('Glycerol','HX14','Col4');

Glydis_12..
         fc('Wa','Col4','HX17') =E= recover_c4_wa*fc('Wa','HX14','Col4');


Glydis_13..
         fc('Wa','Col4','HX18') =E= (1-recover_c4_wa)*fc('Wa','HX14','Col4');

Glydis_14..
         F('Col4','HX18') =E= fc('Wa','Col4','HX18') + fc('MeOH','Col4','HX18') + fc('Glycerol','Col4','HX18');

Glydis_15..
         x_gly =E= fc('Glycerol','Col4','HX18')/F('Col4','HX18');
Glydis_27..
         x_gly =G= 0.92;

Glydis_16..
         F('Col4','HX17') =E= fc('Wa','Col4','HX17')+ fc('MeOH','Col4','HX17') + fc('Glycerol','Col4','HX17');


Glydis_17(J)$((ord(J) eq 1) and (ord(J) eq 2) and (ord(J) eq 3))..
         m_frac_c4(J,'Col4','HX17') =E= ((fc(J,'Col4','HX17')/F('Col4','HX17')/MW(J)))
         /(((fc('Wa','Col4','HX17')/F('Col4','HX17'))/MW('Wa'))
         + ((fc('MeOH','Col4','HX17')/F('Col4','HX17'))/MW('MeOH'))
         + ((fc('Glycerol','Col4','HX17')/F('Col4','HX17'))/MW('Glycerol')));

Glydis_18..
         P_Col4_Temp_1 =E= (1-0.5*dp)*P_Col4;



Glydis_19..
         P_Col4_Temp_1 =E= (m_frac_c4('Wa','Col4','HX17')*((10**(coef_p('Wa','1')-(coef_p('Wa','2')/
         (coef_p('Wa','3')+(Col4HX17Temp+273)))))*750.06))   +
         (m_frac_c4('Glycerol','Col4','HX17')*10**(coef_p('Glycerol','1')-(coef_p('Glycerol','2')/Col4HX17Temp+273) -
           (coef_p('Glycerol','3')*log((Col4HX17Temp+273))))) +
         (m_frac_c4('MeOH','Col4','HX17')*((10**(coef_p('MeOH','1')-(coef_p('MeOH','2')/
         (coef_p('MeOH','3')+(Col4HX17Temp+273)))))*750.06));


Glydis_20..
         T('Col4','HX17') =E= Col4HX17Temp - 273;



Glydis_21(J)$((ord(J) eq 1) and (ord(J) eq 2) and (ord(J) eq 3))..
         m_frac_c4(J,'Col4','HX18') =E= ((fc(J,'Col4','HX18')/F('Col4','HX18')/MW(J)))
         /(((fc('Wa','Col4','HX18')/F('Col4','HX18'))/MW('Wa'))
         + ((fc('MeOH','Col4','HX18')/F('Col4','HX18'))/MW('MeOH'))
         + ((fc('Glycerol','Col4','HX18')/F('Col4','HX18'))/MW('Glycerol')));


Glydis_22..
           P_Col4_Temp_2 =E= (1+0.5*dp)*P_Col4;


Glydis_23..
          P_Col4_Temp_2 =E= (m_frac_c4('Wa','Col4','HX18')*(10**(coef_p('Wa','1')-(coef_p('Wa','2')/
         (coef_p('Wa','3')+(Col4HX18Temp+273)))))*750.06)   +
         (m_frac_c4('Glycerol','Col4','HX18')*10**(coef_p('Glycerol','1')-(coef_p('Glycerol','2')/Col4HX18Temp + 273) -
           (coef_p('Glycerol','3')*log((Col4HX18Temp+273))))) +
         (m_frac_c4('MeOH','Col4','HX18')*(10**(coef_p('MeOH','1')-(coef_p('MeOH','2')/
         (coef_p('MeOH','3')+(Col4HX18Temp+273)))))*750.06);


Glydis_24..
         T('Col4','HX18') =E=   Col4HX18Temp - 273;


**This is a constraint given in the problem **
Glydis_25..

         T('Col4','HX18') =l= 150;
*In degree C

*Heat balance to condenser and reboiler
Glydis_26..
         R_Col4 =E= 1.5*(1/(alpha_c4 - 1))*((0.9999/(m_frac_c4('MeOH','HX14','Col4')+0.001))- (alpha_c4*((1-0.999)/(1-m_frac_c4('MeOH','HX14','Col4')))));


Glydis_28..
     fc('Glycerol','Col4','HX17') =L= 0.003;



R_Col4.lo = 0.75;
R_Col4.up = 3;
*R_Col4.l = (R_Col4.up + R_Col4.lo)/2;

HX15_1..
         Q('HX15') =E= -F('Col4','HX17')*(R_Col4 + 1)*(((fc('Wa','Col4','HX17')/F('Col4','HX17'))*dH_Vap_0('Wa')) +
                 ((fc('Glycerol','Col4','HX17')/F('Col4','HX17'))*dH_Vap_0('Glycerol')) +
                 ((fc('MeOH','Col4','HX17')/F('Col4','HX17'))*dH_Vap_0('MeOH')));
HX16_1..
         Q('HX16') =E= F('Col4','HX18')*(R_Col4+1)*(((fc('Wa','Col4','HX18')/F('Col4','HX18'))*dH_Vap_0('Wa')) +
                 ((fc('Glycerol','Col4','HX18')/F('Col4','HX18'))*dH_Vap_0('Glycerol')) +
                 ((fc('MeOH','Col4','HX18')/F('Col4','HX18'))*dH_Vap_0('MeOH')));

HX17_1..
         Q('HX17') =E= (fc('Wa','Col4','HX17')*c_p_ind('Wa') + fc('Glycerol','Col4','HX17')*c_p_ind('Glycerol') + fc('MeOH','Col4','HX17')*c_p_ind('MeOH'))*(T_amb-T('Col4','HX17'));


HX18_1..
         Q('HX18') =E= (fc('Wa','Col4','HX17')*c_p_ind('Wa') + fc('Glycerol','Col4','HX17')*c_p_ind('Glycerol') + fc('MeOH','Col4','HX17')*c_p_ind('MeOH'))*(T_amb-(T('Col4','HX18')));



**Calculating heat requirements ****
Variables
         Z       Objective function variables
         QS_Max  Heat required in process;


Equations
         Heatreq1;
*, Heatreq2, Heatreq3, Heatreq4, Heatreq5, Heatreq6, Heatreq7, Heatreq8;
*Heatreq1..
*         QS_Max =G= Q('HX45');
*Heatreq2..
*         QS_Max =G= Q('HX10');
*Heatreq3..
*         QS_Max =G= Q('HX12');
*Heatreq4..
*         QS_Max =G= Q('HX14');
*Heatreq5..
*         QS_Max =G= Q('HX16');
*Heatreq6..
*         QS_Max =G= Q('HX19');
*Heatreq7..
*         QS_Max =G= Q('Reactor3');
*Heatreq8..
*         QS_Max =G= Q('HX21');

Heatreq1..
                 QS_Max =E=  Q('HX12') + Q('Reactor3') + Q('HX45') + Q('HX10') +  Q('HX16') +  Q('HX21');
*QS_Max =E=  Q('HX12') + Q('Reactor3') + Q('HX45');

**************************
*** Coal CHP unit ********
**************************
Variables
         Elec_consumed           Electricity consumed
         Steam_consumed          Steam consumed
         Fuel                     Fuel consumption based on steam
         Elec_Prod               Electricity Produced
         Elec_Sold               Electricity Sold
         Cap_CHP                 Capacity of CHP
         ksens                   fraction of heat load rejected through sensible heat transfer
         mw_evap                 Water lost through Evaporation
         mw_drift                Drift loss
         mw_bd                   Water to blowdown
         T_mw                    Total make up water required
         kbd                     Fraction of water to blowdown
         A_ct_w                  Water withdrawn per unit of energy rejected
         A_ct_c                  Water consumed per unit of energy rejected
         I_withdrawn             Water withdrawn per Kwh
         I_consumed              Water consumed per Kwh
         WaterW                  Total water withdrawn
         WaterC                  Total water consumed
         Water_makeup            Makeup water
         Total_demand            Total water demand
         CO2_D                   CO2 demand
         SO2_D                   SO2 demand
         NO2_D                   NO2 demand
         PM10_D                  PM demand
         O3_D                    Ozone demand;
Scalar
         PtoH           Power to heat ratio /0.246/
         DeltaT  Change in water temperature /25/
         Circuflow       Circulation flow of water /3.22/
         n_cc    number of cycles of concentration /15/
         dens_water      Water density /0.954/
         HR Heat rate of coal CHP /7435/
         B       Heat flow out of power plant except heat flow of cooling water /3900/
         C       Water used in other processes not related to cooling    /0.3/;
*Change this to 2991 based on PAR_sun>0 value obtained from weather station data, when running these models seasonally
*Change this to 5793 based on PAR value from weather station data
Equations
         CHP_1, CHP_2, CHP_3, CHP_4, CHP_5, CHP_6, CHP_7, CHP_8, CHP_9, CHP_10, CHP_11, CHP_13, CHP_14,CHP_15, CHP_16, CHP_17, CHP_18,
         CHP_19, CHP_24, CHP_25, CHP_26, CHP_27;

T_mw.lo = 0.00001;

CHP_1..
         Elec_consumed =E= ((fc('FAME','Col5','HX22')*24*3600))*0.0502*3600*245;
*Ref: http://www.usda.gov/oce/reports/energy/EnergyLifeCycleSoybeanBiodieseI6-11.pdf
*Units in KJ/year.Conversion from Kwh/year to Kj/year


CHP_2..
         Steam_consumed =E= (QS_Max)*5880;
*Units in Kj/year

CHP_3..
         Fuel =E= Steam_consumed*3600/(0.80*30450.38);
*Fuel requirement in Kg per year
CHP_4..
         Elec_Prod =E= (Steam_consumed*3600*PtoH)*0.60;
*Elec_Prod =E=(Fuelreq)*30450.38*0.10;
*Electricity produced, in Kj/year. Excess is sold to grid
*Electrical efficiency of 10%
CHP_5..
          Elec_Sold =E= Elec_Prod -(Elec_consumed + 0.10*Elec_Prod);
*In Kj/year
CHP_6..
         Cap_CHP =E= (Elec_Prod/5880)*2.777777777778E-7;
*Converting Kj/year to Kj/hour and Kj/hour to Mw
CHP_7..
          ksens =E= ((-0.0002793*DeltaT) + (0.000109*(DeltaT)**2) - (0.345*DeltaT) + 26.7)/100;
CHP_8..
         mw_evap =E= (0.85*(1/100)*DeltaT)/(10)*Circuflow;
*Units in gpm
CHP_9..
          mw_drift =E= 0.00200*Circuflow;
*Drift loss in gpm
CHP_10..
          mw_bd =E= mw_evap/(n_cc - 1);
*Units in gpm

CHP_11..
         T_mw =E= mw_evap + mw_bd + mw_drift;
*Units in gpm
*CHP_12..
*         kbd*T_mw =E= mw_bd;
*Units in gpm
CHP_13..
         A_ct_w =E=  ((1 - ksens)/(dens_water*2300))*(1+(1/(n_cc - 1)));
CHP_14..
         A_ct_c =E= ((1 - ksens)/(dens_water*2300))*(1+((1-(mw_bd/T_mw))/(n_cc - 1)));

*A is the water needed per unit of energy rejected. A_ct_c and A_ct_w are in units of L/Kj

CHP_15..
         I_withdrawn =E=  A_ct_w*(HR-B)  + C ;
*In L/kwh
CHP_16..
         I_consumed =E=  A_ct_c*(HR-B)  + C ;
*In L/kwh
CHP_17..
         WaterW =E= I_withdrawn*(Elec_Prod/3600);
*In liters per year or Kg/year
CHP_18..
         WaterC =E=  I_consumed*(Elec_Prod/3600);
*In kg/year or liters / year
CHP_19..
         Water_makeup =E= T_mw*0.063;
*In Kg/s

*CHP_20..
*         fc('Wa','Src6','CHP') =E= Water_makeup;
* + (WaterC/Annual);

*Total_withdrawn =E= (fc('Wa','Src6','Mix2')  + WaterW*0.00028)*Annual;
*In Kg/year
*CHP_21..
*         fc('Wa','SrcWa','Src6') =E= fc('Wa','Src6','CHP') + fc('Wa','Src6','Sep2');

*CHP_22..
*         Total_demand =E= fc('Wa','SrcWa','Src6') + (WaterW/Annual);
*In Kg/s
*CHP_23..
*         CO2_D =E= 402.989*2.77778e-7*(Elec_Prod - Elec_Sold)*453.592 + 0.00391845*2.20462*1000*fc('FAME','Col5','HX22')*Annual;
*Emissions in gms/year
*CHP emissions factor of 1178 lb/Mwh and a MEA efficiency of 65%
CHP_24..
         NO2_D =E=0.3693*2.77778e-7*(Elec_Prod - Elec_Sold)*453.592;
*0.3693
*Assuming that NO2 is converted to Ozone during summer months when PAR>0, and stays as NO2 in winter months
*NO2 emissions factor was calcualted from CHP emissions calculator with default NO2 emissions rate as 2.462 lb/Mwh and a 85% SCR efficiency = 0.3693 lb/Mwh
CHP_25..
         SO2_D =E= 0.1122*2.77778e-7*(Elec_Prod - Elec_Sold)*453.592;
*0.29256
CHP_26..
         PM10_D =E= 0.14*2.77778e-7*(Elec_Prod - Elec_Sold)*453.592;
*0.6826 lb/Mwh
*Value obtained from:  https://www.epa.gov/sites/production/files/2015-07/documents/output-based_regulations_a_handbook_for_air_regulators.pdf
*Current output based emissions factor is 0.7 lb/Mwh and with a baghouse filter of efficiency 0.80 - emissions factor is 0.14 lb/mwh
CHP_27..
         O3_D =E= 2.8381*2.77778e-7*(Elec_Prod - Elec_Sold)*453.592;
*Assuming that Ozone is formed from the conversion of NO2 during summer months (Leaf on dates from March to December, and when PAR >0)

*CHP_28..
*         fc('Wa','SrcWa','Src6') + (WaterW/Annual) =L= groundwater_supply + External_supply;



Positive Variables
         x1      Forest ecosystem
         x2      Control equipments factor;
x1.lo = 0.00000001;
x1.up = 1;
x1.l = (x1.lo + x1.up)/2;
x2.lo = 0.00000001;
x2.up = 1;
x2.l = (x2.lo + x2.up)/2;



**********************************
*** Forest Ecosystem *************
**********************************
*Tree 1: White Oak
*Tree 2: American Elm
*Tree 3: Eastern hemlock


*12 years
Table
         coef_flux(K,p)
                      1         2            3
         NO2         0.7627   0.8311     1.0297
         PM10        0.2043   0.2043     0.2043
         SO2         0.3395   0.3421     0.354
         O3          3.2846   3.341     3.4969
         CO2         65.878  23.621     24.5672
         GWater      0.0435  0.0346     0.01610;
Positive variables
         Area_oak                        Area of oak
         Area_elm                        Area of elm trees
         Area_hem                        Area of hemlock trees
         SO2_S                           SO2 Supply
         NO2_S                           NO2 Supply
         PM10_S                          PM10 Supply
         O3_S                            O3 Supply
         CO2_S                           CO2 supply
         Total_area                      Total tree area
         groundwater_supply       Ground water supply from restoration;
Variables
         net_supply              Net supply of air quality regulaiton ES;


SO2_D.lo = 0.01;
NO2_D.lo = 0.01;
PM10_D.lo = 0.01;
O3_D.lo = 0.01;
CO2_D.lo = 0.01;
Total_area.up =161874;
*Max area is 39 acres
Area_oak.up = 161874;
*53958;
*Max ins 33% of 40 acres, in m2
Area_elm.up =161874;
*Max ins 33% of 40 acres, in m2
Area_hem.up = 161874;
*Max ins 34% of 40  acres, in m2
Area_hem.lo = 0;
Area_elm.lo = 0;
Area_oak.lo = 0;
Area_hem.l = (Area_hem.lo + Area_hem.up)/2;
Area_oak.l = (Area_oak.lo + Area_oak.up)/2;
Area_elm.l = (Area_elm.lo + Area_elm.up)/2;




Equations
        Tree5, Tree13, Tree14, Tree15, Tree16, Tree21, Tree22;
*, Tree17, Tree18, Tree19, Tree20,Tree21, Tree22;

*Tree1..
*         PM10_D =E= Area_oak*coef_flux('PM10','2') + Area_elm*coef_flux('PM10','1') + Area_hem*coef_flux('PM10','3');
Tree5..
         Total_area =E= (Area_oak + Area_elm + Area_hem)*x1;
Tree13..
         SO2_S =E= x1*((Area_elm*coef_flux('SO2','2') + Area_oak*coef_flux('SO2','1') + Area_hem*coef_flux('SO2','3')));
Tree14..
         NO2_S =E= x1*((Area_elm*coef_flux('NO2','2') + Area_oak*coef_flux('NO2','1') + Area_hem*coef_flux('NO2','3')));
Tree15..
         PM10_S =E= x1*((Area_elm*coef_flux('PM10','2') + Area_oak*coef_flux('PM10','1') + Area_hem*coef_flux('PM10','3')));
Tree16..
         O3_S =E= x1*((Area_elm*coef_flux('O3','2') + Area_oak*coef_flux('O3','1') + Area_hem*coef_flux('O3','3')));
Tree21..
        CO2_S =E=x1*((Area_elm*coef_flux('CO2','2') + Area_oak*coef_flux('CO2','1') + Area_hem*coef_flux('CO2','3')));
Tree22..
         groundwater_supply =E= x1*((Area_elm*coef_flux('GWater','2') + Area_oak*coef_flux('GWater','1') + Area_hem*coef_flux('GWater','3'))* 0.00003);
*Conversion from m3/year to Kg/s


*groundwater_supply =E= 0.0077*Total_area*0.012*0.0625;
*Value in (m/day)*m2 converted to Kg/s

*Assume that tree intercepts 25% of total rainfall, and 25% of intercepted rainfall reaches the ground and goes as infiltrated water.
*Tree23..
*       net_supply  =e= (PM10_S - PM10_D);
*This is to match the conventional case design


***************************
*** Wastewater streams ***
***************************
Variables
         out_MeOH                Methanol wastewater
         out_Wa                  Outlet wastewater stream
         out_FFA                 Outlet FFA stream
         out_oil                 Outlet oil stream
         out_wa_flow             Wastewater flowrate
         out_glycerol            Outlet of glycerol;


out_wa_flow.lo = 0.001;
out_Wa.lo = 0.0001;
mw_bd.lo = 0.0001;
Equations
         WW1, WW2, WW4, WW5, WW6;

WW1..
         out_MeOH =E= (fc('MeOH','Col4','HX17'))*24*3600;
*In Kg/day

WW2..
         out_Wa =E=  (fc('Wa','Col4','HX17'))*24*3600;
*In Kg/day
*Removed fc('Wa','Col5','HX23') stream


*WW3..
*         out_FFA =E= fc('FFA','HX14','Col4')*24*3600;
*In Kg/day

WW4..
         out_oil =E= (fc('Oil','Reactor4','HX14'))*24*3600;
*In Kg/day


WW6..
         out_glycerol =E= fc('Glycerol','Col4','HX17')*24*3600;
*In kg/day
WW5..
         out_wa_flow =E= (out_MeOH/792) +  (out_Wa/1000) + (out_oil/926.1) + mw_bd*5.5555 + (out_glycerol/1126);
*This is in units of m3/day





Positive Variables
        vol_basin_floc   Flowrate in m3 per hour and detention time in hour. Vol of basin converted to ft3
        cc_coag          Capital cost of coagulation unit
        cc_floc          Capital cost of flocculation unit
        oc_CF            Operation cost of coag floc unit
        cc_CF            Capital cost of coag floc unit
        Amount_alum      Amount of alum needed for removing oil and grease
        vol_basin_coag   Volume of coagulaiton unit
        coag_electricity electricity consumption in a coagulaiton unit
        floc_electricity electricity consumption in flocculation unit;


**** Coagulation flocculation Unit *****
***************************************


Scalar
         CF_coag          Weight of coagulant required in kg per gal of water treated  /0.000094625/
         dose             Dose of alum for removing oil and grease /100/;
*Units in Kg/l

Equations
         CF_1, CF_2, CF_3, CF_4, CF_5, CF_7, CF_8, CF_9, CF_6;

CF_6..
         Amount_alum =E= out_wa_flow*dose/1e3;
*Units in Kg/day Ref: http://www.ecs.umass.edu/cee/reckhow/courses/371/371hw06/371hw06s.pdf

CF_7..
         vol_basin_coag =E= out_wa_flow*(1/1440);
*Detention time is 1 minute, converted to days. Units in m3

CF_8..
         coag_electricity =E= 1.307e-3*vol_basin_coag*(600*600)*3600;
*Units in J/hour

CF_1..
         vol_basin_floc =E= (out_wa_flow/24)*(45/60)*35.3147;
*This is the vol for one basin. Assuming there are 3 totally.Flowrate in m3/day and detention time in min converted to days. Vol of basin converted in m3

CF_9..
         floc_electricity =E= 1.307e-3*(3*vol_basin_floc/35.3147)*3600*(50*50 + 20*20 + 10*10);
*Calculated as the power required in each basin per hour. G1 = 50/sec, G2 = 20/sec, G3 = 10/sec
*1.3e-3 is the absolute viscosity of water
*Units in J/hour

CF_2..
        cc_coag =E= 240.78*out_wa_flow*91.86 + 71071;
*cc_coag =E= 148.5*vol_basin_coag*35.3147;
*cc_coag =E= 240.78*out_wa_flow*91.86 + 71071;
*;
*$148.5 per ft3. Tank volume from m3 to ft3
*
*outlet water flow converted from m3/day to lbs/hour (91.86) and 240.78 is the capital cost equation

CF_3..
         cc_floc =E= 5*10**(-12)*((vol_basin_floc)**3)- 7*10**(-6)*((vol_basin_floc)**2) + 7.7415*(vol_basin_floc) + 158139;
*Ref: https://uta-ir.tdl.org/uta-ir/bitstream/handle/10106/4924/Sharma_uta_2502M_10652.pdf?sequence=1&isAllowed=y
*G = 50/sec
CF_4..
         oc_CF =E= (3E-13*(vol_basin_floc*35.3147)**3)  - (4E-17*(vol_basin_floc*35.3147)**2) + (0.318*(vol_basin_floc*35.3147)) + 6040
         + 0.425*Amount_alum*245 - (3E-8*(vol_basin_coag*35.3147)**3) + (0.0008*(vol_basin_coag)**2) + (7.8303*vol_basin_coag*35.3147) +
         22588;
*oc_CF =E= out_wa_flow*(264.17/1000)*25000 + (CF_coag*out_wa_flow*264.17*365*0.4);
*

*Conversion of water flowrateto gals per day assuming each Floc tank can hold 1000 gals.
*Total operating cost of CF unit
CF_5..
         cc_CF =E= cc_coag + cc_floc;
*Total capital cost of CF unit




********************
** Wetland Design **
********************

Positive variables
         C_in            Concentration of water stream from CF unit
         Cin             Concentration of water stream to wetland
         KV              Rate constant at a diff temp in per day
         alpha           Precip  - Evapotrans
         Cprime          Background conc for KV
         Power           Power term
         Wet_Area        Area of wetland in m2
         Wet_Temp        Temporary - concentration calculation
         LHSdr           Wetland denominator
         LHSnr           Wetland numerator
         LHSW             Left hand side term
         RHS             Right hand side
         cc_wetland      Capital cost of wetland
         oc_wetland      Operating cost of wetland
         External_supply External water withdrawal
         Supply          Wetland water supply
         Delta           Net water supply (S-D);

Positive variables
           HRT             Hydraulic retention time;

Scalar
         n               Average bed porosity /0.42/
         d               Bed depth in meters     /0.6/
         Kv_20           Volumteric rate constant at 20 C  in per day /0.031/
         theta           Tempeature factor /1.06/
         Precip          Precipitation in m per day /0.0077/
         Evaptran        Evapotranspiration in m per day /0.00308/
         Cb              Background COD concentration /30/
         Ce              Effluent concentration /75/;
HRT.lo = 0.001;
Kv.lo = 0.001;
Power.lo = 0.001;
alpha.lo = 0.001;
Wet_Temp.lo = 0.001;
LHSdr.lo = 0.001;
LHSnr.lo = 0.001;
LHSW.lo = 0.001;
HRT.up = 5;
HRT.l = (HRT.up + HRT.lo)/2;


Equations
         Wetland1, Wetland2, Wetland3,  Wetland4,  Wetland5, Wetland6,  Wetland7,  Wetland8,  Wetland9,  Wetland10, Wetland12, Wetland13, Wetland14,
         Wetland15, Wetland16, Wetland18, Wetland19, Wetland20, Wetland21, Wetland17, Wetland22;

Wetland1..
       C_in =E= (1.5*out_MeOH + 1.26*out_glycerol)*1000/((out_Wa/1000) +mw_bd*5.5555);
*COD: Units in mg/liter of water
Wetland2..
         Cin =E= C_in*0.47;
*Coagulation Flocculation unit removed 53% of COD
Wetland3..
         alpha =E= Precip - Evaptran;
*Units in m/day
Wetland4..
          KV =E= Kv_20*theta**(Temp_cooldown - 20);
Wetland5..
         Cprime*((KV*n*D + alpha)) =E= Cb*((KV*n*d));
*Units in mg/liter

Wetland6..
         LHSdr =E= (Cin - Cprime);
Wetland7..
         LHSnr =E= (Ce - Cprime);
Wetland8..
         LHSW*LHSdr =E=  LHSnr;
Wetland9..
         Power =E= (1 + (KV*n*d*216.45));
Wetland10..
         (1 + (alpha/HRT))**(-Power) =E= LHSW;
Wetland12..
         Wet_Area*HRT =E= out_wa_flow;
*Wetland15..
*        Wet_Area =L= 60702.8;

Wetland13..
         cc_wetland =E= 20.6234*Wet_Area;
*Capital cost of $20.4 per m2
Wetland14..
         oc_wetland =E= 1.44*Wet_Area;
*O&M cost of $1.44 per m2
*Water to CHP from wetland in Kg/s.
Wetland15..
         fc('Wa','Src6','CHP') =E= Water_makeup  ;
*In Kg/s
Wetland16..
         fc('Wa','Wetland','Src6') =G= 0.0001*Wet_Area*(d-0.2)*0.000032;
*m3/year converted to Kg/s of water. Water withdrawal from wetland can be only 50% of the max water available.

Wetland17..
         fc('Wa','Wetland','Src6') =L= 0.50*Wet_Area*(d-0.2)*0.000032;
*fc('Wa','SrcWa','Src6') =E= (1-x2)*(out_wa_flow*0.012);
Wetland18..
*Total_demand =E= External_supply + Supply;
              fc('Wa','SrcWa','Src6') + fc('Wa','Wetland','Src6') =E= fc('Wa','Src6','CHP') + fc('Wa','Src6','Sep2') + (WaterW/Annual);

Wetland19..
         Total_demand =E= fc('Wa','Src6','CHP') + fc('Wa','Src6','Sep2') + (WaterW/Annual);
*In Kg/s

Wetland20..
         External_supply =E= fc('Wa','SrcWa','Src6');
*In Kg/s
Wetland21..
         Supply =E= fc('Wa','Wetland','Src6') + groundwater_supply;
Wetland22..
         Delta =E= Total_demand - fc('Wa','Wetland','Src6');




*** Technological systems design ************
********************************************
****** New FGD Design ***********
******************************
Variables
         uncontrolled_so2                Uncontrolled SO2 emissions in lb per MMBtu
         controlled_so2                  Target emissions rate in lb per MMbtu
         SO2_removal                     SO2 removal efficiency
         Scaling_factor_FGD              FGD scaling factor
         TCI_FGD                         Capital cost FGD
         Fixed_cost_FGD                  Fixed cost FGD
         Variable_cost_FGD               Variable cost FGD
         SO2_removed                     Amount of SO2 removed
         FGD_efficiency                  FGD efficiency;
uncontrolled_so2.lo = 0.0001;
SO2_removal.lo = 0.00000001;
Cap_CHP.lo = 0.00001;
Elec_consumed.lo = 0.000001;
SO2_removed.lo = 0.0001;



Equations
         FGD1, FGD2, FGD3,FGD4,FGD5,FGD6,FGD7, FGD8;



FGD1..
*uncontrolled_so2*Elec_consumed =E= SO2_D*0.000108;
                                 uncontrolled_so2 =E= ((((SO2_D*1e-3)/(Elec_consumed*1e-3/3600))*0.60)/3.413)*2.204;

FGD2..
*controlled_so2*Elec_consumed =E= SO2_removed*0.000108;
                                 controlled_so2 =E=((((SO2_removed)*1e-3)/(Elec_consumed*1e-3/3600))*0.60/3.413)*2.204;

FGD3..                           SO2_removal =E= (controlled_so2)/uncontrolled_so2;

FGD4..
                                   Scaling_factor_FGD =E= (500/Cap_CHP)**0.6;
FGD8..
                                  FGD_efficiency=e=(SO2_removal)*x2;

FGD5..
         TCI_FGD*0.9 =E= 149*Scaling_factor_FGD*Cap_CHP*1000*1.76*FGD_efficiency;
*1.76 is the cost inflation value FROM 1990
FGD6..
         Fixed_cost_FGD =E= 5.4*Cap_CHP*1000*1.76;
FGD7..
         Variable_cost_FGD =E= 0.83*Cap_CHP*0.65*5880;
*FGD9..
*         SO2_removed =E= SO2_D;




*** SCR Design ***
*******************
Variables
         uncontrolled_nox        Uncontrolled Nox emissions
         controlled_nox          Controlled NOx emissions
         NO2_Seq                 NO2 sequestration
         NO2_D                   NO2 demand
         Scaling_factor          Scaling factor
         TCI_SCR                 Capital cost SCR
         Fixed_cost_SCR           Fixed cost SCR
         Variable_cost_SCR       Variable cost SCR
         NOx_removal             NOx removal efficiency
         SCR_efficiency          SCR efficiency
         Net_removal_NO2         Net removal of NO2;

Elec_consumed.lo = 0.0001;
Cap_CHP.lo = 0.0001;

uncontrolled_nox.lo = 0.0001;
NOx_removal.lo = 0.00001;
NO2_Seq.lo = 0.0001;

Equations
         SCR1, SCR2, SCR3, SCR4, SCR5, SCR6, SCR7, SCR8;


SCR1..
         uncontrolled_nox =E= ((((NO2_D*1e-3)/(Elec_consumed*1e-3/3600))*0.60)/3.413)*2.204;
*Units in Kg/MMBtu. Conversion of NO2 emissions to Kg, elec consumed from KJ to Kwh to Mwh and Kg/Mwh to Kg/mmBtu to lb/MMbtu

SCR2..
         controlled_nox =E=((((NO2_Seq)*1e-3)/(Elec_consumed*1e-3/3600))*0.60/3.413)*2.204;
*Target NOx emissions which will vary depending on supply

SCR3..
         NOx_removal =E= (controlled_nox)*x2/uncontrolled_nox;

SCR4..
         Scaling_factor =E= (243/Cap_CHP)**0.27;

SCR5..
         (TCI_SCR*0.9) =E= 100*Scaling_factor*Cap_CHP*1000*1.38*NOx_removal;
*1.38 is the cost inflation value
SCR6..
         Fixed_cost_SCR =E= 0.066*Cap_CHP*1000*1.38*0.9;
SCR7..
         Variable_cost_SCR =E= 0.6*Cap_CHP*0.65*5880*0.9;
SCR8..
         SCR_efficiency =E= (NOx_removal);
*SCR9..
*         NO2_Seq =E= NO2_D;


**** Baghouse design ***
**********************
Variables
         heat_input              System heat
         stack_flow              Flow rate of flue gas in stack
         Capital_cost_baghouse Capital cost of baghouse
         OM_cost_baghouse        OM cost
         PM10_seq                PM10 sequestered
         PM10_removal            Efficiency
         Fuel_consumption        Fuel consumption
         uncontrolled_pm10       Uncontrolled pm10 emissions
         controlled_pm10         PM10 removal
         Net_removal_PM10        Net PM10 removal;


uncontrolled_pm10.lo = 0.0001;
PM10_removal.lo = 0.00001;
PM10_seq.lo=0.0001;
Equations
         Baghouse1, Baghouse2, Baghouse3, Baghouse4, Baghouse5, Baghouse6, Baghouse7, Baghouse8, Baghouse9;

Baghouse1..
         Fuel_consumption =E=(Fuel/5880)*2.204;

Baghouse2..
          heat_input=E= (Fuel_consumption*0.453*30450.8)*9.4781e-7;
*Conversion from Kj/hour to MMBtu/hour
Baghouse3..
         stack_flow =e= 126.966*heat_input;
*Stack flow in dscf/hour
Baghouse4..
         Capital_cost_baghouse*0.9 =e= 29*stack_flow*0.016*1.41*PM10_removal;
*Conversion to cubic ft per min, 1.41 is the conversion value to 2012
Baghouse5..
         OM_cost_baghouse =e= 11*stack_flow*0.017*1.41;
Baghouse6..
         uncontrolled_pm10 =E= ((((PM10_D*1e-3)/(Elec_consumed*1e-3/3600))*0.60)/3.413)*2.204;
*Units in Kg/MMBtu. Conversion of NO2 emissions to Kg, elec consumed from KJ to Kwh to Mwh and Kg/Mwh to Kg/mmBtu to lb/MMbtu

Baghouse7..
         controlled_pm10 =E=((((PM10_Seq)*1e-3)/(Elec_consumed*1e-3/3600))*0.60/3.413)*2.204;
Baghouse8..
         PM10_removal =e= controlled_pm10*x2/uncontrolled_pm10;
Baghouse9..
         PM10_Seq + PM10_S =E= PM10_D;





** Objective function calculations **
Positive Variable
         Manuf_cost      Total manufacturing cost
         Revenue         Revenue from sale of products
         Cash_flow       Total cash flow
         CI             Total capital investment
         PV              Present value
         NPV             Net present value
         Dep             Depreciation
         Cw              Working capital;

Scalar
         FCI             Fixed Capital Investment /3219315/
         tax             tax rate                /0.05/;





Cash_flow.lo =1E-6;
*PV.up = 1E10;
*NPV.lo = 850000;


Equations
         Cost1, Cost2, Cost3, Cost4, Cost5, Cost6, Cost7;
Cost1..
         Manuf_cost =E= (347608)*0.75   + Cost_KOH*fc('KOH','Src5','Mix3')*Annual + Cost_MeOH*fc('MeOH','Src3','Mix2')*Annual + Cost_PhosA*fc('PhosA','Src7','Reactor4')*Annual
                          + Cost_Wa*Delta*3600*24*245  + (Fuel)*0.6916 + 10*Cap_CHP*1000 + Cost_Oil*fc('Oil','Src','Mix3')*Annual + oc_CF*0.75 + oc_Wetland*0.75
         + x1*(26.8*(Area_oak)*0.000247105*0.75 + 25.6*(Area_elm)*0.000247105*0.75 + 30.08*(Area_hem)*0.000247105*0.75) + (x2)*0.75*(Fixed_cost_FGD + Variable_cost_FGD + Fixed_cost_SCR + Variable_cost_SCR +OM_cost_baghouse);
*Here WaterC is the water consumed in the coal CHP in Kg/hour
*+
Cost2..
         Revenue =E= Cost_FAME*fc('FAME','Col5','HX22')*Annual + Cost_Glycerol*fc('Glycerol','Col4','HX18')*Annual +
                 Cost_PotPhos*fc('PotPhos','Reactor4','Snk3')*Annual + (Elec_Sold/3600)*0.0667;
*Electricity sold in one year
Cost3..
         Cash_flow =E= Revenue - Manuf_cost;

Cost4..
         CI =E= 0.75*(FCI  +  1536*1000*Cap_CHP + Total_area*0.000247105*Cost_land  + x1*(57.6*(Area_oak)*0.000247105 + 56.32*(Area_elm)*0.000247105 + 67.84*(Area_hem)*0.000247105) + cc_CF + cc_Wetland + TCI_FGD + TCI_SCR + Capital_cost_baghouse);
Cost5..
         Cw =e= (0.15*(FCI + 1536*1000*Cap_CHP))*0.75;
*Working capital for plant

Cost6..
         Dep =e= 0.75*((FCI  +  1536*1000*Cap_CHP - 39600)/20);
*Depreciation calculated based on straight line depreciation method. Salvage cost is assumed to be zero. Since land value does not depreciate, depreciation should be calculated without including land cost in fixed capital investment costs
*Depreciation is calculated for all "Technological" equipments


Cost7..
         NPV =E= 0.75*((-(CI + Cw) + (((Revenue - Manuf_cost)*(1-tax))*(1 - (1.07)**(-20)))/0.07 + (Dep*tax)*(1 - (1.07)**(-20))/0.07 + Cw/(1.07)**20)/20);



*Cost4..
*         TCI =E= 3219315  +  1536*1000*Cap_CHP + cc_CF + cc_Wetland + Cost_tree*(Total_area)*0.000247105  + Total_area*0.000247105*Cost_land;
*This is the water withdrawn in one day; Capital cost for CHP in Mw using $60.4/Mwh Ref: https://www.eia.gov/forecasts/aeo/electricity_generation.cfm
*Capital cost of CHP is 2595 per Kwh
*Cost5..
*         PV =E= ((Cash_flow)/0.07)*(1 - (1/((1.07)**20)));
*Cost6..

*         NPV =E= (PV - TCI)/20;
*Cost7..
*         NPV =G= 3.7E7;




Equation
         Obj;

*Obj..         (Cost_FAME*fc('FAME','Col5','HX22') + Cost_Glycerol*fc('Glycerol','Col4','HX18') + Cost_PotPhos*fc('PotPhos','Reactor4','Snk3')
*                        - Cost_KOH*fc('KOH','Src5','Mix3') - Cost_MeOH*fc('MeOH','Src3','Mix2')
*                         - Cost_PhosA*fc('PhosA','Src7','Reactor4') - (Cost_Steam/dH_Vap_0('Wa'))*QS_Max) =E= Z;

Equations
Obj;


Obj..
*Z=E= (Supply - Total_demand);
*Z=E= (groundwater_supply - Total_demand)
*Z=E= (PM10_S - PM10_D)
                                                                                                                 Z=E=NPV;
Model Biodiesel   /ALL/;

Option NLP =BARON;


*Biodiesel.optfile=1;
Solve Biodiesel Using NLP Maximizing Z;
$ontext

set
                 count   counter
                 /1*18/;

parameter lim(count)
/1                -29721.245
2                -29471.245
3                -29221.245
4                -28971.245
5                -28721.245
6                -28471.245
7                -28221.245
8                -27971.245
9                -27721.245
10                -27471.245
11                -27221.245
12                -26971.245
13                -26721.245
14                -26471.245
15                -26221.245
16                -25971.245
17                -25721.245
18                -25471.245/;


parameter lim(count)
/1        -61282.349
2        -59782.349
3        -58782.349
4        -57782.349
5        -56782.349
6        -55782.349
7        -54782.349
8        -53782.349
9        -52782.349
10        -51782.349
11        -50782.349
12        -49782.349
13        -48782.349
14        -47782.349
15        -46782.349
16        -45782.349
17        -44782.349
18        -43782.349
19        -42782.349
20        -41782.349
21        -40782.349
22        -39782.349
23        -38782.349
24        -37782.349
25        -36782.349
26        -35782.349
27        -34782.349
28        -33782.349
29        -32782.349
30        -31782.349
31        -30782.349
32        -29782.349
33        -28782.349
34        -27782.349
35        -26782.349
36        -25782.349
37        -24782.349
38        -23782.349
39        -22782.349
40        -21782.349
41        -20782.349 /;



parameter data(count);


Loop (count,
net_supply.fx = lim(count);
Solve Biodiesel Using NLP Maximizing Z;
data(count) = Obj.l;
);

display data;
$offtext
