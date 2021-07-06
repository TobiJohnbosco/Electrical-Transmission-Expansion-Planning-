$title Transmission Expansion Planning

Set
   bus        / 1*14   /
   slack(bus) / 1     /
   Gen        / g1*g8 /
   k          / k1*k4 /;

Scalar
   Sbase /  100 /
   M     / 1000 /;

Alias (bus,node);

Table GenData(Gen,*) 'generating units characteristics'
       b   pmin  pmax
   g1  20  0     400
   g2  30  0     400
   g3  10  0     600
   g4  20  0     400
   g5  10  0     400
* -----------------------------------------------------

Set GBconect(bus,Gen) 'connectivity index of each generating unit to each bus' / 1.g1, 2.g2, 3.g3, 6.g4, 8.g5 /;

Table BusData(bus,*) 'demands of each bus in MW'
      Pd
   1  80
   2  240
   3  40
   4  160
   5  240
   6  360
   7  100
   8  200
   9  220
   10 360
   11 80
   12 40
   13 160
   14 180;

Table branch(bus,node,*) 'network technical characteristics'
            X    LIMIT  Cost  stat
   1.2      0.059  100    40    1
   1.5      0.223  80     60    1
   2.3      0.197  100    20    1
   2.4      0.176  100    20    1
   2.5      0.174  100    40    1
   3.4      0.171  100    30    1
   4.5      0.042  100    20    1
   4.7      0.209  100    30    0
   4.9      0.556  100    30    0
   5.6      0.250  100    20    0
   6.11     0.199  100    40    1
   6.12     0.256  80     20    1
   6.13     0.130  100    30    1
   7.8      0.176  100    40    0
   7.9      0.110  100    30    0
   9.10     0.085  100    20    1
   9.14     0.270  100    30    1
   10.11    0.192  100    20    1
   12.13    0.199  100    30    1
   13.14    0.348  100    40    1;
Set conex(bus,node) 'Bus connectivity matrixl';
conex(bus,node)$(branch(bus,node,'x')) = yes;
conex(bus,node)$conex(node,bus)        = yes;

branch(bus,node,'x')$branch(node,bus,'x')             =  branch(node,bus,'x');
branch(bus,node,'cost')$branch(node,bus,'cost')       =  branch(node,bus,'cost');
branch(bus,node,'stat')$branch(node,bus,'stat')       =  branch(node,bus,'stat');
branch(bus,node,'Limit')$(branch(bus,node,'Limit')=0) =  branch(node,bus,'Limit');
branch(bus,node,'bij')$conex(bus,node)                =1/branch(bus,node,'x');
M = smax((bus,node)$conex(bus,node),branch(bus,node,'bij')*3.14*2);
*****************************************************

Variable OF, Pij(bus,node,k), Pg(Gen), delta(bus), LS(bus);
Binary Variable alpha(bus,node,k);
alpha.l(bus,node,k) = 1;
alpha.fx(bus,node,k)$(conex(bus,node) and ord(k)=1 and branch(node,bus,'stat')) = 1;

Equation const1A, const1B, const1C, const1D, const1E, const2, const3;
***********************************************************************

const1A(bus,node,k)$conex(node,bus)..
   Pij(bus,node,k) - branch(bus,node,'bij')*(delta(bus) - delta(node)) =l=  M*(1 - alpha(bus,node,k));

const1B(bus,node,k)$conex(node,bus)..
   Pij(bus,node,k) - branch(bus,node,'bij')*(delta(bus) - delta(node)) =g= -M*(1 - alpha(bus,node,k));

const1C(bus,node,k)$conex(node,bus)..
   Pij(bus,node,k)   =l= alpha(bus,node,k)*branch(bus,node,'Limit')/Sbase;

const1D(bus,node,k)$conex(node,bus)..
   Pij(bus,node,k)   =g=-alpha(bus,node,k)*branch(bus,node,'Limit')/Sbase;

const1E(bus,node,k)$conex(node,bus)..
   alpha(bus,node,k) =e= alpha(node,bus,k);

const2(bus)..
   LS(bus) + sum(Gen$GBconect(bus,Gen), Pg(Gen))-BusData(bus,'pd')/Sbase =e= sum((k,node)$conex(node,bus), Pij(bus,node,k));

const3..
   OF =g= 10*8760*(sum(Gen, Pg(Gen)*GenData(Gen,'b')*Sbase)+100000*sum(bus ,LS(bus)))
       +  1e6*sum((bus,node,k)$conex(node,bus), 0.5*branch(bus,node,'cost')*alpha(bus,node,k)$(ord(k)>1 or branch(node,bus,'stat')=0));

Model loadflow / all /;

LS.up(bus) = BusData(bus,'pd')/Sbase;
LS.lo(bus) = 0;
Pg.lo(Gen) = GenData(Gen,'Pmin')/Sbase;
Pg.up(Gen) = GenData(Gen,'Pmax')/Sbase;

delta.up(bus)   = pi/3;
delta.lo(bus)   =-pi/3;
delta.fx(slack) = 0;
Pij.up(bus,node,k)$((conex(bus,node))) = 1*branch(bus,node,'Limit')/Sbase;
Pij.lo(bus,node,k)$((conex(bus,node))) =-1*branch(bus,node,'Limit')/Sbase;

option optCr = 0;
solve loadflow minimizing OF using mip;