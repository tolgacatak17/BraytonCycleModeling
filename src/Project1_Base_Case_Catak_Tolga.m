clc 
clear all

Tref1 = 0; % in kelvin
Pref = 101.325; % in kilo pascal

Tref2 = 300; %in kelvin

nsc = 0.8;
nst = 0.9;

alpha = 3.653;
beta = -1.337*10^-3;
gamma = 3.294*10^-6;
lambda = -1.913*10^-9;
epsilon = 0.2763*10^-12;

Rprime = 8.314; %in kj/kmole*K
R = Rprime/28.97;

%cp = (alpha + beta*T + gamma*(T^2) + lambda*(T^3) + epsilon*(T^4))*R;

%% STATE 1
T1 = 300; %in Kelvin
P1 = 100; %in kilo pascal

cp1_prime = zeros(1,300);
cp1_prime (1,1) = alpha;

cp1 = zeros(1,300);
cp1(1,1) = cp1_prime (1,1)/28.97;

h1 = 0;
for i= 1 : 1 : 299
    cp1_prime (1,i+1) = (alpha + beta*(0+i) + gamma*((0+i)^2) + lambda*((0+i)^3) + epsilon*((0+i)^4))*Rprime;
    cp1(1,i+1) = cp1_prime(1,i+1)/28.97;

    temp_area = (cp1(1,i)+cp1(1,i+1))/2;
    h1 = h1 + temp_area; 
end

s1 = - R*log(P1/Pref) + 1.70203;

k1= cp1(1,300)/(cp1(1,300)-R);
%% STATE 2
P2 = 4*P1;
s2s = s1;

cp2s_prime = zeros(1,3000001);
cp2s_prime (1,1) = alpha;

cp2s = zeros(1,3000001);
cp2s(1,1) = cp2s_prime (1,1)/28.97;

temp_area = 0;

s2s_int = 1.70203-R*log(P2/Pref);
s2s_exp = s2s_int;
for i= 3000000 : 1 : 6000000
    cp2s_prime (1,i+1) = (alpha + beta*(i/10000) + gamma*((i/10000)^2) + lambda*((i/10000)^3) + epsilon*((i/10000)^4))*Rprime;
    cp2s(1,i+1) = cp2s_prime(1,i+1)/28.97;

    temp_area = (((cp2s(1,i)/(i/10000))+((cp2s(1,i+1))/(i/10000)))*0.0001)/2;
    s2s_exp = s2s_exp + temp_area;
    delta = s2s_exp - s2s ;

    if (delta > 0) 
        T2s = i/10000;
        break
    end
end

temp_area = 0;
h2s = 0;
for i= 1 : 1 : 4448397
    cp2s_prime (1,i+1) = (alpha + beta*(i/10000) + gamma*((i/10000)^2) + lambda*((i/10000)^3) + epsilon*((i/10000)^4))*Rprime;
    cp2s(1,i+1) = cp2s_prime(1,i+1)/28.97;
    temp_area = ((cp2s(1,i)+cp2s(1,i+1))*0.0001)/2;
    h2s = h2s + temp_area; 
end

h2 = h1 + (h2s-h1)/nsc; 

cp2_prime = zeros(1,6000001);
cp2_prime (1,1) = alpha;

cp2 = zeros(1,6000001);
cp2(1,1) = cp2_prime (1,1)/28.97;

h2_exp = 0;

for i= 1 : 1 : 6000000
    cp2_prime (1,i+1) = (alpha + beta*(i/10000) + gamma*((i/10000)^2) + lambda*((i/10000)^3) + epsilon*((i/10000)^4))*Rprime;
    cp2(1,i+1) = cp2_prime(1,i+1)/28.97;

    temp_area = ((cp2(1,i)+cp2(1,i+1))*0.0001)/2;
    h2_exp = h2_exp + temp_area;
    delta = h2_exp - h2 ;

    if (delta > 0) 
        T2 = i/10000;
        break
    end
end

temp_area = 0;
s2_int = 1.70203;
s2 = s2_int;
for i= 3000000 : 1 : 4806172 
    cp2_prime (1,i+1) = (alpha/(i/10000) + beta + gamma*((i/10000)^1) + lambda*((i/10000)^2) + epsilon*((i/10000)^3))*Rprime;
    cp2(1,i+1) = cp2_prime(1,i+1)/28.97;
    temp_area = (((cp2(1,i)+cp2(1,i+1))*0.0001)/2);
    s2 = s2 + temp_area;
end
s2 = s2 - R*log(P2/Pref); 

k2=cp2(1,T2*10000)/(cp2(1,T2*10000)-R); 
%% STATE 3

T3 = T1; %in Kelvin
P3 = P2; %in kilo pascal

cp3_prime = zeros(1,300);
cp3_prime (1,1) = alpha;

cp3 = zeros(1,300);
cp3(1,1) = cp3_prime (1,1)/28.97;

h3 = 0;
for i= 1 : 1 : 299
    cp3_prime (1,i+1) = (alpha + beta*(0+i) + gamma*((0+i)^2) + lambda*((0+i)^3) + epsilon*((0+i)^4))*Rprime;
    cp3(1,i+1) = cp3_prime(1,i+1)/28.97;

    temp_area = (cp3(1,i)+cp3(1,i+1))/2;
    h3 = h3 + temp_area; 
end

s3 = - R*log(P3/Pref) + 1.70203;

k3= cp3(1,300)/(cp3(1,300)-R);
V_dot3=0; %hesaplayip duzelt
%% STATE 4

P4 = 4*P3;
s4s = s3;

cp4s_prime = zeros(1,3000001);
cp4s_prime (1,1) = alpha;

cp4s = zeros(1,3000001);
cp4s(1,1) = cp4s_prime (1,1)/28.97;

temp_area = 0;

s4s_int = 1.70203-R*log(P4/Pref);
s4s_exp = s4s_int;
for i= 3000000 : 1 : 6000000
    cp4s_prime (1,i+1) = (alpha + beta*(i/10000) + gamma*((i/10000)^2) + lambda*((i/10000)^3) + epsilon*((i/10000)^4))*Rprime;
    cp4s(1,i+1) = cp4s_prime(1,i+1)/28.97;

    temp_area = (((cp4s(1,i)/(i/10000))+((cp4s(1,i+1))/(i/10000)))*0.0001)/2;
    s4s_exp = s4s_exp + temp_area;
    delta = s4s_exp - s4s ;

    if (delta > 0) 
        T4s = i/10000;
        break
    end
end

temp_area = 0;
h4s = 0;

for i= 1 : 1 : T4s*10000
    cp4s_prime (1,i+1) = (alpha + beta*(i/10000) + gamma*((i/10000)^2) + lambda*((i/10000)^3) + epsilon*((i/10000)^4))*Rprime;
    cp4s(1,i+1) = cp4s_prime(1,i+1)/28.97;
    temp_area = ((cp4s(1,i)+cp4s(1,i+1))*0.0001)/2;
    h4s = h4s + temp_area; 
end

h4 = h3 + (h4s-h3)/nsc; 

cp4_prime = zeros(1,6000001);
cp4_prime (1,1) = alpha;

cp4 = zeros(1,6000001);
cp4(1,1) = cp4_prime (1,1)/28.97;

h4_exp = 0;

for i= 1 : 1 : 6000000
    cp4_prime (1,i+1) = (alpha + beta*(i/10000) + gamma*((i/10000)^2) + lambda*((i/10000)^3) + epsilon*((i/10000)^4))*Rprime;
    cp4(1,i+1) = cp4_prime(1,i+1)/28.97;

    temp_area = ((cp4(1,i)+cp4(1,i+1))*0.0001)/2;
    h4_exp = h4_exp + temp_area;
    delta = h4_exp - h4 ;

    if (delta > 0) 
        T4 = i/10000;
        break
    end
end

k4=cp4(1,T4*10000)/(cp4(1,T4*10000)-R);

temp_area = 0;
s4_int = 1.70203;
s4 = s4_int;
for i= 3000000 : 1 : T4*10000 
    cp4_prime (1,i+1) = (alpha/(i/10000) + beta + gamma*((i/10000)^1) + lambda*((i/10000)^2) + epsilon*((i/10000)^3))*Rprime;
    cp4(1,i+1) = cp4_prime(1,i+1)/28.97;
    temp_area = (((cp4(1,i)+cp4(1,i+1))*0.0001)/2);
    s4 = s4 + temp_area;
end
s4 = s4 - R*log(P4/Pref);

%% STATE 6

T6 = 1000; %in Kelvin
P6 = P4*0.95; %in kilo pascal

cp6_prime = zeros(1,1000);
cp6_prime (1,1) = alpha;

cp6 = zeros(1,1000);
cp6(1,1) = cp6_prime (1,1)/28.97;

h6 = 0;
for i= 1 : 1 : 999
    cp6_prime (1,i+1) = (alpha + beta*(0+i) + gamma*((0+i)^2) + lambda*((0+i)^3) + epsilon*((0+i)^4))*Rprime;
    cp6(1,i+1) = cp6_prime(1,i+1)/28.97;

    temp_area = (cp6(1,i)+cp6(1,i+1))/2;
    h6 = h6 + temp_area; 
end

s6 = 1.70203;
for i= 3000000 : 1 : 10000000 
    cp6_prime (1,i+1) = (alpha/(i/10000) + beta + gamma*((i/10000)^1) + lambda*((i/10000)^2) + epsilon*((i/10000)^3))*Rprime;
    cp6(1,i+1) = cp6_prime(1,i+1)/28.97;
    temp_area = (((cp6(1,i)+cp6(1,i+1))*0.0001)/2);
    s6 = s6 + temp_area;
end
s6 = s6 - R*log(P6/Pref);

k6=cp6(1,1000)/(cp6(1,1000)-R);

%% STATE 7

P7 = P1;
s7s = s6;

cp7s_prime = zeros(1,5000001);
cp7s_prime (1,1) = alpha;

cp7s = zeros(1,5000001);
cp7s(1,1) = cp7s_prime (1,1)/28.97;

temp_area = 0;

s7s_int = 1.70203-R*log(P7/Pref);
s7s_exp = s7s_int;
for i= 3000000 : 1 : 8000000
    cp7s_prime (1,i+1) = (alpha + beta*(i/10000) + gamma*((i/10000)^2) + lambda*((i/10000)^3) + epsilon*((i/10000)^4))*Rprime;
    cp7s(1,i+1) = cp7s_prime(1,i+1)/28.97;

    temp_area = (((cp7s(1,i)/(i/10000))+((cp7s(1,i+1))/(i/10000)))*0.0001)/2;
    s7s_exp = s7s_exp + temp_area;
    delta = s7s_exp - s7s ;

    if (delta > 0) 
        T7s = i/10000;
        break
    end
end

temp_area = 0;
h7s = 0;

for i= 1 : 1 : T7s*10000
    cp7s_prime (1,i+1) = (alpha + beta*(i/10000) + gamma*((i/10000)^2) + lambda*((i/10000)^3) + epsilon*((i/10000)^4))*Rprime;
    cp7s(1,i+1) = cp7s_prime(1,i+1)/28.97;
    temp_area = ((cp7s(1,i)+cp7s(1,i+1))*0.0001)/2;
    h7s = h7s + temp_area; 
end

h7 = h6 - ((h6-h7s)*nst);

cp7_prime = zeros(1,8000001);
cp7_prime (1,1) = alpha;

cp7 = zeros(1,8000001);
cp7(1,1) = cp7_prime (1,1)/28.97;

h7_exp = 0;

for i= 1 : 1 : 8000000
    cp7_prime (1,i+1) = (alpha + beta*(i/10000) + gamma*((i/10000)^2) + lambda*((i/10000)^3) + epsilon*((i/10000)^4))*Rprime;
    cp7(1,i+1) = cp7_prime(1,i+1)/28.97;

    temp_area = ((cp7(1,i)+cp7(1,i+1))*0.0001)/2;
    h7_exp = h7_exp + temp_area;
    delta = h7_exp - h7 ;

    if (delta > 0) 
        T7 = i/10000;
        break
    end
end

k7=cp7(1,T7*10000)/(cp7(1,T7*10000)-R);

temp_area = 0;
s7_int = 1.70203;
s7 = s4_int;
for i= 3000000 : 1 : T7*10000 
    cp7_prime (1,i+1) = (alpha/(i/10000) + beta + gamma*((i/10000)^1) + lambda*((i/10000)^2) + epsilon*((i/10000)^3))*Rprime;
    cp7(1,i+1) = cp7_prime(1,i+1)/28.97;
    temp_area = (((cp7(1,i)+cp7(1,i+1))*0.0001)/2);
    s7 = s7 + temp_area;
end
s7 = s7 - R*log(P7/Pref);

%% STATE 5
nsr = 0.88;
h5 = nsr*(h7-h4)+h4;
P5 = P4 ;

cp5_prime = zeros(1,6000001);
cp5_prime (1,1) = alpha;

cp5 = zeros(1,6000001);
cp5(1,1) = cp5_prime (1,1)/28.97;

h5_exp = 0;

for i= 1 : 1 : 6000000
    cp5_prime (1,i+1) = (alpha + beta*(i/10000) + gamma*((i/10000)^2) + lambda*((i/10000)^3) + epsilon*((i/10000)^4))*Rprime;
    cp5(1,i+1) = cp5_prime(1,i+1)/28.97;

    temp_area = ((cp5(1,i)+cp5(1,i+1))*0.0001)/2;
    h5_exp = h5_exp + temp_area;
    delta = h5_exp - h5 ;

    if (delta > 0) 
        T5 = i/10000;
        break
    end
end

k5=cp5(1,T5*10000)/(cp5(1,T5*10000)-R);

temp_area = 0;
s5_int = 1.70203;
s5 = s5_int;
for i= 3000000 : 1 : 5465529 
    cp5_prime (1,i+1) = (alpha/(i/10000) + beta + gamma*((i/10000)^1) + lambda*((i/10000)^2) + epsilon*((i/10000)^3))*Rprime;
    cp5(1,i+1) = cp5_prime(1,i+1)/28.97;
    temp_area = (((cp5(1,i)+cp5(1,i+1))*0.0001)/2);
    s5 = s5 + temp_area;
end
s5 = s5 - R*log(P5/Pref); 

%% STATE 8

h8 = -(h5-h4-h7);
P8 = P7;

cp8_prime = zeros(1,6000001);
cp8_prime (1,1) = alpha;

cp8 = zeros(1,6000001);
cp8(1,1) = cp8_prime (1,1)/28.97;

h8_exp = 0;

for i= 1 : 1 : 6000000
    cp8_prime (1,i+1) = (alpha + beta*(i/10000) + gamma*((i/10000)^2) + lambda*((i/10000)^3) + epsilon*((i/10000)^4))*Rprime;
    cp8(1,i+1) = cp8_prime(1,i+1)/28.97;

    temp_area = ((cp8(1,i)+cp8(1,i+1))*0.0001)/2;
    h8_exp = h8_exp + temp_area;
    delta = h8_exp - h8 ;

    if (delta > 0) 
        T8 = i/10000;
        break
    end
end

k8=cp8(1,T8*10000)/(cp8(1,T8*10000)-R);

temp_area = 0;
s8_int = 1.70203;
s8 = s8_int;
for i= 3000000 : 1 : 4896563 
    cp8_prime (1,i+1) = (alpha/(i/10000) + beta + gamma*((i/10000)^1) + lambda*((i/10000)^2) + epsilon*((i/10000)^3))*Rprime;
    cp8(1,i+1) = cp8_prime(1,i+1)/28.97;
    temp_area = (((cp8(1,i)+cp8(1,i+1))*0.0001)/2);
    s8 = s8 + temp_area;
end
s8 = s8 - R*log(P8/Pref);

%% Work and Heat Transfer Calculations

t_w = (h6 - h7)*100/105; 
c_w1 = (h2 - h1);
c_w2 = (h4 - h3);

w_cycle_hv = 110*10^3; % given from teacher

mdot = w_cycle_hv / (t_w - c_w1 - c_w2)

W_dot1 = -mdot*(h2-h1);
W_dot2 = 0;
W_dot3 = -mdot*(h4-h3);
W_dot4 = 0;
W_dot5 = 0;
W_dot6 = mdot*(h6-h7)*100/105;
W_dot7 = 0;
W_dot8 = 0;

W_cycle = W_dot1+W_dot3+W_dot6;
Q_dot_in = h6 - h5; 

Q_dot1 = 0;
Q_dot2 = -mdot*(h2-h3);
Q_dot3 = 0;
Q_dot4 = +mdot*(h5-h4);
Q_dot5 = +mdot*(h6-h5);
Q_dot6 = -mdot*((h6 - h7)*(5/105));
Q_dot7 = -mdot*(h7-h8);
Q_dot8 = -mdot*(h8-h1);
Q_cycle = Q_dot1+Q_dot2+Q_dot3+Q_dot4+Q_dot5+Q_dot6+Q_dot7+Q_dot8;

V_dot1 = R*T1*mdot/P1;
V_dot2 = R*T2*mdot/P2;
V_dot3 = R*T3*mdot/P3;
V_dot4 = R*T4*mdot/P4;
V_dot5 = R*T5*mdot/P5;
V_dot6 = R*T6*mdot/P6;
V_dot7 = R*T7*mdot/P7;
V_dot8 = R*T8*mdot/P8;

sigma_gen_dot1 = - (mdot*s1) + (mdot*s2) - (Q_dot1/T1);
sigma_gen_dot2 = - (mdot*s2) + (mdot*s3) - (Q_dot2/T1);
sigma_gen_dot3 = - (mdot*s3) + (mdot*s4) - (Q_dot3/T1);
sigma_gen_dot4 = - (mdot*s4) + (mdot*s5) - (Q_dot4*2/(T5+T8));
sigma_gen_dot5 = - (mdot*s5) + (mdot*s6) - (Q_dot5/T6);
sigma_gen_dot6 = - (mdot*s6) + (mdot*s7) - (Q_dot6/T7);
sigma_gen_dot7 = - (mdot*s7) + (mdot*s8) - (Q_dot7*2/(T5+T8));
sigma_gen_dot8 = - (mdot*s8) + (mdot*s1) - (Q_dot8/T1);

sigma_gen_dot_cycle = sigma_gen_dot1+sigma_gen_dot2+sigma_gen_dot3+sigma_gen_dot4+sigma_gen_dot5+sigma_gen_dot6+sigma_gen_dot7+sigma_gen_dot8;

eff = (t_w - c_w1 -c_w2)/Q_dot_in;

bwr = (c_w1+c_w2)/t_w;
spo = W_cycle / mdot;
%% STATE TABLE
State = [1;2;3;4;5;6;7;8];
Description = {'Low Pressure Compressor Inlet';'Intercooler Inlet';'High Pressure Compressor Inlet';'Regeneration Inlet';'Combustor Inlet';'Turbine Inlet';'Turbine Outlet';'Regeneration Outlet'};
T = [T1;T2;T3;T4;T5;T6;T7;T8];
P = [P1;P2;P3;P4;P5;P6;P7;P8];
m_dot = [mdot;mdot;mdot;mdot;mdot;mdot;mdot;mdot];
V_dot = [V_dot1;V_dot2;V_dot3;V_dot4;V_dot5;V_dot6;V_dot7;V_dot8];
h = [h1;h2;h3;h4;h5;h6;h7;h8];
s = [s1;s2;s3;s4;s5;s6;s7;s8];
k = [k1;k2;k3;k4;k5;k6;k7;k8];

TA1 = table(State,Description,T,P,m_dot,V_dot,h,s,k)

%% PROCESS TABLE
Process = {'1->2';'2->3';'3->4';'4->5';'5->6';'6->7';'7->8';'8->1';'Cycle'};
Description = {'Low Pressure Compressor';'Intercooler';'High Pressure Compressor';'Regeneration (cold)';'Combustor';'Turbine';'Regeneration (hot)';'Theoretical HX';'Entire Cycle'};
Q_dot = [Q_dot1;Q_dot2;Q_dot3;Q_dot4;Q_dot5;Q_dot6;Q_dot7;Q_dot8;Q_cycle];
W_dot = [W_dot1;W_dot2;W_dot3;W_dot4;W_dot5;W_dot6;W_dot7;W_dot8;W_cycle];
sigma_gen_dot = [sigma_gen_dot1;sigma_gen_dot2;sigma_gen_dot3;sigma_gen_dot4;sigma_gen_dot5;sigma_gen_dot6;sigma_gen_dot7;sigma_gen_dot8;sigma_gen_dot_cycle];
eta = {'-';'-';'-';'-';'-';'-';'-';'-';eff};
TA2 = table(Process,Description,Q_dot,W_dot,sigma_gen_dot,eta)