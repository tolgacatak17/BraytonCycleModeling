clc 
clear all

Tref1 = 0; % in kelvin
Pref = 101.325; % in kilo pascal
Tref2 = 300;
isregen = 0;
T_initial = [255 265 275 285 295 305 315]; %in kelvin
n=length(T_initial);

h1=zeros(1,n);
h2s=zeros(1,n);
h2=zeros(1,n);
h3=zeros(1,n);
h4s=zeros(1,n);
h4=zeros(1,n);
h5=zeros(1,n);
h6=zeros(1,n);
h7s=zeros(1,n);
h7=zeros(1,n);
h8=zeros(1,n);

s1=zeros(1,n);
s2s=zeros(1,n);
s2=zeros(1,n);
s3=zeros(1,n);
s4s=zeros(1,n);
s4=zeros(1,n);
s5=zeros(1,n);
s6=zeros(1,n);
s7s=zeros(1,n);
s7=zeros(1,n);
s8=zeros(1,n);


T1=zeros(1,n);
T2s=zeros(1,n);
T2=zeros(1,n);
T3=zeros(1,n);
T4s=zeros(1,n);
T4=zeros(1,n);
T5=zeros(1,n);
T6=zeros(1,n);
T7s=zeros(1,n);
T7=zeros(1,n);
T8=zeros(1,n);

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
%cv = cp - R ,thus we can obtain the k via k = cp/cv

for j = 1:1:n
    %% STATE 1
    T1(j) = T_initial(j); %in Kelvin
    P1 = 100; %in kilo pascal
    
    %k1 = cp1/cv1;
    cp1_prime = zeros(1,T_initial(j));
    cp1_prime (1,1) = alpha;
    
    cp1 = zeros(1,T_initial(j));
    cp1(1,1) = cp1_prime (1,1)/28.97;
    
    cv1 = zeros(1,T_initial(j)); 
    k1 = zeros(1,T_initial(j));
    
    h1(j) = 0;
    for i= 1 : 1 : (T_initial(j)-1)
        cp1_prime (1,i+1) = (alpha + beta*(0+i) + gamma*((0+i)^2) + lambda*((0+i)^3) + epsilon*((0+i)^4))*Rprime;
        cp1(1,i+1) = cp1_prime(1,i+1)/28.97;
    
        temp_area = (cp1(1,i)+cp1(1,i+1))/2;
        h1(j) = h1(j) + temp_area; 
    end
    
    index=round(T_initial(j));
    k1(j)=cp1(1,T_initial(j))/(cp1(1,T_initial(j))-R); % hesaplayip duzelt

    temp_area = 0;
    s1_int = 1.70203;
    s1(j) = s1_int;
    
    if T_initial(j) >= Tref2
        for i= Tref2*10000 : 1 : T_initial(j)*10000 
            cp1_prime (1,i+1) = (alpha/(i/10000) + beta + gamma*((i/10000)^1) + lambda*((i/10000)^2) + epsilon*((i/10000)^3))*Rprime;
            cp1(1,i+1) = cp1_prime(1,i+1)/28.97;
            temp_area = (((cp1(1,i)+cp1(1,i+1))*0.0001)/2);
            s1(j) = s1(j) + temp_area;
        end
        s1(j) = s1(j) - R*log(P1/Pref); 
    else 
        for i= T_initial(j)*10000 : 1 : Tref2*10000 
            cp1_prime (1,i+1) = (alpha/(i/10000) + beta + gamma*((i/10000)^1) + lambda*((i/10000)^2) + epsilon*((i/10000)^3))*Rprime;
            cp1(1,i+1) = cp1_prime(1,i+1)/28.97;
            temp_area = (((cp1(1,i)+cp1(1,i+1))*0.0001)/2);
            s1(j) = s1(j) - temp_area;
        end
        s1(j) = s1(j) - R*log(P1/Pref); 
    end
    
    
    %% STATE 2
    P2 = 4*P1;
    s2s(j) = s1(j);
    
    cp2s_prime = zeros(1,T_initial(j)*10000+1);
    cp2s_prime (1,1) = alpha;
    
    cp2s = zeros(1,T_initial(j)*10000+1);
    cp2s(1,1) = cp2s_prime (1,1)/28.97;
    
    temp_area = 0;
    
    s2s_int(j) = 1.70203-R*log(P2/Pref);
    s2s_exp(j) = s2s_int(j);
    for i= Tref2*10000 : 1 : 600*10000
        cp2s_prime (1,i+1) = (alpha + beta*(i/10000) + gamma*((i/10000)^2) + lambda*((i/10000)^3) + epsilon*((i/10000)^4))*Rprime;
        cp2s(1,i+1) = cp2s_prime(1,i+1)/28.97;
    
        temp_area = (((cp2s(1,i)/(i/10000))+((cp2s(1,i+1))/(i/10000)))*0.0001)/2;
        s2s_exp(j) = s2s_exp(j) + temp_area;
        delta = s2s_exp(j) - s2s(j) ;
    
        if (delta > 0) 
            T2s(j) = i/10000;
            break
        end
    end
    
    temp_area = 0;
    h2s(j) = 0;
    for i= 1 : 1 : T2s(j)*10000
        cp2s_prime (1,i+1) = (alpha + beta*(i/10000) + gamma*((i/10000)^2) + lambda*((i/10000)^3) + epsilon*((i/10000)^4))*Rprime;
        cp2s(1,i+1) = cp2s_prime(1,i+1)/28.97;
        temp_area = ((cp2s(1,i)+cp2s(1,i+1))*0.0001)/2;
        h2s(j) = h2s(j) + temp_area; 
    end
    
    h2(j) = h1(j) + (h2s(j)-h1(j))/nsc; 
    
    cp2_prime = zeros(1,6000001);
    cp2_prime (1,1) = alpha;
    
    cp2 = zeros(1,6000001);
    cp2(1,1) = cp2_prime (1,1)/28.97;
    
    h2_exp(j) = 0;
    
    for i= 1 : 1 : 600*10000
        cp2_prime (1,i+1) = (alpha + beta*(i/10000) + gamma*((i/10000)^2) + lambda*((i/10000)^3) + epsilon*((i/10000)^4))*Rprime;
        cp2(1,i+1) = cp2_prime(1,i+1)/28.97;
    
        temp_area = ((cp2(1,i)+cp2(1,i+1))*0.0001)/2;
        h2_exp(j) = h2_exp(j) + temp_area;
        delta = h2_exp(j) - h2(j) ;
    
        if (delta > 0) 
            T2(j) = i/10000;
            break
        end
    end
    
    index=round(T2(j)*10000);
    k2(j)=cp2(1,index)/(cp2(1,index)-R); % hesaplayip duzelt

    temp_area = 0;
    s2_int(j) = 1.70203;
    s2(j) = s2_int(j);


    if T2(j) >= Tref2
        for i= Tref2*10000 : 1 : T2(j)*10000 
            cp2_prime (1,i+1) = (alpha/(i/10000) + beta + gamma*((i/10000)^1) + lambda*((i/10000)^2) + epsilon*((i/10000)^3))*Rprime;
            cp2(1,i+1) = cp2_prime(1,i+1)/28.97;
            temp_area = (((cp2(1,i)+cp2(1,i+1))*0.0001)/2);
            s2(j) = s2(j) + temp_area;
        end
        s2(j) = s2(j) - R*log(P2/Pref);  
    else 
        for i= T2(j)*10000 : 1 : Tref2*10000 
            cp2_prime (1,i+1) = (alpha/(i/10000) + beta + gamma*((i/10000)^1) + lambda*((i/10000)^2) + epsilon*((i/10000)^3))*Rprime;
            cp2(1,i+1) = cp2_prime(1,i+1)/28.97;
            temp_area = (((cp2(1,i)+cp2(1,i+1))*0.0001)/2);
            s2(j) = s2(j) - temp_area;
        end
        s2(j) = s2(j) - R*log(P2/Pref);  
    end
    
    

    %% STATE 3
    
    T3(j) = T1(j); %in Kelvin
    P3 = P2; %in kilo pascal
    
    cp3_prime = zeros(1,T3(j));
    cp3_prime (1,1) = alpha;
    
    cp3 = zeros(1,T3(j));
    cp3(1,1) = cp3_prime (1,1)/28.97;
    
    h3(j) = 0;
    for i= 1 : 1 : T3(j)-1
        cp3_prime (1,i+1) = (alpha + beta*(0+i) + gamma*((0+i)^2) + lambda*((0+i)^3) + epsilon*((0+i)^4))*Rprime;
        cp3(1,i+1) = cp3_prime(1,i+1)/28.97;
    
        temp_area = (cp3(1,i)+cp3(1,i+1))/2;
        h3(j) = h3(j) + temp_area; 
    end
    
    index=round(T_initial);
    k3(j)=cp3(1,T_initial(j))/(cp3(1,T_initial(j))-R);
    
    temp_area = 0;
    s3_int(j) = 1.70203;
    s3(j) = s3_int(j);


    if T3(j) >= Tref2
        for i= Tref2*10000 : 1 : T3(j)*10000 
            cp3_prime (1,i+1) = (alpha/(i/10000) + beta + gamma*((i/10000)^1) + lambda*((i/10000)^2) + epsilon*((i/10000)^3))*Rprime;
            cp3(1,i+1) = cp3_prime(1,i+1)/28.97;
            temp_area = (((cp3(1,i)+cp3(1,i+1))*0.0001)/2);
            s3(j) = s3(j) + temp_area;
        end
        s3(j) = s3(j) - R*log(P3/Pref);  
    else 
        for i= T3(j)*10000 : 1 : Tref2*10000 
            cp3_prime (1,i+1) = (alpha/(i/10000) + beta + gamma*((i/10000)^1) + lambda*((i/10000)^2) + epsilon*((i/10000)^3))*Rprime;
            cp3(1,i+1) = cp3_prime(1,i+1)/28.97;
            temp_area = (((cp3(1,i)+cp3(1,i+1))*0.0001)/2);
            s3(j) = s3(j) - temp_area;
        end
        s3(j) = s3(j) - R*log(P2/Pref);  
    end

     % hesaplayip duzelt
    

    %% STATE 4
    
    P4 = 4*P3;
    s4s(j) = s3(j);
    
    cp4s_prime = zeros(1,(600-Tref2)*10000+1);
    cp4s_prime (1,1) = alpha;
    
    cp4s = zeros(1,(600-Tref2)*10000+1);
    cp4s(1,1) = cp4s_prime (1,1)/28.97;
    
    temp_area = 0;
    
    s4s_int(j) = 1.70203-R*log(P4/Pref);
    s4s_exp(j) = s4s_int(j);
    for i= Tref2*10000 : 1 : 600*10000
        cp4s_prime (1,i+1) = (alpha + beta*(i/10000) + gamma*((i/10000)^2) + lambda*((i/10000)^3) + epsilon*((i/10000)^4))*Rprime;
        cp4s(1,i+1) = cp4s_prime(1,i+1)/28.97;
    
        temp_area = (((cp4s(1,i)/(i/10000))+((cp4s(1,i+1))/(i/10000)))*0.0001)/2;
        s4s_exp(j) = s4s_exp(j) + temp_area;
        delta = s4s_exp(j) - s4s(j) ;
    
        if (delta > 0) 
            T4s(j) = i/10000;
            break
        end
    end
    
    temp_area = 0;
    h4s(j) = 0;
    
    for i= 1 : 1 : T4s(j)*10000
        cp4s_prime (1,i+1) = (alpha + beta*(i/10000) + gamma*((i/10000)^2) + lambda*((i/10000)^3) + epsilon*((i/10000)^4))*Rprime;
        cp4s(1,i+1) = cp4s_prime(1,i+1)/28.97;
        temp_area = ((cp4s(1,i)+cp4s(1,i+1))*0.0001)/2;
        h4s(j) = h4s(j) + temp_area; 
    end
    
    h4(j) = h3(j) + (h4s(j)-h3(j))/nsc; 
    
    cp4_prime = zeros(1,6000001);
    cp4_prime (1,1) = alpha;
    
    cp4 = zeros(1,6000001);
    cp4(1,1) = cp4_prime (1,1)/28.97;
    
    h4_exp(j) = 0;
    
    for i= 1 : 1 : 600*10000
        cp4_prime (1,i+1) = (alpha + beta*(i/10000) + gamma*((i/10000)^2) + lambda*((i/10000)^3) + epsilon*((i/10000)^4))*Rprime;
        cp4(1,i+1) = cp4_prime(1,i+1)/28.97;
    
        temp_area = ((cp4(1,i)+cp4(1,i+1))*0.0001)/2;
        h4_exp(j) = h4_exp(j) + temp_area;
        delta = h4_exp(j) - h4(j) ;
    
        if (delta > 0) 
            T4(j) = i/10000;
            break
        end
    end
    
    index=round(T4(j)*10000);
    k4(j)=cp4(1,index)/(cp4(1,index)-R); % hesaplayip duzelt

    temp_area = 0;
    s4_int(j) = 1.70203;
    s4(j) = s4_int(j);


    if T4(j) >= Tref2
        for i= Tref2*10000 : 1 : T4(j)*10000 
            cp4_prime (1,i+1) = (alpha/(i/10000) + beta + gamma*((i/10000)^1) + lambda*((i/10000)^2) + epsilon*((i/10000)^3))*Rprime;
            cp4(1,i+1) = cp4_prime(1,i+1)/28.97;
            temp_area = (((cp4(1,i)+cp4(1,i+1))*0.0001)/2);
            s4(j) = s4(j) + temp_area;
        end
        s4(j) = s4(j) - R*log(P4/Pref);  
    else 
        for i= T4(j)*10000 : 1 : Tref2*10000 
            cp4_prime (1,i+1) = (alpha/(i/10000) + beta + gamma*((i/10000)^1) + lambda*((i/10000)^2) + epsilon*((i/10000)^3))*Rprime;
            cp4(1,i+1) = cp4_prime(1,i+1)/28.97;
            temp_area = (((cp4(1,i)+cp4(1,i+1))*0.0001)/2);
            s4(j) = s4(j) - temp_area;
        end
        s4(j) = s4(j) - R*log(P4/Pref);  
    end
    
 
    

    %% STATE 6
    
    T6(j) = 1000; %in Kelvin
    P6 = P4*0.95; %in kilo pascal
    
    cp6_prime = zeros(1,1000);
    cp6_prime (1,1) = alpha;
    
    cp6 = zeros(1,1000);
    cp6(1,1) = cp6_prime (1,1)/28.97;

    h6(j) = 0;
    for i= 1 : 1 : 999
        cp6_prime (1,i+1) = (alpha + beta*(0+i) + gamma*((0+i)^2) + lambda*((0+i)^3) + epsilon*((0+i)^4))*Rprime;
        cp6(1,i+1) = cp6_prime(1,i+1)/28.97;
    
        temp_area = (cp6(1,i)+cp6(1,i+1))/2;
        h6(j) = h6(j) + temp_area; 
    end
    
    s6(j) = 1.70203;
    for i= Tref2*10000 : 1 : 1000*10000 
        cp6_prime (1,i+1) = (alpha/(i/10000) + beta + gamma*((i/10000)^1) + lambda*((i/10000)^2) + epsilon*((i/10000)^3))*Rprime;
        cp6(1,i+1) = cp6_prime(1,i+1)/28.97;
        temp_area = (((cp6(1,i)+cp6(1,i+1))*0.0001)/2);
        s6(j) = s6(j) + temp_area;
    end
    s6(j) = s6(j) - R*log(P6/Pref);
    
    k6(j)=cp6(1,1000)/(cp6(1,1000)-R); % hesaplayip duzelt
    

    %% STATE 7
    
    P7 = P1;
    s7s(j) = s6(j);
    
    cp7s_prime = zeros(1,5000001);
    cp7s_prime (1,1) = alpha;
    
    cp7s = zeros(1,5000001);
    cp7s(1,1) = cp7s_prime (1,1)/28.97;
    
    temp_area = 0;
    
    s7s_int(j) = 1.70203-R*log(P7/Pref);
    s7s_exp(j) = s7s_int(j);
    for i= Tref2*10000 : 1 : 800*10000
        cp7s_prime (1,i+1) = (alpha + beta*(i/10000) + gamma*((i/10000)^2) + lambda*((i/10000)^3) + epsilon*((i/10000)^4))*Rprime;
        cp7s(1,i+1) = cp7s_prime(1,i+1)/28.97;
    
        temp_area = (((cp7s(1,i)/(i/10000))+((cp7s(1,i+1))/(i/10000)))*0.0001)/2;
        s7s_exp(j) = s7s_exp(j) + temp_area;
        delta = s7s_exp(j) - s7s(j) ;
    
        if (delta > 0) 
            T7s(j) = i/10000;
            break
        end
    end
    
    temp_area = 0;
    h7s(j) = 0;
    
    for i= 1 : 1 : T7s(j)*10000
        cp7s_prime (1,i+1) = (alpha + beta*(i/10000) + gamma*((i/10000)^2) + lambda*((i/10000)^3) + epsilon*((i/10000)^4))*Rprime;
        cp7s(1,i+1) = cp7s_prime(1,i+1)/28.97;
        temp_area = ((cp7s(1,i)+cp7s(1,i+1))*0.0001)/2;
        h7s(j) = h7s(j) + temp_area; 
    end
    
    h7(j) = h6(j) - ((h6(j)-h7s(j))*nst);
    
    cp7_prime = zeros(1,8000001);
    cp7_prime (1,1) = alpha;
    
    cp7 = zeros(1,8000001);
    cp7(1,1) = cp7_prime (1,1)/28.97;
    
    h7_exp(j) = 0;
    
    for i= 1 : 1 : 800*10000
        cp7_prime (1,i+1) = (alpha + beta*(i/10000) + gamma*((i/10000)^2) + lambda*((i/10000)^3) + epsilon*((i/10000)^4))*Rprime;
        cp7(1,i+1) = cp7_prime(1,i+1)/28.97;
    
        temp_area = ((cp7(1,i)+cp7(1,i+1))*0.0001)/2;
        h7_exp(j) = h7_exp(j) + temp_area;
        delta = h7_exp(j) - h7(j) ;
    
        if (delta > 0) 
            T7(j) = i/10000;
            break
        end
    end
    
    index=round(T7(j)*10000);
    k7(j)=cp7(1,index)/(cp7(1,index)-R);

    temp_area = 0;
    s7_int(j) = 1.70203;
    s7(j) = s4_int(j);
    for i= Tref2*10000 : 1 : T7(j)*10000 
        cp7_prime (1,i+1) = (alpha/(i/10000) + beta + gamma*((i/10000)^1) + lambda*((i/10000)^2) + epsilon*((i/10000)^3))*Rprime;
        cp7(1,i+1) = cp7_prime(1,i+1)/28.97;
        temp_area = (((cp7(1,i)+cp7(1,i+1))*0.0001)/2);
        s7(j) = s7(j) + temp_area;
    end
    s7(j) = s7(j) - R*log(P7/Pref);
    

    %% STATE 5
    if isregen  
        nsr = 0.88;
        h5(j) = nsr*(h7(j)-h4(j))+h4(j);
        P5 = P4 ;
        
        cp5_prime = zeros(1,6000001);
        cp5_prime (1,1) = alpha;
        
        cp5 = zeros(1,6000001);
        cp5(1,1) = cp5_prime (1,1)/28.97;
        
        h5_exp(j) = 0;
        
        for i= 1 : 1 : 600*10000
            cp5_prime (1,i+1) = (alpha + beta*(i/10000) + gamma*((i/10000)^2) + lambda*((i/10000)^3) + epsilon*((i/10000)^4))*Rprime;
            cp5(1,i+1) = cp5_prime(1,i+1)/28.97;
        
            temp_area = ((cp5(1,i)+cp5(1,i+1))*0.0001)/2;
            h5_exp(j) = h5_exp(j) + temp_area;
            delta = h5_exp(j) - h5(j) ;
        
            if (delta > 0) 
                T5(j) = i/10000;
                break
            end
        end
        
        index=round(T5(j)*10000);
        k5(j)=cp5(1,index)/(cp5(1,index)-R); % hesaplayip duzelt

        temp_area = 0;
        s5_int(j) = 1.70203;
        s5(j) = s5_int(j);
        for i= Tref2*10000 : 1 : T5(j)*10000 
            cp5_prime (1,i+1) = (alpha/(i/10000) + beta + gamma*((i/10000)^1) + lambda*((i/10000)^2) + epsilon*((i/10000)^3))*Rprime;
            cp5(1,i+1) = cp5_prime(1,i+1)/28.97;
            temp_area = (((cp5(1,i)+cp5(1,i+1))*0.0001)/2);
            s5(j) = s5(j) + temp_area;
        end
        s5(j) = s5(j) - R*log(P5/Pref); 
     
    else
        T5(j)=T4(j);
        P5=P4;
        h5(j)=h4(j);
        s5(j)=s4(j);
        k5(j)=k4(j); % hesaplayip duzelt 4u duzeltsen yeter
    end
    %% STATE 8
    if isregen
        h8(j) = -(h5(j)-h4(j)-h7(j));
        P8 = P7;
        
        cp8_prime = zeros(1,6000001);
        cp8_prime (1,1) = alpha;
        
        cp8 = zeros(1,6000001);
        cp8(1,1) = cp8_prime (1,1)/28.97;
        
        h8_exp(j) = 0;
        
        for i= 1 : 1 : 600*10000
            cp8_prime (1,i+1) = (alpha + beta*(i/10000) + gamma*((i/10000)^2) + lambda*((i/10000)^3) + epsilon*((i/10000)^4))*Rprime;
            cp8(1,i+1) = cp8_prime(1,i+1)/28.97;
        
            temp_area = ((cp8(1,i)+cp8(1,i+1))*0.0001)/2;
            h8_exp(j) = h8_exp(j) + temp_area;
            delta = h8_exp(j) - h8(j) ;
        
            if (delta > 0) 
                T8(j) = i/10000;
                break
            end
        end
        
        index=round(T8(j)*10000);
        k8(j)=cp8(1,index)/(cp8(1,index)-R);

        temp_area = 0;
        s8_int(j) = 1.70203;
        s8(j) = s8_int(j);
        for i= Tref2*10000 : 1 : T8(j)*10000 
            cp8_prime (1,i+1) = (alpha/(i/10000) + beta + gamma*((i/10000)^1) + lambda*((i/10000)^2) + epsilon*((i/10000)^3))*Rprime;
            cp8(1,i+1) = cp8_prime(1,i+1)/28.97;
            temp_area = (((cp8(1,i)+cp8(1,i+1))*0.0001)/2);
            s8(j) = s8(j) + temp_area;
        end
        s8(j) = s8(j) - R*log(P8/Pref);

         % hesaplayip duzelt

    else
        T8(j)=T7(j);
        P8=P7;
        h8(j)=h7(j);
        s8(j)=s7(j);
        k8(j)=k7(j); % hesaplayip duzelt 7yi duzeltsen yeter
    end

    
    t_w(j) = (h6(j) - h7(j))*0.95; 
    c_w1(j) = (h2(j) - h1(j));
    c_w2(j) = (h4(j) - h3(j));

    w_cycle = 110*10^3;
    Q_dot_in(j) = h6(j) - h5(j);
   
    m_dot(j)=w_cycle / (t_w(j) - c_w1(j) - c_w2(j)); %hesaplatip duzelt mdot hepsinde ayni
    V_dot1(j)=R*T1(j)*m_dot(j)/P1; %hesaplayip duzelt
    V_dot2(j)=R*T2(j)*m_dot(j)/P2; %hesaplayip duzelt
    V_dot3(j)=R*T3(j)*m_dot(j)/P3; %hesaplayip duzelt
    V_dot4(j)=R*T4(j)*m_dot(j)/P4; %hesaplayip duzelt
    V_dot6(j)=R*T6(j)*m_dot(j)/P6; %hesaplayip duzelt
    V_dot7(j)=R*T7(j)*m_dot(j)/P7; %hesaplayip duzelt

    if isregen
        V_dot5(j)=R*T5(j)*m_dot(j)/P5; %hesaplayip duzelt
        V_dot8(j)=R*T8(j)*m_dot(j)/P8; %hesaplayip duzelt
    else
        V_dot5(j)=V_dot4(j); %hesaplayip duzelt 4u duzeltsen yeter
        V_dot8(j)=V_dot7(j); %hesaplayip duzelt 7yi duzeltsen yeter
    end
    
    W_dot1(j) = -m_dot(j)*(h2(j)-h1(j));
    W_dot2(j) = 0;
    W_dot3(j) = -m_dot(j)*(h4(j)-h3(j));
    W_dot4(j) = 0;
    W_dot5(j) = 0;
    W_dot6(j) = m_dot(j)*(h6(j)-h7(j))*100/105;
    W_dot7(j) = 0;
    W_dot8(j) = 0;
    
    W_cycle(j) = W_dot1(j)+W_dot3(j)+W_dot6(j);
    Q_dot_in(j) = h6(j) - h5(j); 
    
    Q_dot1(j) = 0;
    Q_dot2(j) = -m_dot(j)*(h2(j)-h3(j));
    Q_dot3(j) = 0;
    Q_dot4(j) = +m_dot(j)*(h5(j)-h4(j));
    Q_dot5(j) = +m_dot(j)*(h6(j)-h5(j));
    Q_dot6(j) = -m_dot(j)*((h6(j) - h7(j))*(5/105));
    Q_dot7(j) = -m_dot(j)*(h7(j)-h8(j));
    Q_dot8(j) = -m_dot(j)*(h8(j)-h1(j));
    Q_cycle(j) = Q_dot1(j)+Q_dot2(j)+Q_dot3(j)+Q_dot4(j)+Q_dot5(j)+Q_dot6(j)+Q_dot7(j)+Q_dot8(j);
    
    sigma_gen_dot1(j) = - (m_dot(j)*s1(j)) + (m_dot(j)*s2(j)) - (Q_dot1(j)/T1(j));
    sigma_gen_dot2(j) = - (m_dot(j)*s2(j)) + (m_dot(j)*s3(j)) - (Q_dot2(j)/T1(j));
    sigma_gen_dot3(j) = - (m_dot(j)*s3(j)) + (m_dot(j)*s4(j)) - (Q_dot3(j)/T1(j));
    sigma_gen_dot4(j) = - (m_dot(j)*s4(j)) + (m_dot(j)*s5(j)) - (Q_dot4(j)*2/(T5(j)+T8(j)));
    sigma_gen_dot5(j) = - (m_dot(j)*s5(j)) + (m_dot(j)*s6(j)) - (Q_dot5(j)/T6(j));
    sigma_gen_dot6(j) = - (m_dot(j)*s6(j)) + (m_dot(j)*s7(j)) - (Q_dot6(j)/T7(j));
    sigma_gen_dot7(j) = - (m_dot(j)*s7(j)) + (m_dot(j)*s8(j)) - (Q_dot7(j)*2/(T5(j)+T8(j)));
    sigma_gen_dot8(j) = - (m_dot(j)*s8(j)) + (m_dot(j)*s1(j)) - (Q_dot8(j)/T1(j));
    
    sigma_gen_dot_cycle = sigma_gen_dot1+sigma_gen_dot2+sigma_gen_dot3+sigma_gen_dot4+sigma_gen_dot5+sigma_gen_dot6+sigma_gen_dot7+sigma_gen_dot8;

    if isregen
        eta(j) = (t_w(j) - c_w1(j) -c_w2(j))/Q_dot_in(j);
        bwr(j) = (c_w1(j)+c_w2(j))/t_w(j);
        spo(j) = 110*10^3 / m_dot(j);
    else
    
    end
    eta(j) = (t_w(j) - c_w1(j) -c_w2(j))/Q_dot_in(j);
    bwr(j) = (c_w1(j)+c_w2(j))/t_w(j);
    spo(j) = 110*10^3 / m_dot(j);


        State = [1;2;3;4;5;6;7;8];
        Description = {'Low Pressure Compressor Inlet';'Intercooler Inlet';'High Pressure Compressor Inlet';'Regeneration Inlet';'Combustor Inlet';'Turbine Inlet';'Turbine Outlet';'Regeneration Outlet'};
        T = [T1(j);T2(j);T3(j);T4(j);T5(j);T6(j);T7(j);T8(j)];
        P = [P1;P2;P3;P4;P5;P6;P7;P8];
        m_dot_table = [m_dot(j);m_dot(j);m_dot(j);m_dot(j);m_dot(j);m_dot(j);m_dot(j);m_dot(j)];
        V_dot = [V_dot1(j);V_dot2(j);V_dot3(j);V_dot4(j);V_dot5(j);V_dot6(j);V_dot7(j);V_dot8(j)];
        h = [h1(j);h2(j);h3(j);h4(j);h5(j);h6(j);h7(j);h8(j)];
        s = [s1(j);s2(j);s3(j);s4(j);s5(j);s6(j);s7(j);s8(j)];
        k = [k1(j);k2(j);k3(j);k4(j);k5(j);k6(j);k7(j);k8(j)];
        
        TA1 = table(State,Description,T,P,m_dot_table,V_dot,h,s,k)

        %% PROCESS TABLE
        Process = {'1->2';'2->3';'3->4';'4->5';'5->6';'6->7';'7->8';'8->1';'Cycle'};
        Description = {'Low Pressure Compressor';'Intercooler';'High Pressure Compressor';'Regeneration (cold)';'Combustor';'Turbine';'Regeneration (hot)';'Theoretical HX';'Entire Cycle'};
        Q_dot = [Q_dot1(j);Q_dot2(j);Q_dot3(j);Q_dot4(j);Q_dot5(j);Q_dot6(j);Q_dot7(j);Q_dot8(j);Q_cycle(j)];
        W_dot = [W_dot1(j);W_dot2(j);W_dot3(j);W_dot4(j);W_dot5(j);W_dot6(j);W_dot7(j);W_dot8(j);W_cycle(j)];
        sigma_gen_dot = [sigma_gen_dot1(j);sigma_gen_dot2(j);sigma_gen_dot3(j);sigma_gen_dot4(j);sigma_gen_dot5(j);sigma_gen_dot6(j);sigma_gen_dot7(j);sigma_gen_dot8(j);sigma_gen_dot_cycle(j)];
        eta_table = {'-';'-';'-';'-';'-';'-';'-';'-';eta(j)};
        TA2 = table(Process,Description,Q_dot,W_dot,sigma_gen_dot,eta_table)
    

%%

figure(1)
plot(T1,eta)
grid on
xlabel("Temperature")
ylabel("Cycle Efficiency")

figure(2)
grid on
plot(T1,bwr)
xlabel("Temperature")
ylabel("Backwork Ratio")

figure(3)
grid on
plot(T1,spo)
xlabel("Temperature")
ylabel("Specific Power Output")

figure(4)
grid on
plot(T1,V_dot1)
xlabel("Temperature")
ylabel("Engine Inlet Volumetric Flow Rate")