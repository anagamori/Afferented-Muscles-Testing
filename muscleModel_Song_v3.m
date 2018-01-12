%--------------------------------------------------------------------------
% muscleModel_Loeb.m
% Author: Akira Nagamori
% Last update: 12/20/207
%--------------------------------------------------------------------------

function output = muscleModel_Song_v3(t,Fs,input,modelParameter)

% model parameters
alpha = modelParameter.pennationAngle; % pennation angle
Lm_initial = modelParameter.muscleInitialLength; % muscle initial length
Lt_initial = modelParameter.tendonInitialLength; % tendon initial length
Lmt = Lm_initial*cos(alpha)+Lt_initial; % intial musculotendon length
L0 = modelParameter.optimalLength; % optimal muscle length
L_tendon = modelParameter.tendonSlackLength; % tendon slack length
L0T = L_tendon*1.05; % optimal tendon length

[Lce,Lse,Lmax] =  InitialLength(modelParameter);

density = 1.06;
mass = modelParameter.mass; % muscle mass [kg]
PCSA = (mass*1000)/(density*L0); % PCSA of muscle
sigma = 31.8; % specific tension
F0 = PCSA * sigma; % maximal force

Ur = 0.8; % fractional activation level at which all motor units for a given muscle are recruited (0-1)
F_pcsa_slow = 0.5; % fractional PSCA of slow-twitch motor units (0-1)
F_pcsa_fast = 1-F_pcsa_slow; % fractional PSCA of fast-twitch motor units (0-1)
U1_th = 0.01; % threshold for slow-twitch fiber
U2_th = Ur*F_pcsa_slow;

% activation-frequency relationship (Brown and Loeb 2000)
f_half = 8.5; % frequency at which the motor unit produces half of its maximal isometric force
fmin = 0.5*f_half; % minimum firing frequency of slow-twitch fiber
fmax = 2*f_half; % maximum firing frequency of slow-twitch fiber

f_half_fast = 34;% frequency at which the motor unit produces half of its maximal isometric force
fmin_fast = 0.5*f_half_fast; % minimum firing frequency of fast-twitch fiber
fmax_fast = 2*f_half_fast; % maximum firing frequency of fast-twitch fiber

% simulation parameters
% 'isometric' = Vce = 0

U = input;

%--------------------------------------------------------------------------
% parameter initilization
U_eff = 0;
f_int_slow = 0;
f_eff_slow = 0;
f_eff_slow_dot = 0;
f_int_fast = 0;
f_eff_fast = 0;
f_eff_fast_dot = 0;

Af_slow = 0;
Af_fast = 0;

Y = 0;
S = 0;

Vce = 0; % muscle excursion velocity

%--------------------------------------------------------------------------
% storing variables
f_env_slow_vec = zeros(1,length(t));
f_int_slow_vec = zeros(1,length(t));
f_eff_slow_vec = zeros(1,length(t));
f_eff_slow_dot_vec = zeros(1,length(t));
f_env_fast_vec = zeros(1,length(t));
f_int_fast_vec = zeros(1,length(t));
f_eff_fast_vec = zeros(1,length(t));
f_eff_fast_dot_vec = zeros(1,length(t));

Af_slow_vec = zeros(1,length(t));
Af_fast_vec = zeros(1,length(t));

Y_vec = zeros(1,length(t));
S_vec = zeros(1,length(t));

Force = zeros(1,length(t));
Force_slow = zeros(1,length(t));
Force_fast = zeros(1,length(t));

OutputLse = zeros(1,length(t));
OutputLce = zeros(1,length(t));
OutputVce = zeros(1,length(t));
OutputAce = zeros(1,length(t));
MuscleAcceleration = zeros(1,length(t));
MuscleVelocity = zeros(1,length(t));
MuscleLength = zeros(1,length(t));
MuscleLength(1) = Lce*L0/100;

%--------------------------------------------------------------------------
h = 1/Fs;
%--------------------------------------------------------------------------
% simulation
for i = 1:length(t)
    if U(i) >= U_eff
        T_U = 0.03;
    elseif U(i) < U_eff
        T_U = 0.15;
    end
    
    U_eff_dot = (U(i) - U_eff)/T_U;
    U_eff = U_eff_dot*1/Fs + U_eff;
    
    if U_eff < U1_th
        W1 = 0;
    elseif U_eff < U2_th
        W1 = (U_eff - U1_th)/(U_eff - U1_th);
    else
        W1 = (U_eff - U1_th)/((U_eff - U1_th) + (U_eff - U2_th));
    end
    if U_eff < U2_th
        W2 = 0;
    else
        W2 = (U_eff - U2_th)/((U_eff - U1_th) + (U_eff - U2_th));
    end
    
    % firing frequency input to second-order excitation dynamics of
    % slow-twitch fiber
    if U_eff >=  U1_th
        f_env_slow = (fmax-fmin)/(1-U1_th).*(U_eff-U1_th)+fmin;
        f_env_slow = f_env_slow/f_half;
    else
        f_env_slow = 0;
    end
    
    [f_int_slow,~] = f_slow_function(f_int_slow,f_env_slow,f_env_slow,f_eff_slow_dot,Af_slow,Lce,Fs);
    [f_eff_slow,f_eff_slow_dot] = f_slow_function(f_eff_slow,f_int_slow,f_env_slow,f_eff_slow_dot,Af_slow,Lce,Fs);
    
    if U_eff >= U2_th
        f_env_fast = (fmax_fast-fmin_fast)/(1-U2_th).*(U_eff-U2_th)+fmin_fast;
        f_env_fast = f_env_fast/f_half_fast;
    else
        f_env_fast = 0;
    end
    
    [f_int_fast,~] = f_fast_function(f_int_fast,f_env_fast,f_env_fast,f_eff_fast_dot,Af_fast,Lce,Fs);
    [f_eff_fast,f_eff_fast_dot] = f_fast_function(f_eff_fast,f_int_fast,f_env_fast,f_eff_fast_dot,Af_fast,Lce,Fs);
    
    Y = yield_function(Y,Vce,Fs);
    S = sag_function(S,f_eff_fast,Fs);
    Af_slow = Af_slow_function(f_eff_slow,Lce,Y);
    Af_fast = Af_fast_function(f_eff_fast,Lce,S);
    
    % activation dependent force of contractile elements
    
    
    
    f_env_slow_vec(i) = f_env_slow;
    f_int_slow_vec(i) = f_int_slow;
    f_eff_slow_vec(i) = f_eff_slow;
    f_eff_slow_dot_vec(i) = f_eff_slow_dot;
    f_env_fast_vec(i) = f_env_fast;
    f_int_fast_vec(i) = f_int_fast;
    f_eff_fast_vec(i) = f_eff_fast;
    f_eff_fast_dot_vec(i) = f_eff_fast_dot;
    
    Af_slow_vec(i) = Af_slow;
    Af_fast_vec(i) = Af_fast;
    
    Y_vec(i) = Y;
    S_vec(i) = S;
    
    if Vce <= 0 % concentric
        FV1 = FVcon(Lce,Vce);
        FV2 = FVcon_fast(Lce,Vce);
    elseif Vce > 0 % eccentric
        FV1 = FVecc(Lce,Vce);
        FV2 = FVecc_fast(Lce,Vce);
    end
    FL1 = FL(Lce);
    FL2 = FL_fast(Lce);
    FP1 = Fpe1(Lce/Lmax,Vce);
    % passive element 2
    FP2 = Fpe2(Lce);
    if FP2 > 0
        FP2 = 0;
    end
    
    Fce = U_eff*((W1*Af_slow*(FL1*FV1+FP2)) + (W2*Af_fast*(FL2*FV2+FP2)));
    if Fce < 0
        Fce = 0;
    elseif Fce > 1
        Fce = 1;
    end
    Fce = Fce + FP1;
    Force(i) = Fce*F0;
    Force_slow(i) = W1*Af_slow*F0;
    Force_fast(i) = W2*Af_fast*F0;
    
    ForceSE = Fse(Lse)*F0;
    ForceCE = Force(i);
    
    [t,y] = ode45(@contractionDynamics,[0 1/Fs],[MuscleLength(i);MuscleVelocity(i)]);
    MuscleVelocity(i+1) = y(end,2);
    MuscleLength(i+1) = y(end,1);
   
    %Ace = MuscleAcceleration(i+1)/(L0/100);
    Vce = MuscleVelocity(i+1)/(L0/100);
    Lce = MuscleLength(i+1)/(L0/100);
    Lse = (Lmt - Lce*L0*cos(alpha))/L0T;
    
    OutoutForceTendon(i) = ForceSE;
    OutputLse(i) = Lse; % normalized tendon length
    OutputLce(i) = Lce; % normalized muscle length
    OutputVce(i) = Vce; % normalized muscle excursion velocity
    %OutputAce(i) = Ace; % normalized muscle excursion acceleration
end

output.Force_tendon = OutoutForceTendon;
output.Force_total = Force;
output.Force_slow = Force_slow;
output.Force_fast = Force_fast;
output.f_env_slow = f_env_slow_vec;
output.f_env_fast = f_env_fast_vec;
output.f_int_slow = f_int_slow_vec;
output.f_int_fast = f_int_fast_vec;
output.f_eff_slow = f_eff_slow_vec;
output.f_eff_fast = f_eff_fast_vec;
output.S = S_vec;
output.Y = Y_vec;
output.Af_slow = Af_slow_vec;
output.Af_fast = Af_fast_vec;
output.Lce = OutputLce;
output.Vce = OutputVce;

%--------------------------------------------------------------------------
% function used in simulation
    function Y = yield_function(Y,V,Fs)
        c_y = 0.35;
        V_y = 0.1;
        T_y = 0.2;
        
        Y_dot = (1-c_y*(1-exp(-abs(V)/V_y))-Y)/T_y;
        Y = Y_dot*1/Fs + Y;
    end

    function S = sag_function(S,f_eff,Fs)
        if f_eff < 0.1
            a_s = 1.76;
        else
            a_s = 0.96;
        end
        T_s = 0.043;
        S_dot = (a_s-S)/T_s;
        S = S_dot*1/Fs + S;
    end

    function Af = Af_slow_function(f_eff,L,Y)
        a_f = 0.56;
        n_f0 = 2.1;
        n_f1 = 5;
        n_f = n_f0 + n_f1*(1/L-1);
        Af = 1 - exp(-(Y*f_eff/(a_f*n_f))^n_f);
    end

    function Af = Af_fast_function(f_eff,L,S)
        a_f = 0.56;
        n_f0 = 2.1;
        n_f1 = 3.3;
        n_f = n_f0 + n_f1*(1/L-1);
        Af = 1 - exp(-((S*f_eff)/(a_f*n_f))^n_f);
    end

    function [f_out,f_out_dot,Tf] = f_slow_function(f_out,f_in,f_env,f_eff_dot,Af,Lce,Fs)
        % T_f1 = 0.0484;
        % T_f2 = 0.032;
        % T_f3 = 0.0664;
        % T_f4 = 0.0356;
        
        T_f1 = 0.0343;
        T_f2 = 0.0227;
        T_f3 = 0.047;
        T_f4 = 0.0252;
        
        if f_eff_dot >= 0
            Tf = T_f1*Lce^2+T_f2*f_env;
        else
            Tf = (T_f3 + T_f4*Af)/Lce;
        end
        f_out_dot = (f_in - f_out)/Tf;
        f_out = f_out_dot*1/Fs + f_out;
        
    end

    function [f_out,f_out_dot,Tf] = f_fast_function(f_out,f_in,f_env,f_eff_dot,Af,Lce,Fs)
        % T_f1 = 0.0121;
        % T_f2 = 0.008;
        % T_f3 = 0.0166;
        % T_f4 = 0.0089;
        
        T_f1 = 0.0206;
        T_f2 = 0.0136;
        T_f3 = 0.0282;
        T_f4 = 0.0151;
        
        if f_eff_dot >= 0
            Tf = T_f1*Lce^2+T_f2*f_env;
        else
            Tf = (T_f3 + T_f4*Af)/Lce;
        end
        f_out_dot = (f_in - f_out)/Tf;
        f_out = f_out_dot*1/Fs + f_out;
        
    end

    function FL = FL(L)
        %---------------------------
        % force length (F-L) relationship for slow-tiwtch fiber
        % input: normalized muscle length and velocity
        % output: F-L factor (0-1)
        %---------------------------
        beta = 2.3;
        omega = 1.12;
        rho = 1.62;
        
        FL = exp(-abs((L^beta - 1)/omega)^rho);
    end

    function FL = FL_fast(L)
        %---------------------------
        % force length (F-L) relationship for fast-twitch fiber
        % input: normalized muscle length and velocity
        % output: F-L factor (0-1)
        %---------------------------
        beta = 1.55;
        omega = 0.75;
        rho = 2.12;
        
        FL = exp(-abs((L^beta - 1)/omega)^rho);
    end

    function FVcon = FVcon(L,V)
        %---------------------------
        % concentric force velocity (F-V) relationship for slow-twitch fiber
        % input: normalized muscle length and velocity
        % output: F-V factor (0-1)
        %---------------------------
        Vmax = -7.88;
        cv0 = 5.88;
        cv1 = 0;
        
        
        FVcon = (Vmax - V)/(Vmax + (cv0 + cv1*L)*V);
    end

    function FVcon = FVcon_fast(L,V)
        %---------------------------
        % concentric force velocity (F-V) relationship for fast-twitch fiber
        % input: normalized muscle length and velocity
        % output: F-V factor (0-1)
        %---------------------------
        Vmax = -9.15;
        cv0 = -5.7;
        cv1 = 9.18;
        
        FVcon = (Vmax - V)/(Vmax + (cv0 + cv1*L)*V);
    end

    function FVecc = FVecc(L,V)
        %---------------------------
        % eccentric force velocity (F-V) relationship for slow-twitch fiber
        % input: normalized muscle length and velocity
        % output: F-V factor (0-1)
        %---------------------------
        av0 = -4.7;
        av1 = 8.41;
        av2 = -5.34;
        bv = 0.35;
        FVecc = (bv - (av0 + av1*L + av2*L^2)*V)/(bv+V);
    end

    function FVecc = FVecc_fast(L,V)
        %---------------------------
        % eccentric force velocity (F-V) relationship for fast-twitch fiber
        % input: normalized muscle length and velocity
        % output: F-V factor (0-1)
        %---------------------------
        av0 = -1.53;
        av1 = 0;
        av2 = 0;
        bv = 0.69;
        FVecc = (bv - (av0 + av1*L + av2*L^2)*V)/(bv+V);
    end

    function Fpe1 = Fpe1(L,V)
        %---------------------------
        % passive element 1
        % input: normalized muscle length
        % output: passive element force (0-1)
        %---------------------------
        c1_pe1 = 23;
        k1_pe1 = 0.046;
        Lr1_pe1 = 1.17;
        eta = 0.01;
        
        Fpe1 = c1_pe1 * k1_pe1 * log(exp((L - Lr1_pe1)/k1_pe1)+1) + eta*V;
        
    end

    function Fpe2 = Fpe2(L)
        %---------------------------
        % passive element 2
        % input: normalized muscle length
        % output: passive element force (0-1)
        %---------------------------
        c2_pe2 = -0.02;
        k2_pe2 = -21;
        Lr2_pe2 = 0.70;
        
        Fpe2 = c2_pe2*exp((k2_pe2*(L-Lr2_pe2))-1);
        
    end

    function Fse = Fse(LT)
        %---------------------------
        % series elastic element (tendon)
        % input: tendon length
        % output: tendon force (0-1)
        %---------------------------
        cT_se = 27.8; %27.8
        kT_se = 0.0047;
        LrT_se = 0.964;
        
        Fse = cT_se * kT_se * log(exp((LT - LrT_se)/kT_se)+1);
        
    end

    function dxdt = contractionDynamics(t,x)
       dxdt = zeros(2,1);
       dxdt(1) = x(2);
       dxdt(2) = (ForceSE*cos(alpha)-ForceCE*(cos(alpha)).^2)/(mass) ...
            + (x(2).^2*tan(alpha).^2)/(x(1));
    end

    function [Lce_initial,Lse_initial,Lmax] =  InitialLength(modelParameter)
        %---------------------------
        % Determine the initial lengths of muscle and tendon and maximal
        % muscle length
        %---------------------------
        
        % serires elastic element parameters
        cT = 27.8;
        kT = 0.0047;
        LrT = 0.964;
        % parallel passive element parameters
        c1 = 23;
        k1 = 0.046;
        Lr1 = 1.17;
        
        % passive force produced by parallel passive element at maximal
        % muscle length
        PassiveForce = c1 * k1 * log(exp((1 - Lr1)/k1)+1);
        % tendon length at the above passive force
        Normalized_SE_Length = kT*log(exp(PassiveForce/cT/kT)-1)+LrT;
        
        % maximal musculotendon length defined by joint range of motion
        Lmt_temp_max = modelParameter.optimalLength*cos(modelParameter.pennationAngle) ...
            +modelParameter.tendonSlackLength + 1;
        
        % optimal muscle length
        L0_temp = modelParameter.optimalLength;
        % optimal tendon length (Song et al. 2008)
        L0T_temp = modelParameter.tendonSlackLength*1.05;
        
        % tendon length at maximal muscle length
        SE_Length =  L0T_temp * Normalized_SE_Length;
        % maximal fasicle length
        FasclMax = (Lmt_temp_max - SE_Length)/L0_temp;
        % maximal muscle fiber length
        Lmax = FasclMax/cos(modelParameter.pennationAngle);
        
        % initial musculotendon length defined by the user input
        Lmt_temp = modelParameter.muscleInitialLength * cos(modelParameter.pennationAngle) + modelParameter.tendonInitialLength;
        
        % initial muscle length determined by passive muscle force and
        % tendon force
        InitialLength =  (Lmt_temp-(-L0T_temp*(kT/k1*Lr1-LrT-kT*log(c1/cT*k1/kT))))/(100*(1+kT/k1*L0T_temp/Lmax*1/L0_temp)*cos(modelParameter.pennationAngle));
        % normalize the muscle legnth to optimal muscle length
        Lce_initial = InitialLength/(L0_temp/100);
        % calculate initial length of tendon and normalize it to optimal
        % tendon length
        Lse_initial = (Lmt_temp - InitialLength*cos(modelParameter.pennationAngle)*100)/L0T_temp;
    end

end
