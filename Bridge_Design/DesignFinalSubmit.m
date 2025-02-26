%% 0. Initialize Parameters
% 

L = 1200; % Length of bridge
n = 1200; % Discretize into 1 mm seg.
P = 135; % Weight of one freight car [N]
x = linspace(0, L, n+1); % x-axis
%% 1. SFD, BMD under train loading
% 

x_train = [52 228 392 568 732 908]; % Train Load Locations
x_train = x_train - x_train(end); % train starts from the left
L_train = x_train(end) - x_train(1);
P_train = [1 1 1 1 1.35 1.35] * P/2;
n_train = 600; % num of train locations
SFDi = zeros(n_train, n+1); % 1 SFD for each train loc.
BMDi = zeros(n_train, n+1); % 1 BMD for each train loc.

% Solve for SFD and BMD with the train at different locations

for i = 1:n_train
    train_pos = x_train + (i-1)*(L_train + L)/n_train;

    % sum of moments at A eqn
    real_train_pos = train_pos(train_pos >= 0 & train_pos <= L); 
    real_P = P_train(train_pos >= 0 & train_pos <= L);
    moment_A = sum(real_train_pos.*real_P);

    % sum of Fy eqn
    By = moment_A / L;
    Ay = sum(real_P) - By;

    SFDi(i, 1) = Ay;
    next_ind = 1;
    for xi = 2:n
        SFDi(i, xi) = SFDi(i, xi-1);
        if next_ind <= length(real_P) & xi >= real_train_pos(next_ind)
            SFDi(i, xi) = SFDi(i, xi) - real_P(next_ind);
            next_ind = next_ind + 1;
        end
    end
    SFDi(i, end) = SFDi(i, end-1) + By;
    BMDi(i,:) = cumsum(-SFDi(i,:).*(L/n));
end

for i = 1:n_train
    plot(x, SFDi(i,:), DisplayName="n\_train="+i)
    hold on
end
title("SFDi")
xlabel("Bridge Location (mm)")
ylabel("Shear Sorce (N)")
hold off

for i = 1:n_train
    plot(x, BMDi(i,:), DisplayName="n\_train="+i)
    hold on
end
title("BMDi")
xlabel("Bridge Location (mm)")
ylabel("Bending Moment (Nmm)")
hold off

SFD = max(abs(SFDi)) % SFD envelope
BMD = min(BMDi) % BMD envelope
min(BMD)
plot(x, SFD)
title("Shear Force Envelope")
xlabel("Bridge Location")
ylabel("Max Shear Force (N)")

plot(x, BMD)
title("Bending Moment Envelope")
xlabel("Bridge Location")
ylabel("Max Bending Moment (Nmm)")
%% 2. Define Bridge Parameters

% = x, b, h

param_tf = [0, 100, 1.27];
param_glue = [0, 80/2, 1.27];
param_side = [0, 1.27, 140-1.27];
param_bot = [0, 80, 1.27];
param_inside = [0, 0, 0;
                400, 80, 1.27;
                800, 0, 0];
param_insideg = [0, 0, 0;
                400, 1.27, 15;
                800, 0, 0];
param_insbot = [0, 0, 0;
                575, 50, 1.27;
                625, 0, 0];

% locations of diaphragms
param_diaph = [0, 50, 250, 950, 1150, 1200];

btf = zeros(1, n+1);
ttf = zeros(1, n+1);
bglue = zeros(1, n+1);
hglue = zeros(1, n+1);
bside = zeros(1, n+1);
hside = zeros(1, n+1);
bbot = zeros(1, n+1);
hbot = zeros(1, n+1);
diaph = zeros(1, n+1);
binside = zeros(1, n+1);
hinside = zeros(1, n+1);
binsideg = zeros(1, n+1);
hinsideg = zeros(1, n+1);
binsbot = zeros(1, n+1);
hinsbot = zeros(1, n+1);

btf(1,1) = param_tf(1,2);
ttf(1,1) = param_tf(1,3);
bglue(1,1) = param_glue(1,2);
hglue(1,1) = param_glue(1,3);
bside(1,1) = param_side(1,2);
hside(1,1) = param_side(1,3);
bbot(1,1) = param_bot(1,2);
hbot(1,1) = param_bot(1,3);
diaph(1,1) = param_diaph(1,2);
binside(1,1) = param_inside(1,2);
hinside(1,1) = param_inside(1,3);
binsideg(1,1) = param_insideg(1,2);
hinsideg(1,1) = param_insideg(1,3);
binsbot(1,1) = param_insbot(1,2);
hinsbot(1,1) = param_insbot(1,3);

next_ind_tf = 2;
next_ind_glue = 2;
next_ind_side = 2;
next_ind_bot = 2;
next_ind_diaph = 2;
next_ind_inside = 2;
next_ind_insbot = 2;

for i = 2 : n+1
    btf(1, i) = btf(1, i-1);
    ttf(1, i) = ttf(1, i-1);
    bglue(1, i) = bglue(1, i-1);
    hglue(1, i) = hglue(1, i-1);
    bside(1, i) = bside(1, i-1);
    hside(1, i) = hside(1, i-1);
    bbot(1, i) = bbot(1, i-1);
    hbot(1, i) = hbot(1, i-1);
    diaph(1, i) = diaph(1, i-1);
    binside(1, i) = binside(1, i-1);
    hinside(1, i) = hinside(1, i-1);
    binsideg(1, i) = binsideg(1, i-1);
    hinsideg(1, i) = hinsideg(1, i-1);
    binsbot(1, i) = binsbot(1, i-1);
    hinsbot(1, i) = hinsbot(1, i-1);

    if next_ind_tf <= length(param_tf(:,1)) & i >= param_tf(next_ind_tf, 1)
        btf(1, i) = param_tf(next_ind_tf, 2);
        ttf(1, i) = param_tf(next_ind_tf, 3);
        next_ind_tf = next_ind_tf + 1;
    end
    if next_ind_glue <= length(param_glue(:,1)) & i >= param_glue(next_ind_glue, 1)
        bglue(1, i) = param_glue(next_ind_glue, 2);
        hglue(1, i) = param_glue(next_ind_glue, 3);
        next_ind_glue = next_ind_glue + 1;
    end
    if next_ind_side <= length(param_side(:,1)) & i >= param_side(next_ind_side, 1)
        bside(1, i) = param_side(next_ind_side, 2);
        hside(1, i) = param_side(next_ind_side, 3);
        next_ind_side = next_ind_side + 1;
    end
    if next_ind_bot <= length(param_bot(:,1)) & i >= param_bot(next_ind_bot, 1)
        bbot(1, i) = param_bot(next_ind_bot, 2);
        hbot(1, i) = param_bot(next_ind_bot, 3);
        next_ind_bot = next_ind_bot + 1;
    end
    if next_ind_diaph < length(param_diaph(1,:)) & i >= param_diaph(1, next_ind_diaph)
        next_ind_diaph = next_ind_diaph + 1;
        diaph(1, i) = param_diaph(1,next_ind_diaph)-param_diaph(1, next_ind_diaph-1);
    end
    if next_ind_bot <= length(param_bot(:,1)) & i >= param_bot(next_ind_bot, 1)
        bbot(1, i) = param_bot(next_ind_bot, 2);
        hbot(1, i) = param_bot(next_ind_bot, 3);
        next_ind_bot = next_ind_bot + 1;
    end
    if next_ind_inside <= length(param_inside(:,1)) & i >= param_inside(next_ind_inside, 1)
        binside(1, i) = param_inside(next_ind_inside, 2);
        hinside(1, i) = param_inside(next_ind_inside, 3);
        binsideg(1, i) = param_insideg(next_ind_inside, 2);
        hinsideg(1, i) = param_insideg(next_ind_inside, 3);
        next_ind_inside = next_ind_inside + 1;
    end
    if next_ind_insbot <= length(param_insbot(:,1)) & i >= param_insbot(next_ind_insbot, 1)
        binsbot(1, i) = param_insbot(next_ind_insbot, 2);
        hinsbot(1, i) = param_insbot(next_ind_insbot, 3);
        next_ind_insbot = next_ind_insbot + 1;
    end
   
end

%% 3. Calculate Sectional Properties
% 

% ybar. location of centroidal axis from the bottom

Aside = bside .* hside;
yside = hside ./ 2;
Iside = bside .* (hside.^3) ./ 12;

Abot = bbot .* hbot;
ybot = hbot ./ 2;
Ibot = bbot .* (hbot.^3) ./ 12;

Atf = btf .* ttf;
ytf = hside + ttf./2;
Itf = btf .* (ttf.^3) ./ 12;

Aglue = bglue .* hglue;
yglue = hside - hglue./2;
Iglue = bglue .* (hglue.^3) ./ 12;

Ainside = binside .*hinside;
yinside = hside - hglue - hinside./2;
Iinside = binside .* (hinside.^3) ./12;

Ainsideg = binsideg .* hinsideg;
yinsideg = hside-hglue-hinside-hinsideg./2;
Iinsideg = binsideg .* hinsideg.^3 ./12;

Ainsbot = binsbot .* hinsbot;
yinsbot = hbot + hinsbot./2;
Iinsbot = binsbot .* hinsbot.^3./12;

Atotal = Atf + 2.*Aside + 2.*Aglue + Ainside + 2.*Ainsideg + Abot + Ainsbot;
numerator = Atf.*ytf + 2.*Aside.*yside + 2.*Aglue.*yglue + Abot.*ybot + ...
     + Ainside.*yinside + 2.*Ainsideg.*yinsideg + Ainsbot.*yinsbot;
ybar = numerator./Atotal

ybartotop = ttf + hside - ybar

I = 2.*Iside + 2.*Iglue + Itf + Ibot + Iinside + 2.*Iinsideg + Iinsbot...
    + 2.*Aside.*(yside-ybar).^2 + 2.*Aglue.*(yglue-ybar).^2 + Atf.*(ytf-ybar).^2 + Abot.*(ybot-ybar).^2 ...
    + Ainside.*(yinside-ybar).*2 + 2.*Ainsideg.*(yinsideg-ybar).^2 + Ainsbot.*(yinsbot-ybar).^2

% Q at centroidal axes
Qcent = 2.*ybar.*bside.*ybar./2 + Abot.*(ybar-hbot./2) + Ainsbot.*(ybar-hbot-hinsbot./2)

% Q at glue location
Qgluetop = Atf.*(ybartotop - ttf./2)
Qglueinside = Atf.*(ybartotop - ttf./2) + 2.*Aglue.*(ybartotop-ttf-hglue./2) + 2*1.27*1.27.*(ybartotop-1.27*1.5)
Qgluebot = 1.27*(80+2*1.27) .* (ybar-1.27/2)
%% 4. Calculate Applied Stress

% S: sigma (flexural stress)   T: tau (shear stress)

S_top = abs(BMD).*(ybartotop)./I
S_bot = abs(BMD).*ybar./I
T_cent = abs(SFD).*Qcent./(I.*(2.*bside))
T_gluetop = abs(SFD).*Qgluetop./(I.*(2.*(bglue+1.27)))
T_glueinside = abs(SFD(binside~=0)).*Qglueinside(binside~=0)./(I(binside~=0).*(binside(binside~=0)))
T_gluebot = abs(SFD(binsbot~=0)).*Qgluebot(binsbot~=0)./(I(binsbot~=0).*(binsbot(binsbot~=0)))
%% 5. Material and Thin Plate Buckling Capacities
% 

E = 4000;
mu = 0.2;
S_tens = 30;
S_comp = 6;
T_max = 4;
T_gmax = 2;
S_buck1 = 4*pi^2*E/(12*(1-mu^2)) * (1.27*3/80).^2
S_buck2 = 0.425*pi^2*E/(12*(1-mu^2)) * (ttf./((btf-bbot-2.*bside)./2)).^2
S_buck3 = 6*pi^2*E/(12*(1-mu^2)) * (bside./(hside-ybar)).^2
T_buck = 5*pi^2*E/(12*(1-mu^2)) * ((bside./(hside-hbot-hglue-hinside)).^2 + (bside./diaph).^2)
%% 6. FOS
% 

FOS_tens = S_tens ./ S_bot
FOS_comp = S_comp ./ S_top
FOS_shear = T_max ./ T_cent
FOS_glue_top = T_gmax ./ T_gluetop
FOS_glue_ins = T_gmax ./ T_glueinside
FOS_glue_bot = T_gmax ./ T_gluebot
FOS_buck1 = S_buck1 ./ S_top
FOS_buck2 = S_buck2 ./ S_top
FOS_buck3 = S_buck3 ./ S_top
FOS_buckV = T_buck ./ T_cent
 
FOS_tensmin = min(FOS_tens)
FOS_compmin = min(FOS_comp)
FOS_shearmin = min(FOS_shear)
FOS_gluemin = min([min(FOS_glue_top), min(FOS_glue_ins), min(FOS_glue_bot)])
FOS_buck1min = min(FOS_buck1)
FOS_buck2min = min(FOS_buck2)
FOS_buck3min = min(FOS_buck3)
FOS_buckVmin = min(FOS_buckV)

%% 7. Min FOS and the failure load Pfail

minFOS = min([FOS_tens, FOS_comp, FOS_shear, FOS_buck1, FOS_buck2, FOS_buck3, FOS_buckV])
Pf = minFOS*P

%% 8. Vfail and Mfail
% 

Mf_tens = FOS_tens.*BMD
Mf_comp = FOS_comp.*BMD
Vf_shear = FOS_shear.*SFD
% Vf_glue = FOS_glue.*SFD
Mf_buck1 = FOS_buck1.*BMD
Mf_buck2 = FOS_buck2.*BMD
Mf_buck3 = FOS_buck3.*BMD
Vf_buckV = FOS_buckV.*SFD
%% 9. Output plots of Vfail and Mfail
% 

subplot(2,3,1)
hold on; grid on; grid minor;
plot(x, Vf_shear, 'r')
plot(x, -Vf_shear, 'r')
plot(x, SFDi, 'k');
plot([0, L], [0, 0], 'k', 'LineWidth', 2)
legend('Matboard Shear Failure')
xlabel('Distance along bridge (mm)')
ylabel('Shear Force (N)')

% subplot(2,3,2)
% hold on; grid on; grid minor;
% plot(x, Vf_glue, 'r')
% plot(x, -Vf_glue, 'r')
% plot(x, SFDi, 'k');
% plot([0, L], [0, 0], 'k', 'LineWidth', 2)
% legend('Glue Shear Failure')
% xlabel('Distance along bridge (mm)')
% ylabel('Shear Force (N)')

subplot(2,3,3)
hold on; grid on; grid minor;
plot(x, Vf_buckV, 'r')
plot(x, -Vf_buckV, 'r')
plot(x, SFDi, 'k');
plot([0, L], [0, 0], 'k', 'LineWidth', 2)
legend('Buckling Shear Failure')
xlabel('Distance along bridge (mm)')
ylabel('Shear Force (N)')

subplot(2,3,4)
hold on; grid on; grid minor;
plot(x, Mf_tens, 'r')
plot(x, Mf_comp, 'b')
plot(x, BMDi, 'k');
plot([0, L], [0, 0], 'k', 'LineWidth', 2)
legend('Matboard Tension Failure', ...
    'Matboard Compression Failure')
xlabel('Distance along bridge (mm)')
ylabel('Bending Moment (Nmm)')

subplot(2,3,5)
hold on; grid on; grid minor;
plot(x, Mf_buck1, 'r')
plot(x, Mf_buck2, 'b')
plot(x, BMDi, 'k');
plot([0, L], [0, 0], 'k', 'LineWidth', 2)
legend('Matboard Buckling Failure, Top Flange - Mid', ...
    'Matboard Buckling Failure, Top Flange - Sides')
xlabel('Distance along bridge (mm)')
ylabel('Bending Moment (Nmm)')

subplot(2,3,6)
hold on; grid on; grid minor;
plot(x, Mf_buck3, 'r')
plot(x, BMDi, 'k');
plot([0, L], [0, 0], 'k', 'LineWidth', 2)
legend('Matboard Buckling Failure, Webs')
xlabel('Distance along bridge (mm)')
ylabel('Bending Moment (Nmm)')

xlim([0 1200])
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96])
