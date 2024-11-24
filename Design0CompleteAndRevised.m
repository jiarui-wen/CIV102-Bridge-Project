%% 0. Initialize Parameters
% 

L = 1200; % Length of bridge
n = 1200; % Discretize into 1 mm seg.
P = 400; % Total weight of train [N]
x = linspace(0, L, n+1); % x-axis
%% 1. SFD, BMD under train loading
% 

x_train = [52 228 392 568 732 908]; % Train Load Locations
x_train = x_train + 120;  % centre the train on the bridge
P_train = [1 1 1 1 1 1] * P/6;
n_train = 600; % num of train locations
SFDi = zeros(n_train, n+1); % 1 SFD for each train loc.
BMDi = zeros(n_train, n+1); % 1 BMD for each train loc.
% Solve for SFD and BMD with the train at different locations

for i = 1:n_train
    % % start location of train
    % train_pos = x_train - (n_train - i) * L / n_train
    % train_pos = x_train - (i-1)*L/(2*n_train)
    i;
    train_pos = x_train - (i-1)*x_train(end)/n_train;

    % sum of moments at A eqn
    real_train_pos = train_pos(train_pos >= 0 & train_pos <= L);
    real_P = P_train(train_pos >= 0 & train_pos <= L);
    moment_A = sum(real_train_pos.*real_P);

    % sum of Fy eqn
    By = moment_A / L;
    Ay = sum(real_P) - By;

    % construct applied loads
    % w(x)
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
% legend()
title("SFDi")
xlabel("Bridge location (mm)")
ylabel("Shear force (N)")
hold off

for i = 1:n_train
    plot(x, BMDi(i,:), DisplayName="n\_train="+i)
    hold on
end
% legend()
title("BMDi")
xlabel("Position of front of train (mm)")
ylabel("Bending moment (Nmm)")
hold off


SFD = max(abs(SFDi)) % SFD envelope
BMD = min(BMDi) % BMD envelope

plot(x, SFD)
title("SFE")
xlabel("Bridge location")
ylabel("Max shear force (N)")

plot(x, BMD)
title("BME")
xlabel("Bridge location")
ylabel("Max bending moment (Nmm)")
%% 2. Define Bridge Parameters

% = xc, bft, tft,
% param_tf = [0, 100, 1.27;
%         400, 100, 2.54;
%         800, 100, 1.27];  this is me experimenting with changing
%         cross-sections; for now we don't need this but the code is
%         written to work with changing cross-sections as well

% = x_c, b, h
param_tf = [0, 100, 1.27];
param_glue = [0, 5, 1.27];
param_side = [0, 1.27, 75];
param_bot = [0, 80-1.27*2, 1.27];

% param_tf = [0, 110, 1.27;
%             400, 110, 1.27;
%             800, 110, 1.27];
% param_glue = [0, 15, 1.27];
% param_side = [0, 1.27, 180-3.81];
% param_bot = [0, 100-1.27*2, 1.27];

%x_c Location, x, of cross-section change
%bft Top Flange Width
%tft Top Flange Thickness

% % Extracting user input assuming linear relationship
% bft = interp1(param(:,1), param(:,2), x)
% tft = interp1(param(:,1), param(:,3), x)

btf = zeros(1, n+1);
ttf = zeros(1, n+1);
bglue = zeros(1, n+1);
hglue = zeros(1, n+1);
bside = zeros(1, n+1);
hside = zeros(1, n+1);
bbot = zeros(1, n+1);
hbot = zeros(1, n+1);

btf(1,1) = param_tf(1,2);
ttf(1,1) = param_tf(1,3);
bglue(1,1) = param_glue(1,2);
hglue(1,1) = param_glue(1,3);
bside(1,1) = param_side(1,2);
hside(1,1) = param_side(1,3);
bbot(1,1) = param_bot(1,2);
hbot(1,1) = param_bot(1,3);

next_ind_tf = 2;
next_ind_glue = 2;
next_ind_side = 2;
next_ind_bot = 2;
for i = 2 : n+1
    btf(1, i) = btf(1, i-1);
    ttf(1, i) = ttf(1, i-1);
    bglue(1, i) = bglue(1, i-1);
    hglue(1, i) = hglue(1, i-1);
    bside(1, i) = bside(1, i-1);
    hside(1, i) = hside(1, i-1);
    bbot(1, i) = bbot(1, i-1);
    hbot(1, i) = hbot(1, i-1);

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
end

%% 3. Calculate Sectional Properties
% 

% ybar. location of centroidal axis from the bottom

Aside = bside .* hside
yside = hside ./ 2
Iside = bside .* (hside.^3) ./ 12

Abot = bbot .* hbot
ybot = hbot ./ 2
Ibot = bbot .* (hbot.^3) ./ 12

Atf = btf .* ttf
ytf = hside + ttf./2
Itf = btf .* (ttf.^3) ./ 12

Aglue = bglue .* hglue
yglue = hside - hglue./2
Iglue = bglue .* (hglue.^3) ./ 12

ybar = (2.*Aside.*yside + Abot.*ybot + Atf.*ytf + 2.*Aglue.*yglue) ./ (2.*Aside + Abot + Atf + 2.*Aglue)
ybartotop = ttf + hside - ybar

% I
I = 2.*Iside + 2.*Iglue + Itf + Ibot ...
    + 2.*Aside.*(yside-ybar).^2 + 2.*Aglue.*(yglue-ybar).^2 + Atf.*(ytf-ybar).^2 + Abot.*(ybot-ybar).^2

% Q at centroidal axes
Qcent = 2.*ybar.*bside.*ybar./2 + Abot.*(ybar-hbot./2)
% Q at glue location
Qglue = Atf.*(ybartotop - ttf./2)

%% 4. Calculate Applied Stress

% S: sigma (flexural stress)   T: tau (shear stress)

S_top = abs(BMD).*(ybartotop)./I
min(BMD)
ybartotop
S_bot = abs(BMD).*ybar./I
T_cent = abs(SFD).*Qcent./(I.*(2.*bside))
max(T_cent)
bglue = bglue+1.27
I
T_glue = abs(SFD).*Qglue./(I.*(2.*(bglue)))
%% 5. Material and Thin Plate Buckling Capacities
% 

E = 4000;
mu = 0.2;
S_tens = 30;
S_comp = 6;
T_max = 4;
T_gmax = 2;
S_buck1 = 4*pi^2*E/(12*(1-mu^2)) * (ttf./(bbot+2.*bside)).^2
S_buck2 = 0.425*pi^2*E/(12*(1-mu^2)) * (ttf./((btf-bbot-2.*bside)./2)).^2
S_buck3 = 6*pi^2*E/(12*(1-mu^2)) * (bside./(hside-ybar)).^2
T_buck = 5*pi^2*E/(12*(1-mu^2)) * ((bside./(hside-1.27*2)).^2 + (bside./400).^2)
%% 6. FOS
% 

FOS_tens = S_tens ./ S_bot
FOS_comp = S_comp ./ abs(S_top)
FOS_shear = T_max ./ T_cent
FOS_glue = T_gmax ./ T_glue
FOS_buck1 = S_buck1 ./ S_top
FOS_buck2 = S_buck2 ./ S_top
FOS_buck3 = S_buck3 ./ S_top
FOS_buckV = T_buck ./ T_cent
 
FOS_tensmin = min(FOS_tens)
FOS_compmin = min(FOS_comp)
FOS_shearmin = min(FOS_shear)
FOS_gluemin = min(FOS_glue)
FOS_buck1min = min(FOS_buck1)
FOS_buck2min = min(FOS_buck2)
FOS_buck3min = min(FOS_buck3)
FOS_buckVmin = min(FOS_buckV)

% FOS_tens = S_tens / max(S_bot)
% FOS_comp = S_comp / max(S_top)
% FOS_shear = T_max / max(T_cent)
% FOS_glue = T_gmax / max(T_glue)
% FOS_buck1 = S_comp / max(S_buck1)
% FOS_buck2 = S_comp / max(S_buck2)
% FOS_buck3 = S_comp / max(S_buck3)
% 
% FOS_buckV = min(min(T_max, min(T_buck))) / max(T_cent)
% 
% FOS_buckV = zeros(1, length(T_buck))
% for i = 1:length(T_buck)
%     if T_buck(1, i) < T_max
%         FOS_buckV(1, i) = T_buck(1, i)/T_cent(1, i);
%     else
%         FOS_buckV(1, i) = T_max/T_cent(1, i);
%     end
% end
% FOS_buckVmin = min(FOS_buckV)


%% 7. Min FOS and the failure load Pfail
% 

minFOS = min([FOS_tens, FOS_comp, FOS_shear, FOS_glue, FOS_buck1, FOS_buck2, FOS_buck3, FOS_buckV])
Pf = minFOS*P
%% 8. Vfail and Mfail
% 

Mf_tens = FOS_tens.*BMD
Mf_comp = FOS_comp.*BMD
Vf_shear = FOS_shear.*SFD
Vf_glue = FOS_glue.*SFD
Mf_buck1 = FOS_buck1.*BMD
Mf_buck2 = FOS_buck2.*BMD
Mf_buck3 = FOS_buck3.*BMD
Vf_buckV = FOS_buckV.*SFD

% Mf_tens = -S_tens.*I./ybar
% Mf_comp = -S_comp.*I./ybartotop
% Vf_shear = T_max.*I.*2.*bside ./ Qcent;
% Vf_glue = T_gmax.*I.*2.*bglue ./ Qglue;
% Mf_buck1 = -S_buck1.*I./ ybartotop;
% Mf_buck2 = -S_buck2.*I./ ybartotop;
% Mf_buck3 = -S_buck3.*I./ ybartotop; % not sure about ybartotop here
% Vf_buckV = min(T_max, T_buck).*I.*bside ./ Qcent
%% 9. Output plots of Vfail and Mfail
% 

subplot(2,3,1)
hold on; grid on; grid minor;
plot(x, Vf_shear, 'r')
plot(x, -Vf_shear, 'r')
plot(x, SFDi, 'k');
plot([0, L], [0, 0], 'k', 'LineWidth', 2)
% legend('Matboard Shear Failure')
xlabel('Distance along bridge (mm)')
ylabel('Shear Force (N)')

subplot(2,3,2)
hold on; grid on; grid minor;
plot(x, Vf_glue, 'r')
plot(x, -Vf_glue, 'r')
plot(x, SFDi, 'k');
plot([0, L], [0, 0], 'k', 'LineWidth', 2)
% legend('Glue Shear Failure')
xlabel('Distance along bridge (mm)')
ylabel('Shear Force (N)')

subplot(2,3,3)
hold on; grid on; grid minor;
plot(x, Vf_buckV, 'r')
plot(x, -Vf_buckV, 'r')
plot(x, SFDi, 'k');
plot([0, L], [0, 0], 'k', 'LineWidth', 2)
% legend('Buckling Shear Failure')
xlabel('Distance along bridge (mm)')
ylabel('Shear Force (N)')

subplot(2,3,4)
hold on; grid on; grid minor;
plot(x, Mf_tens, 'r')
plot(x, Mf_comp, 'b')
plot(x, BMDi, 'k');
plot([0, L], [0, 0], 'k', 'LineWidth', 2)
% legend('Matboard Tension Failure')
% legend('Matboard Compression Failure')
xlabel('Distance along bridge (mm)')
ylabel('Bending Moment (Nmm)')

subplot(2,3,5)
hold on; grid on; grid minor;
plot(x, Mf_buck1, 'r')
plot(x, Mf_buck2, 'b')
plot(x, BMDi, 'k');
plot([0, L], [0, 0], 'k', 'LineWidth', 2)
% legend('Matboard Buckling Failure, Top Flange - Mid')
% legend('Matboard Buckling Failure, Top Flange - Sides')
xlabel('Distance along bridge (mm)')
ylabel('Bending Moment (Nmm)')

subplot(2,3,6)
hold on; grid on; grid minor;
plot(x, Mf_buck3, 'r')
plot(x, BMDi, 'k');
plot([0, L], [0, 0], 'k', 'LineWidth', 2)
% legend('Matboard Buckling Failure, Webs')
xlabel('Distance along bridge (mm)')
ylabel('Bending Moment (Nmm)')