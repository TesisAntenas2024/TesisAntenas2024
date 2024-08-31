% Antena de microstrip con inset fed - Calculo y Corner analysis
% MATLAB(R) y Antenna Toolbox 23.2.
% Uso, ejemplo:
%   f = 2.4e9 // frecuencia de resonancia
%   er = 2.2  // permitividad del sustrato
%   h = 1.5e-3 // altura del sustrato
% Sebastian Fraga, Lucas Reyes

clear;

%Inpunts
disp('PARAMETROS DE ENTRADA:')
f = 2.4e9
er = 2.2
h =  1.5e-3

%Constantes
c = 3.0e8;
Zo = 50;
lambda0 = c/f;
k0 = 2*pi/lambda0;
setWavenumber(k0);
% Verificar que la altura del sustrato sea lo suficientemente grande
% para que el efecto del los fringing fields no sea significante
% Consultar en la hoja de datos del sustrato los tamaños estandar

%Calculo para parche
[Wp, Lp, ereff] = geometry_Patch(f, er, h, c);
setGlobalAncho(Wp);
setGlobalLargo(Lp);
%Calculo para microstrip
[Za, W, g_1, g_2, d_1, d_2, Lf] = geometry_Notch(er, h, Zo, c, f);
[dr, ref_cof,VSWR,return_loss,Antenna_Bandwidth, Efficiency] = ideal_Parameters(Zo, Za, Wp, Lp, c, er, h);
%Realizar la antena en PCB
%corner_Analysis(Wp, Lp, W, g_1, g_2, d_1, d_2, f, er, h, Zo);

% RES: Calcula los parametros de ancho y largo de la antena
function [W, L, erff] = geometry_Patch(fr, err, hr, c)

    % Calcular el ancho del parche
    W = (c/(2*fr))*(sqrt(2/(err+1)));
    % Calcular el valor de la permitividad efectiva dada la ecuacion:
    erff = 0.5*(err+1) + 0.5*(err-1)*(1+12*hr/W)^(-0.5);
    % Calcular el largo efectivo
    Leff = c/(2*fr*sqrt(erff));
    % Calcular la extension normalizada del largo
    DL = hr*0.412*((erff+0.3)*(W/hr + 0.264))/((erff-0.258)*(W/hr + 0.8));
    % Calcular el largo del parche
    L = Leff - 2*DL;
end

%RES: Calcula el largo en funcion del 

% RES: Calcula los parametros de ancho y largo del notch 
% (apertura en la anten de parche por la microstrip)
function [Rin, W, gr1, gr2, Lf1, Lf2, L] = geometry_Notch(err, hr, Zo, c, fr)
    
    % Calcular el ancho de la microstrip
    Wl = claculate_Roots(Zo, err);
    W = double(Wl)*hr;
    erff = 0.5*(err+1) + 0.5*(err-1)*(1+12*hr/W)^(-0.5);
    % Calculo el ancho del notch
    t = 0:pi;
    g1(t);
    G1 = integral(@g1,0,pi);
    g12(t);
    G12 = integral(@g12,0,pi);
    Rin = 1/(2*(G1+G12));
    % Se estudian varios valores para ver parametros como el return loss,
    % frecuencia de resonancia, etc.
    Lf1 = (getGlobalLargo/pi)*(acos(sqrt(Zo/Rin)));
    Lf2 = (10^-4)*(0.001699*err^7 + 0.13671*err^6 - 6.1783*err^5 + 93.187*err^4 - 682.69*err^3 + 2561.9*err^2 - 4031*err + 6697)*getGlobalLargo/2;
    DL = hr*0.412*((erff+0.3)*(W/hr + 0.264))/((erff-0.258)*(W/hr + 0.8));
    Leff = c/(2*fr*sqrt(erff));
    L = Leff - 2*DL;
    gr1 = 3*W; 
    gr2 = 1.5*W;
end

% RES: Devuelve los parametros ideales de la antena. Estos mismos, despues
% se veran afectados por la variacion de la frecuencia de resonancia y las
% dimensiones del notch
function [dr, rho,ROE,RL,B, E] = ideal_Parameters(Zo, Rin, Wp, Lp, c0, er, h)

    dr = (Lp/pi)*(asin(sqrt(Zo/Rin)));
    rho = (Rin-50)/(Rin+50);
    ROE = (1+abs(rho))/(1-abs(rho));
    RL = -20*log10(rho);
    B = 0.91*(er-1)/er^3*Wp*h*c0/Lp^3/1e6;
    lambda0 = 2*pi/getWavenumber;
    Ratio = 3/32*sqrt(er)*lambda0^2*0.0009/(h*Wp);
    E = 1/(1 + Ratio);
end

% RES: Dibuja la antena en la PCB y realiza los calculos pertinantes.
% Quitar comenrarios de las funciones que se quieran probar
function corner_Analysis(Wp, Lp, W, g_1, g_2, d_1, d_2, f, er, h, Zo)

    Wg = Wp + 6*h;
    Lg = Lp + 6*h;

    % Wg = Wp + 12*h;
    % Lg = Lp + 12*h;

    anchos = [g_1, g_2, (g_1 + g_2)/2];
    largos = [d_1, d_2, (d_1 + d_2)/2];

    i=1;
    while i<3
     k=1;
        while k<3
    
        %Creo la geometria de la antena    
        ground = antenna.Rectangle('Length', Lg, 'Width', Wg);
        patch = antenna.Rectangle('Length', Lp, 'Width', Wp);
        notch = antenna.Rectangle('Length', largos(i), 'Width', anchos(k), 'Center', [(Lp/2)-(largos(i)/2),0]);
        feed = antenna.Rectangle('Length', Lg/2, 'Width', W, 'Center', [(Lg/4),0]);
        substrait = dielectric('Name','Custom','EpsilonR', er, 'Thickness', h); %dielectric('Rt-Duroid');
        %conductor = metal('Name','Copper','Conductivity', 5.8001e07, 'Thickness', 1e-3); %dielectric('Rt-Duroid');
        
        geometry = (patch-notch) + feed;
    
        % Defino los materiales en la PCB
        pcbAnt = pcbStack;
        pcbAnt.Name = 'Inset-Fed Patch Antenna';
        pcbAnt.BoardThickness = h;
        pcbAnt.BoardShape = ground;
        pcbAnt.Layers = {geometry,substrait,ground};
        pcbAnt.FeedLocations = [Lg/2 0 1 3];
        pcbAnt.FeedDiameter = W/2;
        % Dibujar PCB
        figure;
        show(pcbAnt);
        mesh(pcbAnt, 'MaxEdgeLength',2.5e-3,'MinEdgeLength',0.8e-3);

        % Analisis de datos 
        freqRange = linspace(f-0.1*f,f + 0.1*f,20);
        % impedancia
        figure;
        impedance(pcbAnt, freqRange)
        % Return Loss / VSWR
        %figure;
        %returnLoss (pcbAnt, freqRange);
        
        k = k+1;
        end
        i = i+1;
    end

    %Creo la geometria de la antena
    ground = antenna.Rectangle('Length', Lg, 'Width', Wg);
    patch = antenna.Rectangle('Length', Lp, 'Width', Wp);
    notch = antenna.Rectangle('Length', largos(3), 'Width', anchos(3), 'Center', [(Lp/2)-(largos(3)/2),0]);
    feed = antenna.Rectangle('Length', Lg/2, 'Width', W, 'Center', [(Lg/4),0]);
    substrait = dielectric('Name','Custom','EpsilonR', er, 'Thickness', h); %dielectric('Rt-Duroid');
    %conductor = metal('Name','Copper','Conductivity', 5.8001e07, 'Thickness', 1e-3); %dielectric('Rt-Duroid');
        
    geometry = (patch-notch) + feed;
    
    % Defino los materiales en la PCB
    pcbAnt = pcbStack;
    pcbAnt.Name = 'Inset-Fed Patch Antenna';
    pcbAnt.BoardThickness = h;
    pcbAnt.BoardShape = ground;
    pcbAnt.Layers = {geometry,substrait,ground};
    pcbAnt.FeedLocations = [Lg/2 0 1 3];
    pcbAnt.FeedDiameter = W/2;
    % Dibujar PCB
    figure;
    show(pcbAnt);
    mesh(pcbAnt, 'MaxEdgeLength',2.5e-3,'MinEdgeLength',0.8e-3);
    
    % Analisis de datos 
    freqRange = linspace(f-0.1*f,f + 0.1*f,20);
    % impedancia
    figure;
    impedance(pcbAnt, freqRange)
    % Return Loss / VSWR
    %figure;
    %returnLoss (pcbAnt, freqRange);


    % Resize
    %scale = 2.4/2.3;
    %antennaObject.Length = antennaObject.Length/(scale);
    %antennaObject.Width = antennaObject.Width/(scale);
    %antennaObject.Height = antennaObject.Height/(scale);
    %antennaObject.PatchCenterOffset = antennaObject.PatchCenterOffset/(scale);
    %antennaObject.StripLineWidth = antennaObject.StripLineWidth/(scale);
    %antennaObject.NotchLength = antennaObject.NotchLength/(scale);
    %antennaObject.NotchWidth = antennaObject.NotchWidth/(scale);
   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Resoluciones Numericas %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% RES: Calcula las raices (W/h) para una impedancia dada
function [e1] = claculate_Roots(Ro, erff)

    K = (120*pi)/(Ro*sqrt(erff)) - 1.393;
    syms x;
    a = x + (2/3)*log(x+1.4444) - K;
    e1 = solve(a);
end

% RES: Calcula la funcion de la conductacia G1
function [f] = g1(t)
    
    k0 = getWavenumber;
    W = getGlobalAncho;
    f = (1/(120*pi^2))*((sin((k0*W*0.5*cos(t))/cos(t))).^2*(sin(t)).^3);
end

% RES: Calcula la funcion de la conductacia G12
function [k] = g12(t)

    k0 = getWavenumber;
    W = getGlobalAncho;
    L = getGlobalLargo;
    k =(1/(120*pi^2))*((((sin(k0*W*0.5*cos(t)))/cos(t)).^2)*(sin(t).^3)).*besselj(0,k0*L*sin(t));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Variables Globales %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
function setGlobalAncho(val)
global ancho
ancho = val;
end

function r = getGlobalAncho
global ancho
r = ancho;
end

function setGlobalLargo(val)
global largo
largo = val;
end

function r = getGlobalLargo
global largo
r = largo;
end

function setWavenumber(val)
global k0
k0 = val;
end

function r = getWavenumber
global k0
r = k0;
end