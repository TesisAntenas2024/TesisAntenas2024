% Antena de microstrip con inset fed
% MATLAB(R) y Antenna Toolbox 23.2.
% Uso, ejemplo:
%   f = 2450e6 % frecuencia de resonancia
%   er = 4.76  % permitividad del sustrato
%   h = 1.6e-3 % altura del sustrato
% Sebastian Fraga, Lucas Reyes
% Con ayuda del código de Tom Schucker https://github.com/Tschucker/matlab_antenna_designs/blob/master/patch_antennas/oshpark_2_45ghz_patch/oshpark_2_45ghz_patch.m

clear;

%% Parámetros de diseño
fc = 2450e6; 
er = 4.76;
h =  1.6e-3;

% Ancho de la antena patch
Wp = 0.0373;

% Longitud de la antena patch
Lp = 0.0289;

% Ancho de microstrip para 50 ohmios
W = 2.3e-3; %m

% Ancho del inset
g = W*3

% Longitud del inset
d = 0.0089; %0.0100 -> 2400 MHz 0.0098 -> 2450 MHz

% Longitud del plano de tierra
%Lg = Lp + 6*h;
Lg = Lp + 12*h;

% Ancho del plano de tierra
Wg = Wp + 6*h;
%Wg = Wp + 12*h;

% Creo la geometría de la antena
ground = antenna.Rectangle('Length', Lg, 'Width', Wg);
patch = antenna.Rectangle('Length', Lp, 'Width', Wp);
notch = antenna.Rectangle('Length', d, 'Width',g, 'Center', [(Lp/2)-(d/2),0]);
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
figure;
mesh(pcbAnt, 'MaxEdgeLength',100.0e-3,'MinEdgeLength',0.8e-3);

%% Graficar la impedancia de la antena patch básica.
% Enumerar frecuencias
freqRange = linspace(fc-0.05*fc,fc + 0.1*fc,22);

%% Graficar impedancia compleja
figure
impedance(pcbAnt,freqRange)

%% Graficar parámetros S de RF
S = sparameters(pcbAnt, freqRange);
figure; 
rfplot(S);

%% Patrón de radiación 2D
figure;
patternAzimuth(pcbAnt, fc, 0, 'Azimuth', 0:5:360)

%% Graficar densidad de corriente
figure;
current(pcbAnt, fc)

%% Graficar el patrón de radiación de la antena patch básica.
figure
pattern(pcbAnt, fc)

%% Generación de GERBER de la PCB

% Conector
c = PCBConnectors.SMAEdge;
c.SignalLineWidth = W;
c.EdgeLocation = 'east';
c.ExtendBoardProfile = false;

% Servicio de PCB
%s = PCBServices.OSHParkWriter;
s = PCBServices.MayhewWriter;
service.Filename = 'Patch_Antenna.zip';

% Escribir GERBER
PW = PCBWriter(pcbAnt,s,c);
gerberWrite(PW);
