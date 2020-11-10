clear all, close all, clc;

N = 8;    % broj korisnika
ro = 1;   % koef. zaobljenja (za 0 se dobija prav. impuls)
Peb_zadato = 0.025;   % zadata ver. greske

Bitski_Protok(1:N,1) =     [64    64    64      64      64    64    64   64   ];
Bitski_Protok'

disp('1. SLUCAJ : idealna kontrola snage, uz ostale neidealnosti')

tau(1:N,1) =               [0     1/4   1/4     1/4     1/2   0     0    0    ];    % Greska sinhronizacije PSS [chip] po korisnicima
theta(1:N,1) =             [0     0     pi/16   pi/8    0     0     0    0    ];    % Greska sinhronizacije faze [rad] po korisnicima 
RaspSnage(1:N,1) =         1*ones(1,8);    % Raspodela snage                             
SIRvektor(1:N,1) =         [Inf   Inf   Inf     Inf     Inf   -12   -6   0    ];    % SIR [dB] po korisnicima  

% Prvi korisnik: potpuno idealan i neometan slucaj
% Drugi-peti korisnik: greska sinhronizacije po postavci projekta
% Sesti-osmi korisnik: uneta interferencija, SIR je -12 dB, -6 dB ili 0 dB

calc_function(ro,Peb_zadato,tau,theta,RaspSnage,SIRvektor,Bitski_Protok,N);



disp('2. SLUCAJ : neidealna kontrola snage')
disp('a) jedan korisnik ima 2x vecu snagu')
tau(1:N,1) =       [0     0     0       0       0     0     0    0    ];     
theta(1:N,1) =     [0     0     0       0       0     0     0    0    ];    
RaspSnage(1:N,1) =         1*ones(1,8); 
RaspSnage(1,1) = 2;                         
SIRvektor(1:N,1) = [Inf   Inf   Inf     Inf     Inf   -12   -6   0    ];

calc_function(ro,Peb_zadato,tau,theta,RaspSnage,SIRvektor,Bitski_Protok,N);    



disp('b) jedan korisnik ima 2x manju snagu')
tau(1:N,1) =       [0     0     0       0       0     0     0    0    ];     
theta(1:N,1) =     [0     0     0       0       0     0     0    0    ];    
RaspSnage(1:N,1) =         1*ones(1,8);      
RaspSnage(1,1) = 0.5;                          
SIRvektor(1:N,1) = [Inf   Inf   Inf     Inf     Inf   -12   -6   0    ];  

calc_function(ro,Peb_zadato,tau,theta,RaspSnage,SIRvektor,Bitski_Protok,N);