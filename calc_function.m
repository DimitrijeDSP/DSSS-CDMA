function calc_function(ro,Peb_zadato,tau,theta,RaspSnage,SIRvektor,Bitski_Protok,N)

pkg load communications;
pkg load signal;

SF = 16;  % procesno pojacanje
M = 4;    % broj simbola QPSK

Peb_izmereno = zeros(N,1);                                          

Vb = Bitski_Protok;   % Bitski protok po kanalu za svih 4 brzina signaliziranja
Vchip = 8 * 256000;   % konstantan za sve korisnike

disp('Iz Peb_zadato se za QPSK nalazi:')
Eb_pN_dB = 10*log10((erfcinv(2*Peb_zadato))^2)
N_bita = 100/Peb_zadato;
while 2 > 1
  if ceil(log2(N_bita)) == floor(log2(N_bita))
    break;
  endif
  N_bita = N_bita + 1;
endwhile
N_bita = round(N_bita)   % duzina koja treba da se uzme

%N_bita = 8;         % radi brze simulacije

tacke = [e^(i*pi/4) e^(i*3*pi/4) e^(i*5*pi/4) e^(i*7*pi/4)];

sim_info_duz = N_bita / log2(M);

N_chipova = SF * sim_info_duz;  % Broj chip-ova    

Gold_Sekvenca = zeros(N_chipova,1);
OVSF_Sekvenca = zeros(N_chipova,2);

% Generisanje OVSF sekvence
ovsf_seq = ovsf16();   % 16 OSVF sekvenci duzine 16

% Generisanje Gold-ove sekvence
gold_seq = gold(N_chipova);  % jedna gold-ova sekvenca duzine N_chipova
Gold_Sekvenca = gold_seq';

GreskaSinhFaze(1:N,1) = exp(-1i*theta(1:N,1)); 

for tr_kor = 1:N

    Matrica_Chipova = zeros(1,N_chipova,2);  

    Signal = zeros(N_chipova,1);
    Simboli = zeros(N,N_chipova);
    Kodirani_Simboli = zeros(N_chipova,1);
    
    % Generisanje Gold-ove sekvence : ista je za sve tr_kore
    
    % Produzavanje OVSF na trajanje signala
    for i_grana = 1:2
        ceo_niz_ovsf = zeros(1,SF);
        ceo_niz_ovsf(1,1:SF) = ovsf_seq(2*tr_kor+i_grana-2,1:SF);
        ceo_niz_ovsf = repmat(ceo_niz_ovsf,sim_info_duz);
        OVSF_Sekvenca(:,i_grana) = ceo_niz_ovsf(1,:);
    endfor 
        
    % Binarna sekvenca
    x = randsrc(N_bita,1,[1 0; 1/2 1/2]);

    % Mapiranje bita u kompleksne simbole za QPSK 
    sim_info = zeros(sim_info_duz,1);
    for br = 1:2:(N_bita-1)
      temp = x(br:br+1);
      if (temp(1) == 1) && (temp(2) == 1)
        faza = pi/4;
      elseif (temp(1) == 1) && (temp(2) == 0)
        faza = 3*pi/4;
      elseif (temp(1) == 0) && (temp(2) == 0)
        faza = 5*pi/4;
      elseif (temp(1) == 0) && (temp(2) == 1)
        faza = 7*pi/4;
      endif
      sim_info((br+1)/2,1) = e^(i*faza);
    endfor

    % Deljenje simbola na duzinu chipova, zbog mnozenja sa PSS
    chip_info = prosiri(sim_info,SF);
    Simboli_ovsf = zeros(N_chipova,1); 
    Simboli_ovsf(1:N_chipova,1) = real(chip_info(1:N_chipova,1)).*OVSF_Sekvenca(1:N_chipova,1) + 1i*imag(chip_info(1:N_chipova,1)).*OVSF_Sekvenca(1:N_chipova,2);
    Kodirani_Simboli(1:N_chipova,1) = Simboli_ovsf(1:N_chipova,1).*Gold_Sekvenca(1:N_chipova,1);
      
    if tau(tr_kor,1) == 0     % Idealna sinhronizacija PSS
        Signal = sqrt(RaspSnage(tr_kor,1))*Kodirani_Simboli*GreskaSinhFaze(tr_kor,1);
        N_odb = N_chipova;    % deljenje chipova na odbirke
    else                      % Neidealna sinhronizacija PSS
        N_odb = N_chipova / tau(tr_kor,1) - 1;
        Signal = zeros(N_odb,1);
        Kodirani_Simboli = prosiri(Kodirani_Simboli,1/tau(tr_kor,1));

        Signal(1:N_odb,1) = sqrt(RaspSnage(tr_kor,1))*Kodirani_Simboli(2:length(Kodirani_Simboli),1)*GreskaSinhFaze(tr_kor,1);
    endif

    % Dodavanje suma i interferencije
    Es_pN_dB = Eb_pN_dB * log2(M);
    Vsim = Vb(tr_kor,1) / log2(M);
    Vchip = Vsim * SF;
    B = Vchip/2*(1+ro);
    SNRdB = Es_pN_dB + 10*log10(Vsim) - 10*log10(B);
    %SNRdB = Eb_pN_dB;
    SNRdB = 20;
    
    Prx = mean(abs(sqrt(RaspSnage(tr_kor,1))*Kodirani_Simboli));
    PrxdB = 10*log10(Prx);
    NoiseLevel_dBw = PrxdB - SNRdB + 10*log10(SF);       % Nivo kompleksnog AWGN na ulazu u prijemnik
    ComplexNoise = zeros(1:N_odb,1);
    ComplexNoise(1:N_odb,1) = wgn(1, N_odb, NoiseLevel_dBw, 'complex');  % Odbirci kompleksnog AWGN na prijemu
    Signal_plus_N = Signal + ComplexNoise;
    %Signal_plus_N = awgn(Signal, SNRdB, 'measured');
    Signal_plus_I = awgn(Signal, SIRvektor(tr_kor,1), 'measured');
    %Signal_plus_N_plus_I = Signal_plus_N + Signal_plus_I - Signal;
    Signal_plus_N_plus_I = ComplexNoise + Signal_plus_I;
    
    
    Gold_Sekvenca_pros = Gold_Sekvenca;
    OVSF_Sekvenca_pros = OVSF_Sekvenca;
    
    if tau(tr_kor,1) != 0    
      Gold_Sekvenca_pros = prosiri(Gold_Sekvenca,1/tau(tr_kor,1));
      OVSF_Sekvenca_pros = zeros(N_chipova/tau(tr_kor,1),2);
      OVSF_Sekvenca_pros(:,1) = prosiri(OVSF_Sekvenca(:,1),1/tau(tr_kor,1));
      OVSF_Sekvenca_pros(:,2) = prosiri(OVSF_Sekvenca(:,2),1/tau(tr_kor,1));
    endif 
    
    % Mnozenje sa Gold na prijemu
    Kompres1(1:N_odb,1) = Signal_plus_N_plus_I(1:N_odb,1).*Gold_Sekvenca_pros(1:N_odb,1);
    Kompres2(1:N_odb,1) = real(Kompres1(1:N_odb,1)).*OVSF_Sekvenca_pros(1:N_odb,1) + 1i*imag(Kompres1(1:N_odb,1)).*OVSF_Sekvenca_pros(1:N_odb,2);
    
    %Mnozenje sa OVSF na prijemu
    Kompres2 = downsample(Kompres2,SF);
    if tau(tr_kor,1) != 0
      Kompres2 = downsample(Kompres2,1/tau(tr_kor,1));
    endif 

    % Odlucivanje : kombinacija znaka realnog i imaginarnog dela odgovara kvadrantu
    sim_info_rx = (sign(real(Kompres2)) + i*sign(imag(Kompres2)))/sqrt(2);

    % Demapiranje simbola u bite
    y = zeros(N_bita,1);
    dva_bita = zeros(1,2);
    
    for k = 1:sim_info_duz
      re = real(sim_info_rx(k));
      im = imag(sim_info_rx(k));
      if (re > 0) && (im > 0)
        dva_bita = [1 1];
      elseif (re < 0) && (im > 0)
        dva_bita = [1 0];
      elseif (re < 0) && (im < 0)
        dva_bita = [0 0];
      elseif (re > 0) && (im < 0)
        dva_bita = [0 1];  
      endif
        y(2*(k-1)+(1:2),1) = dva_bita;
    endfor

    Peb_izmereno(tr_kor,1) = mean(x~=y);
    
endfor %tr_kor

Snage_na_prijemu = RaspSnage'                             

Bitski_Protoci = Vb'

disp('Verovatnoca greske po korisnicima:')
Peb_izmereno'

disp('Verovatnoca greske za konv. QPSK:')
Peb_konv = Peb_zadato

disp("\n\n\n")