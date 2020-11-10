
%..........................................................................
% This Program generates 10 Gold code sequences, each 31 bits long. These
% codes are outputted as coloumns of matrix Co_Mat.
% Programed by Imran Ali (Mehran University, Jamshoro, Sindh, Pakistan)
%..........................................................................


function gold_seq = gold(G)

%.................Generation of first perferred PN sequence................
sd1 =randsrc(1,5,[0 1]);      % First user's seed.
PN1=[];                         % Spreading code vector of user-1
for j=1:G        
    PN1=[PN1 sd1(1)];
    if sd1(1)==sd1(4)
        temp1=0;
    else temp1=1;
    end
    sd1(1)=sd1(2);
    sd1(2)=sd1(3);
    sd1(3)=sd1(4);
    sd1(4)=sd1(5);
    sd1(5)=temp1;
end
%..........................................................................
%.................Generation of Second perferred PN sequence...............
sd2 =randsrc(1,5,[0 1]);      
PN2=[];                        
for j=1:G        
    PN2=[PN2 sd2(1)];
    if sd2(1)==sd2(2)
        temp1=0;
    else temp1=1;
    end
    if sd2(4)==temp1
        temp2=0;
    else temp2=1;
    end
    if sd2(5)==temp2
        temp3=0;
    else temp3=1;
    end
    sd2(1)=sd2(2);
    sd2(2)=sd2(3);
    sd2(3)=sd2(4);
    sd2(4)=sd2(5);
    sd2(5)=temp3;
end
%..........................................................................
%.........................Generation of Gold Codes.........................
Co_Mat=[];

    code=[];
    PN2(31)=PN2(1);
    for k=1:G-1
        PN2(k)=PN2(k+1);
    end
    for j=1:G
        code=[code xor(PN1(j),PN2(j))];
    end
    Co_Mat=[Co_Mat code'];    

for row=1:G
        if Co_Mat(row,1)==0
            Co_Mat(row,1)=-1;
        end
end            

%..........................................................................
%....................Chechking corelation performance......................
A=[]; AA=0;B=[];
  
        for k=1:G
            AA=AA+Co_Mat(k,1)*Co_Mat(k,1);
        end
        A=[A AA];
        AA=0;
 
    B=[B A'];
    A=[];

gold_seq = Co_Mat';

%..........................................................................
% FINAL COMMENTS: The matrix B will be 10-by-10 whose diagonal values will
% be 31 always, showing correlation peak of same codes. Elsewhere in the
% matix B, you will see small values, showing good correlation performace
% among diffrent codes. HAVE A NICE DAY.
%..........................................................................
