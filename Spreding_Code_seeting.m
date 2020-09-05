%-------------STBC-MIMO-GFDM---------------%
clear all;
%****************************************************************
fprintf( 'STBC-MIMO-GFDM 2x2 \n\n') ;
%****************************************************************
global Num_Sub_Per_Group
global Num_Bit_Per_Group
global Num_Space_Bit_Per_Group
global Walsh_Metric
global index
%---------IM-SS setting------------------%
Num_Sub_Per_Group = 4;
index = 2;
Num_Space_Bit_Per_Group = floor(log2(Num_Sub_Per_Group));
Num_Bit_Per_Group = Num_Space_Bit_Per_Group + index;
%---------Spreading Code-----------------%
%Walsh_Metric = hadamard(Num_Sub_Per_Group);
Chirp = exp(1i*pi/Num_Sub_Per_Group*[0:Num_Sub_Per_Group-1].^2);
Walsh_Metric = zeros(Num_Sub_Per_Group);
for ii = 0:Num_Sub_Per_Group-1
    Walsh_Metric(:,ii+1) = circshift(Chirp.',ii).';
end