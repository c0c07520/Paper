function Trans_Sym = Generate_Trans_Sym_GFDM(Random_Bit)
global Num_Sub_Per_Group
global Num_Space_Bit_Per_Group
global Walsh_Metric
global mod_obj
global index
Space_Bit = Random_Bit(1:Num_Space_Bit_Per_Group,:);
Symbol_Bit = Random_Bit(Num_Space_Bit_Per_Group+1:end,:);
Symbol_Bit2 = Symbol_Bit(:).';
Walsh_Mode_Loc = (bi2de(Space_Bit.','left-msb') + 1).';
Information_Symbol = modulation(Symbol_Bit2,index);
%Information_Symbol = modulate(mod_obj,Symbol_Bit);
Trans_Information_Sym = repmat(Information_Symbol,Num_Sub_Per_Group,1);
Trans_Sym = Walsh_Metric(Walsh_Mode_Loc,:).'.*Trans_Information_Sym;
end