function [Weights] = Func_Length2Weight(Lengths,y,z)

% Weights of length i
%
% units of weight will depend on lenght: if length is cm, weights are kg,
% if length is mm weigth is mm
%
% Inputs: Lenghts = vector of lengths 
%         y = length-weight parameter
%         z = length-weight parameter
% Outputs: Weights = vector of weights corresponding to the input lengths

Weights = NaN(length(Lengths),1);

     for i=1:length(Lengths)
        Weights(i) = y*(Lengths(i)).^z;
     end

end
