function [NLL, propCW] = giveNLL(params, oriVals, nCCW, nCW, model_type)

for iPF = 1:2   % There are two psychometric functions to be fit
    
    % Decode function arguments
    if model_type(2) == '1'
        lapseRate = params(iPF);
        noiseInt  = params(3);
        decCrit   = params(3+iPF);
    elseif model_type(2) == '2'
        lapseRate = params(iPF);
        noiseInt  = params(2+iPF);
        decCrit   = params(4+iPF);
    end
    propCW(iPF,:) = lapseRate + (1 - 2*lapseRate)*normcdf((oriVals - decCrit)/(2*noiseInt));
%     NLL(iPF)      = -sum(log(binopdf(nCW(iPF,:), nCCW(iPF,:) + nCW(iPF,:), propCW(iPF,:))));
    NLL(iPF) = -sum(log(max(1e-300, binopdf(nCW(iPF,:), nCCW(iPF,:) + nCW(iPF,:), propCW(iPF,:)))));
end
<<<<<<< HEAD
% ==> Robbe's comment here!
NLL = sum(NLL);
end
=======
% ==> sum NLL
NLL = sum(NLL);
end
>>>>>>> master
