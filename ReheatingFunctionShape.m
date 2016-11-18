clear all, close all, clf, clc


%******************************************************************************
% PARAMETER DEFINITIONS *******************************************************
%*****************************************************************************/

YS = 5.4e8;
E = 1.96e11;
eps_yield = YS/E;
cte = 1.358e-5;

solidusT = 1873;
onsetT = 0.4*solidusT;
dT_ramp = 0.25*(solidusT - onsetT);

eps_saturation = -1.75e-1*eps_yield;

minT = 0.3*solidusT;
maxT = 0.8*solidusT;

critDT = eps_yield / cte;

dTs = [0:dT_ramp/15:(solidusT-onsetT)];

%******************************************************************************
% SUPERPOSED FUNCTIONS ********************************************************
%*****************************************************************************/

slope_a = eps_yield / critDT;
slope_c = eps_saturation / dT_ramp;

funcA = zeros(length(dTs), 1);
funcB = zeros(length(dTs), 1);
funcC = zeros(length(dTs), 1);
for i=1:length(dTs)
  minDT = min(dT_ramp, critDT);
  maxDT = max(dT_ramp, critDT);
  if (dTs(i) <= critDT)
    funcA(i) = slope_a*dTs(i) - eps_yield;
  else
    funcA(i) = 0;
  end

  if (dTs(i) <= minDT)
    funcB(i) = dTs(i) * (slope_c - slope_a) + eps_yield;
  elseif (dTs(i) > minDT && dTs(i) <= maxDT)
    if (dT_ramp <= critDT)
      funcB(i) = eps_saturation - slope_a*dTs(i) + eps_yield;
    else
      funcB(i) = slope_c*dTs(i);
    end
  else
    funcB(i) = eps_saturation;
  end
end

funcC = funcA + funcB;
figure(1)
plot(dTs, funcA, dTs, funcB, dTs, funcC, '--')
legend('a', 'b', 'c')
