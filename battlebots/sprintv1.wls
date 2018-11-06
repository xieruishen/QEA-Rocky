#!/usr/bin/env wolframscript
(* ::Package:: *)

eq1 = theta[s] == Subscript[e, v][s]*Subscript[G, VELOCITYCONTROLLER]*Subscript[G, MOTOR]*Subscript[G, VELOCITYTOANGLE];
eq2 = Subscript[e, v][s] == Subscript[K, p] * (Subscript[e, theta][s]- (0.01/s - v[s]) *Subscript[K, d]/s) + Subscript[e, theta][s]*Subscript[K, i]/s - v[s]
eq3 = v[s] == Subscript[e, v][s]*Subscript[G, VELOCITYCONTROLLER]*Subscript[G, MOTOR]
eq4 = Subscript[e, theta][s] == Subscript[theta, d][s] - theta[s]


sol = Solve[{eq1, eq2,eq3,eq4},{theta[s], Subscript[e, v][s], v[s],Subscript[e, theta][s]}][[1]]
{Subscript[G, thetatotheta][s] -> Simplify[theta[s]/Subscript[theta, d][s]/.sol]}


Subscript[G, TOTALSYSTEM] = Factor[theta[s]/Subscript[theta, d][s]/.sol/.{Subscript[G, VELOCITYCONTROLLER] -> Subscript[J, p]+Subscript[J, i]/s, Subscript[G, MOTOR] -> \[Beta] \[Alpha] /(s+\[Alpha]), Subscript[G, VELOCITYTOANGLE] -> -s/(L s^2- g)}]


constants = {g-> 98/10, L-> 0.04, \[Alpha]-> 10, \[Beta]-> 1/400};
Subscript[theta, d][s]=LaplaceTransform[(UnitStep[t]-UnitStep[t-1/4])/15,t,s];
f[kpval_, kival_, jpval_, jival_, kdval_]:=N[Max[Re[Values[Solve[Denominator[Factor[Subscript[G, TOTALSYSTEM]/.constants/.{Subscript[K, p]->kpval,Subscript[K, i]-> kival,Subscript[J, p]-> jpval,Subscript[J, i]->jival,Subscript[K, d]-> kdval}]]==0,s]]]]]
minimizedsol = NMinimize[f[-6,-36,500,jival,kdval],{kpval,kival,jpval,jival,kdval}]
goodParams = minimizedsol[[2]]/.{kpval-> Subscript[K, p], kival-> Subscript[K, i], jpval -> Subscript[J, p], jival -> Subscript[J, i], kdval -> Subscript[K, d]}


poles = Values[NSolve[Denominator[Subscript[G, TOTALSYSTEM]]==0,s]/.constants/.goodParams]


goodParams = {Subscript[K, p] -> -6, Subscript[K, i] -> -36, Subscript[J, p] -> 500, Subscript[J, i] -> 5000, Subscript[K, d] -> -1}


Subscript[G, fullsystemsubbed] = Subscript[G, TOTALSYSTEM]/.constants/.goodParams


Subscript[G, velocitytoangle] = -s/(L s^2- g)/.constants;
Subscript[theta, d][s]=LaplaceTransform[(UnitStep[t]-UnitStep[t-1/4])/15,t,s]
Subscript[G, fullsystemsubbed]
thetaResponse =InverseLaplaceTransform[Factor[Subscript[G, fullsystemsubbed]],s,t]
vResponse =InverseLaplaceTransform[Factor[Subscript[G, fullsystemsubbed]/ Subscript[G, velocitytoangle]],s,t]
(*Plot[thetaResponse,{t,0,10},PlotRange-> All,AxesLabel->{"time (s)", "angle (radians)"}]
Plot[vResponse,{t,0,10},PlotRange-> All,AxesLabel->{"time (s)", "velocity (m/s)"}]*)


Denominator[Factor[Subscript[G, TOTALSYSTEM]]]/.constants/.{Subscript[J, p] -> 500, Subscript[J, i] -> 5000, Subscript[K, d] -> 1}
poles = Values[NSolve[Denominator[Subscript[G, TOTALSYSTEM]]==0,s]/.constants/.{Subscript[J, p] -> 500, Subscript[J, i] -> 5000, Subscript[K, d] -> 1}];
Manipulate[ListPlot[ReIm[poles/.{Subscript[K, p] -> Kp, Subscript[K, i] -> Ki}],PlotRange-> {{-50,50},{-100,100}},PlotMarkers-> {Graphics[{Blue,Disk[]}],.025},AxesLabel->{"Re","Im"}],{Kp,-1000,1000},{Ki,-1000,1000}]


denominatorCoefs = CoefficientList[Denominator[Factor[Subscript[G, TOTALSYSTEM]/.constants]],s]
denominatorCoefs[[6]]
denominatorCoefs = denominatorCoefs/denominatorCoefs[[6]]


Subscript[p, 1] = -3 - 2\[Pi] I;
Subscript[p, 2] = -3 + 2\[Pi] I;
Subscript[p, 3]= -2 + I;
Subscript[p, 4]= -2 + I;
Subscript[p,5] = -10;
targetDenominator = (s-Subscript[p, 1])(s-Subscript[p, 2])(s-Subscript[p, 3])(s-Subscript[p, 4])(s-Subscript[p,5])
targetDenominatorCoefficients = CoefficientList[targetDenominator,s]
controlParameters = Solve[denominatorCoefs[[1]]==targetDenominatorCoefficients[[1]] && denominatorCoefs[[2]] == targetDenominatorCoefficients[[2]] && denominatorCoefs[[3]] == targetDenominatorCoefficients[[3]] && denominatorCoefs[[4]] == targetDenominatorCoefficients[[4]]&& denominatorCoefs[[5]] == targetDenominatorCoefficients[[5]],{Subscript[J, p],Subscript[J, i],Subscript[K, p], Subscript[K, i]}]
parameters = N[controlParameters]






