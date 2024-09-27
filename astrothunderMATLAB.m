function[totalImpulse, specificImpulse, MEOP, avgChamPressure, maxThrust, avgThrust, exhaustVel, burnTime, maxMassFlux, rocketImpulseClass, designPressureRatio, portThroatAreaRatio, propWeight, volumetricLoadingFraction, optimumAreaPerfExpansionMEOP, ratioInnerGrainAreaToThroatArea, exitMachNum] = astrothunderMATLAB(grainOuterDiameter, grainInnerDiameter,grainLength, numGrains, throatDiameter)

%Constants
charVel = 5020; %ft/s
burnRateCoef = 0.0375;
burnRateExp = 0.3;
propDensity = 0.0018; %slugs/in^3
stagTemp = 4700; %R
propSpecGasConstant = 200; %psi*ft^3/slugs*R
ratioSpecHeats = 1.2;
g = 32.2; %ft/s^2
ambientPressure = 14.7; %psi

deltat = 0.001; %s




end
