function[totalImpulse, specificImpulse, MEOP, avgChamPressure, maxThrust, avgThrust, exhaustVel, burnTime, maxMassFlux, rocketImpulseClass, designPressureRatio, portThroatAreaRatio, propWeight, volumetricLoadingFraction, optimumAreaPerfExpansionMEOP, ratioInnerGrainAreaToThroatArea, exitMachNum] = astrothunderMATLAB(grainOuterDiameter, grainInnerDiameter, grainLength, numGrains, throatDiameter)

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
time = 0;
i = 1;
pie = 3.14159265359;


throatArea = pie*((throatDiameter/2)^2);
%recursion: program runs until propellant weight = 0 (initPropWeight - sum(burnTime.*deltat)) or (propDensity*volume)


%Functions
function [exposedBurnArea] = exposedBurnAreaFunc(grainOuterDiameter, grainInnerDiameter, grainLength)
    capArea = 0.25*pie*(grainOuterDiameter^2) - 0.25*pie*(grainInnerDiameter^2);
    innerGrainArea = pie*grainInnerDiameter*grainLength;
    exposedBurnArea = (2*capArea + innerGrainArea)*numGrains;
end

function [ratioBurnAreaToThroatArea] = burnToThroatFunc(exposedBurnArea, throatArea)
    ratioBurnAreaToThroatArea = exposedBurnArea./throatArea;
    end

function [chamberPressure] = chamberPressureFunc(ratioBurnAreaToThroatArea)
    chamberPressure = (burnRateCoef*charVel*propDensity*ratioBurnAreaToThroatArea)^(1/(1-burnRateExp));
end

function [burnSurfaceRegressionRate] = burnSurfaceRegressionRateFunc(chamberPressure)
    burnSurfaceRegressionRate = burnRateCoef*(chamberPressure^burnRateExp);
end

function [massFlowRate] = massFlowRateFunc(burnSurfaceRegressionRate, exposedBurnArea)
    massFlowRate = propDensity*exposedBurnArea*burnSurfaceRegressionRate;
end

%Stored Data
initialSize = 1000;
chamberPressureVec = zeros(1, initialSize);

%Iteration Loop
while grainLength && grainWidth > 0

    %Grain geometry calculations
    exposedBurnArea = exposedBurnAreaFunc(grainOuterDiameter, grainInnerDiameter, grainLength);
    ratioInnerGrainAreaToThroatArea = burnToThroatFunc(exposedBurnArea, throatArea);
    
    %Chamber Pressure
    chamberPressure = chamberPressureFunc(ratioBurnAreaToThroatArea);
    chamberPressureVec(i) = chamberPressure;

    %Burn Surface Regression Rate
    burnSurfaceRegressionRate = burnSurfaceRegressionRateFunc(chamberPressure);

    %Mass Flow Rate
    massFlowRate = massFlowRateFunc(burnSurfaceRegressionRate, exposedBurnArea);

    %New Grain Geometry Calculations
    grainInnerDiameter = grainInnerDiameter + burnSurfaceRegressionRate*deltat;
    grainLength = grainLength - burnSurfaceRegressionRate*deltat;

    %Iteration
    time = time + deltat;
    i = i + 1;

    %Checks if the number of iterations is greater than the number of preallocated values in the vectors, and if needed doubles the lengths of those vectors
    if i > length(chamberPressureVec)
        chamberPressureVec = [chamberPressureVec, zeroes(1, initialSize)];
    end

    MEOP  = max(chamberPressureVec);
end
end