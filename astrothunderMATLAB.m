function[totalImpulse, specificImpulse, MEOP, avgChamPressure, maxThrust, avgThrust, exhaustVel, burnTime, maxMassFlux, rocketImpulseClass, designPressureRatio, portThroatAreaRatio, propWeight, volumetricLoadingFraction, optimumAreaPerfExpansionMEOP, ratioInnerGrainAreaToThroatArea, exitMachNum, exitTemp] = astrothunderMATLAB(grainInnerDiameter, grainOuterDiameter, grainLength, numGrains, throatDiameter)

%Constants
charVel = 5020; %ft/s
burnRateCoef = 0.0375;
burnRateExp = 0.3;
propDensity = 0.0018; %slugs/in^3
stagTemp = 4700; %R
propSpecGasConstant = 2000; %psi*ft^3/slugs*R
ratioSpecHeats = 1.2;
g = 32.2; %ft/s^2
ambientPressure = 14.7; %psi
deltat = 0.01; %s
time = 0;
i = 1;



%Initial Calculations
throatArea = pi*((throatDiameter/2)^2);
grainWidth = grainOuterDiameter./2 - grainInnerDiameter./2;   

%Functions
function [exposedBurnArea] = exposedBurnAreaFunc(grainOuterDiameter, grainInnerDiameter, grainLength, numGrains)
    capArea = 0.25*pi*(grainOuterDiameter^2) - 0.25*pi*(grainInnerDiameter^2);
    innerGrainArea = pi*grainInnerDiameter*grainLength;
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

function [exitMachNum] = exitMachNumFunc(MEOP)
    pressureRatio = ambientPressure/MEOP;
    exitMachNum = flowisentropic(ratioSpecHeats, pressureRatio, 'pres');
end

function [exhaustTemp] = exhaustTempFunc(exitMachNum)
    [~, tempRatio] = flowisentropic(ratioSpecHeats, exitMachNum, 'mach');
    exhaustTemp = stagTemp.*tempRatio;
end

function [localSpeedofSound] = localSpeedOfSoundFunc(exhaustTemp)
    localSpeedofSound = sqrt(ratioSpecHeats*propSpecGasConstant*exhaustTemp);
end

function [exhaustVel] = exhaustVelFunc(exitMachNum, localSpeedOfSound)
    exhaustVel = exitMachNum*localSpeedOfSound;
end

function [specificImpulse] = specificImpulseFunc(exhaustVel)
    specificImpulse = exhaustVel/g;
end



%Stored Data
initialSize = 1000;
chamberPressureVec = zeros(1, initialSize);
massFlowVec = zeros(1, initialSize);
thrustVec = zeros(1, initialSize);
timeVec = zeros(1, initialSize);
%timeVec(1) = 0;

%Iteration Loop
while grainWidth > 0  && grainLength > 0

    %Grain geometry calculations
    exposedBurnArea = exposedBurnAreaFunc(grainOuterDiameter, grainInnerDiameter, grainLength, numGrains);
    ratioInnerGrainAreaToThroatArea = burnToThroatFunc(exposedBurnArea, throatArea);
    
    %Chamber Pressure
    chamberPressure = chamberPressureFunc(ratioInnerGrainAreaToThroatArea);
    chamberPressureVec(i) = chamberPressure;

    %Burn Surface Regression Rate
    burnSurfaceRegressionRate = burnSurfaceRegressionRateFunc(chamberPressure);

    %Mass Flow Rate
    massFlowRate = massFlowRateFunc(burnSurfaceRegressionRate, exposedBurnArea);
    massFlowVec(i) = massFlowRate;

    %New Grain Geometry Calculations
    grainInnerDiameter = grainInnerDiameter + 2*burnSurfaceRegressionRate*deltat;
    grainLength = grainLength - 2*burnSurfaceRegressionRate*deltat;
    grainWidth = grainOuterDiameter - grainInnerDiameter;

    %Iteration
    time = time + deltat;
    i = i + 1;
    timeVec(i) = time;
    
    %Checks if the number of iterations is greater than the number of preallocated values in the vectors, and if needed doubles the lengths of those vectors
    if i > length(chamberPressureVec)
        chamberPressureVec = [chamberPressureVec, zeros(1, initialSize)];
        massFlowVec = [massFlowVec, zeros(1, initialSize)];
        thrustVec = [thrustVec, zeros(1, initialSize)];
        timeVec = [timeVec, zeros(1, initialSize)];
    end

    
end

%Final Calculations
avgChamberPressure = sum(chamberPressureVec)./length(chamberPressureVec);
MEOP = max(chamberPressureVec);
exitMachNum = exitMachNumFunc(MEOP);
exhaustTemp = exhaustTempFunc(exitMachNum);
localSpeedOfSound = localSpeedOfSoundFunc(exhaustTemp);
exhaustVel = exhaustVelFunc(exitMachNum, localSpeedOfSound);
specificImpulse = specificImpulseFunc(exhaustVel);
plot(timeVec, chamberPressureVec);
%plot(timeVec, massFlowVec)
%plot(timeVec, thrustVec)

%Outputs
totalImpulse = 'test';
%specificImpulse found above
avgChamPressure = avgChamberPressure;
maxThrust = 'test';
avgThrust = 'test';
%exhaustVel found above 
burnTime = 'test';
maxMassFlux = 'test';
rocketImpulseClass = 'test';
designPressureRatio = 'test';
portThroatAreaRatio = 'test';
propWeight = 'test';
volumetricLoadingFraction = 'test';
optimumAreaPerfExpansionMEOP = 'test';
%ratioInnerGrainAreaToThroatArea = ratioInnerGrainAreaToThroatArea;
%exitMachNum found above
exitTemp = exhaustTemp;

end       


%Test Case
%[totalImpulse, specificImpulse, MEOP, avgChamPressure, maxThrust, avgThrust, exhaustVel, burnTime, maxMassFlux, rocketImpulseClass, designPressureRatio, portThroatAreaRatio, propWeight, volumetricLoadingFraction, optimumAreaPerfExpansionMEOP, ratioInnerGrainAreaToThroatArea, exitMachNum, exitTemp] = astrothunderMATLAB(1.8, 3.239, 7.5, 5, 1.25)