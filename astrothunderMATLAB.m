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
deltat = 0.001; %s
time = 0;
i = 1;



%Initial Calculations
throatArea = pi*((throatDiameter/2)^2);
grainWidth = grainOuterDiameter./2 - grainInnerDiameter./2; 
capArea = 0.25*pi*(grainOuterDiameter^2) - 0.25*pi*(grainInnerDiameter^2);
propVolume = capArea*grainLength*numGrains;

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

function [exitArea] = exitAreaFunc(exitMachNum)
    [~, ~, ~, ~, areaRatio] = flowisentropic(ratioSpecHeats, exitMachNum, 'mach')
    exitArea = areaRatio.*throatArea;
end

function [exitPressure] = exitPressureFunc(exitMachNum)
    [~, ~, pressureRatioExit] = flowisentropic(ratioSpecHeats, exitMachNum, 'mach');
    exitPressure = MEOP.*pressureRatioExit;
end

%Stored Data
initialSize = 1000;
chamberPressureVec = zeros(1, initialSize);
massFlowVec = zeros(1, initialSize);
timeVec = zeros(1, initialSize);

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
        timeVec = [timeVec, zeros(1, initialSize)];
    end

    
end

%Removes zeros at the end of the vectors
chamberPressureVecMask = chamberPressureVec ~= 0;
chamberPressureVec = chamberPressureVec(chamberPressureVecMask);
massFlowVecMask = massFlowVec ~= 0;
massFlowVec = massFlowVec(massFlowVecMask);
timeVecMask = timeVec ~= 0;
timeVec = timeVec(timeVecMask);


%Final Calculations
avgChamberPressure = sum(chamberPressureVec)./length(chamberPressureVec);
MEOP = max(chamberPressureVec);
exitMachNum = exitMachNumFunc(MEOP);
exhaustTemp = exhaustTempFunc(exitMachNum);
localSpeedOfSound = localSpeedOfSoundFunc(exhaustTemp)
exhaustVel = exhaustVelFunc(exitMachNum, localSpeedOfSound);
specificImpulse = specificImpulseFunc(exhaustVel);
propWeight = propDensity.*propVolume*32.174;
exitArea = exitAreaFunc(exitMachNum);
exitPressure = exitPressureFunc(exitMachNum);
thrustVec = (exhaustVel.*massFlowVec) + (exitArea.*(exitPressure - ambientPressure));
averageThrust = sum(thrustVec*deltat)./time;
maxThrust = max(thrustVec);



%Graphing
figure;
plot(timeVec, chamberPressureVec);
title('Pressure Chamber vs. Time');
figure;
plot(timeVec, massFlowVec);
title('Mass Flow Rate vs. Time');
figure;
plot(timeVec, thrustVec);
title('Thrust vs. Time');

%Outputs
totalImpulse = 'test';
specificImpulse = [num2str(specificImpulse) ' s'];
avgChamPressure = [num2str(avgChamberPressure) ' psi'];
maxThrust = [num2str(maxThrust) ' lbs'];
avgThrust = [num2str(averageThrust) ' lbs'];
exhaustVel = [num2str(exhaustVel) ' ft / s'];
burnTime = [num2str(time) ' s'];
maxMassFlux = 'test';
rocketImpulseClass = 'test';
designPressureRatio = 'test';
portThroatAreaRatio = 'test';
propWeight = [num2str(propWeight) ' lbs'];
volumetricLoadingFraction = 'test';
optimumAreaPerfExpansionMEOP = 'test';
ratioInnerGrainAreaToThroatArea = 'test';
exitMachNum = num2str(exitMachNum);
exitTemp = [num2str(exhaustTemp) ' R'];

end       


%Test Case
%[totalImpulse, specificImpulse, MEOP, avgChamPressure, maxThrust, avgThrust, exhaustVel, burnTime, maxMassFlux, rocketImpulseClass, designPressureRatio, portThroatAreaRatio, propWeight, volumetricLoadingFraction, optimumAreaPerfExpansionMEOP, ratioInnerGrainAreaToThroatArea, exitMachNum, exitTemp] = astrothunderMATLAB(1.8, 3.239, 7.5, 5, 1.25)  