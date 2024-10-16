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
exitArea = 7.1967; %in



%Initial Calculations
throatArea = pi*((throatDiameter/2)^2);
grainWidth = grainOuterDiameter./2 - grainInnerDiameter./2; 
capArea = 0.25*pi*(grainOuterDiameter^2) - 0.25*pi*(grainInnerDiameter^2);
propVolume = capArea*grainLength*numGrains;
areaRatio = exitArea./throatArea;

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

function [exitMachNum, exhaustTemp, exitPressure] = isentropicFunc(areaRatio)
    [exitMachNum, tempRatio, pressureRatioExit] = flowisentropic(ratioSpecHeats, areaRatio, 'sup');
    exhaustTemp = stagTemp.*tempRatio;
    exitPressure = chamberPressure.*pressureRatioExit;
end

%function [pressureRatioExit, exitPressure] = exitPressureFunc(exitMachNum)
%    [~, ~, pressureRatioExit] = flowisentropic(ratioSpecHeats, exitMachNum, 'mach');
%    exitPressure = chamberPressure.*pressureRatioExit;
%end

function [localSpeedofSound] = localSpeedOfSoundFunc(exhaustTemp)
    localSpeedofSound = sqrt(ratioSpecHeats*propSpecGasConstant*exhaustTemp);
end

function [exhaustVel] = exhaustVelFunc(exitMachNum, localSpeedOfSound)
    exhaustVel = exitMachNum*localSpeedOfSound;
end

function [specificImpulse] = specificImpulseFunc(thrustVec, massFlowVec)
    specificImpulse = thrustVec/(massFlowVec*g);
end

%function [exitArea] = exitAreaFunc(exitMachNum)
 %   [~, ~, ~, ~, areaRatio] = flowisentropic(ratioSpecHeats, exitMachNum, 'mach');
 %   exitArea = areaRatio.*throatArea;
%end

function [thrust] = thrustFunc(massFlowRate, exhaustVel, exitArea, exitPressure)
    thrust = (exhaustVel.*massFlowRate) + (exitArea.*(exitPressure - ambientPressure));
end

%Stored Data
initialSize = 1000;
chamberPressureVec = zeros(1, initialSize);
massFlowVec = zeros(1, initialSize);
timeVec = zeros(1, initialSize);
thrustVec = zeros(1, initialSize);
exhaustVelVec = zeros(1, initialSize);
exitPressureVec = zeros(1, initialSize);

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

    %Exhaust Velocity
    [exitMachNum, exhaustTemp, exitPressure] = isentropicFunc(areaRatio);
    exitPressureVec(i) = exitPressure;
    
    localSpeedOfSound = localSpeedOfSoundFunc(exhaustTemp);
    exhaustVel = exhaustVelFunc(exitMachNum, localSpeedOfSound);    
    exhaustVelVec(i) = exhaustVel;

    %exitArea = exitAreaFunc(exitMachNum);
    
    %Thrust
    thrust = thrustFunc(massFlowRate, exhaustVel, exitArea, exitPressure);
    thrustVec(i) = thrust;

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
        thrustVec = [thrustVec, zeros(1, initialSize)];
        exhaustVelVec = [exhaustVelVec, zeros(1, initialSize)];
        exitPressureVec = [exitPressureVec, zeros(1, initialSize)];
    end

    
end

%Removes zeros at the end of the vectors
chamberPressureVecMask = chamberPressureVec ~= 0;
chamberPressureVec = chamberPressureVec(chamberPressureVecMask);
massFlowVecMask = massFlowVec ~= 0;
massFlowVec = massFlowVec(massFlowVecMask);
timeVecMask = timeVec ~= 0;
timeVec = timeVec(timeVecMask);
thrustVecMask = thrustVec ~= 0;
thrustVec = thrustVec(thrustVecMask);
exhaustVelVecMask = exhaustVelVec ~= 0;
exhaustVelVec = exhaustVelVec(exhaustVelVecMask);

%Final Calculations
avgChamberPressure = sum(chamberPressureVec)./length(chamberPressureVec);
avgExhaustVel = sum(exhaustVelVec)./length(exhaustVelVec);
MEOP = max(chamberPressureVec);
[exitMachNum, exhaustTemp, exitPressure] = isentropicFunc(areaRatio);
localSpeedOfSound = localSpeedOfSoundFunc(exhaustTemp);
exhaustVel = exhaustVelFunc(exitMachNum, localSpeedOfSound);
specificImpulse = specificImpulseFunc(thrustVec, massFlowVec);
propWeight = propDensity.*propVolume*32.174;
averageThrust = sum(thrustVec)./length(thrustVec);
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
%MAKE AN EXIT PRESSURE VECTOR TO REPLACE EXIT VELOCITY VECTOR