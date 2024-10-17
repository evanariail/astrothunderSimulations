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
chamberVolume = (0.25*pi*(grainOuterDiameter^2))*(grainLength*numGrains);
portArea = pi*(grainInnerDiameter./2)^2;

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

function [exitMachNum, exhaustTemp, exitPressure, pressureRatioExit] = isentropicFunc(areaRatio)
    [exitMachNum, tempRatio, pressureRatioExit] = flowisentropic(ratioSpecHeats, areaRatio, 'sup');
    exhaustTemp = stagTemp.*tempRatio;
    exitPressure = chamberPressure.*pressureRatioExit;
end

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
    %[~, ~, ~, ~, areaRatio] = flowisentropic(ratioSpecHeats, exitMachNum, 'mach');
    %exitArea = areaRatio.*throatArea;
%end

function [thrust] = thrustFunc(massFlowRate, exhaustVel, exitArea, exitPressure)
    thrust = (exhaustVel.*massFlowRate) + (exitArea.*(exitPressure - ambientPressure));
end

function [totalImpulse] = totalImpulseFunc(thrustVec)
    totalImpulse = sum(thrustVec).*deltat;
end

%Stored Data
initialSize = 1000;
chamberPressureVec = zeros(1, initialSize);
massFlowVec = zeros(1, initialSize);
timeVec = zeros(1, initialSize);
thrustVec = zeros(1, initialSize);
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
    [exitMachNum, exhaustTemp, exitPressure, ~] = isentropicFunc(areaRatio);
    exitPressureVec(i) = exitPressure;
    
    localSpeedOfSound = localSpeedOfSoundFunc(exhaustTemp);
    exhaustVel = exhaustVelFunc(exitMachNum, localSpeedOfSound);    
    
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

%Final Calculations
avgChamberPressure = sum(chamberPressureVec)./length(chamberPressureVec);
MEOP = max(chamberPressureVec);
[exitMachNum, exhaustTemp, ~, pressureRatioExit] = isentropicFunc(areaRatio); 
localSpeedOfSound = localSpeedOfSoundFunc(exhaustTemp);
exhaustVel = exhaustVelFunc(exitMachNum, localSpeedOfSound);
specificImpulse = specificImpulseFunc(thrustVec, massFlowVec);
propWeight = propDensity.*propVolume*32.174;
averageThrust = sum(thrustVec)./length(thrustVec);
maxThrust = max(thrustVec);
totalImpulse = totalImpulseFunc(thrustVec);
volumetricLoadingFraction = propVolume./chamberVolume;
portThroatAreaRatio = portArea/throatArea;
%massFluxVec = (massFlowVec)/(deltat*pi*(grainInnerDiameter./2)^2);

%If re-calculating optimal exit Area ====> Uncomment the exit area function and the line below.
%exitArea = exitAreaFunc(exitMachNum);

%Motor Classification
if totalImpulse >= 576.01 && totalImpulse <= 1151
    rocketImpulseClass = 'L';
elseif totalImpulse >= 1151.01 && totalImpulse <= 2302
    rocketImpulseClass = 'M';
elseif totalImpulse >= 2302.01 && totalImpulse <= 4604
    rocketImpulseClass = 'N';
elseif totalImpulse >= 4604.01 && totalImpulse <= 9208
    rocketImpulseClass = 'O';
elseif totalImpulse >= 9208.01 && totalImpulse <= 18400
    rocketImpulseClass = 'P';
else
    rocketImpulseClass = 'Classification not in range L < motor < P.';
end

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
totalImpulse = sprintf('Total Impulse: %.4f lbs * s', totalImpulse);
specificImpulse = sprintf('Specific Impulse: %.4f R', specificImpulse);
MEOP = sprintf('Maximum Expected Engine Operating Pressure: %.4f psi', MEOP);
avgChamPressure = sprintf('Average Chamber Pressure: %.4f psi', avgChamberPressure);
maxThrust = sprintf('Maximum Thrust: %.4f lbs', maxThrust);
avgThrust = sprintf('Average Thrust: %.4f lbs', averageThrust);
exhaustVel = sprintf('Maximum Exhaust Velocity: %.4f ft / s', exhaustVel);
burnTime = sprintf('Burn Time: %.4f s', time);
maxMassFlux = sprintf('Maximum Mass Flux: %.4f lbs / s / in ^ 2', maxMassFlux);
rocketImpulseClass = sprintf('Motor Class: %s', rocketImpulseClass);
designPressureRatio = sprintf('Design Pressure Ratio: %.4f', pressureRatioExit);
portThroatAreaRatio = sprintf('Port Area to Throat Area Ratio: %.4f', portThroatAreaRatio);
propWeight = sprintf('Propellant Weight: %.4f lbs', propWeight);
volumetricLoadingFraction = sprintf('Volumetric Loading Fraction: %.4f', volumetricLoadingFraction);
optimumAreaPerfExpansionMEOP = sprintf('Optimum Exit Area for Perfect Expansion: %.4f in ^ 2', exitArea);
exitMachNum = sprintf('Exit Mach Number: %.4f', exitMachNum);
exitTemp = sprintf('Exhaust Temperature: %.4f R', exhaustTemp);

end       


%Test Case
%[totalImpulse, specificImpulse, MEOP, avgChamPressure, maxThrust, avgThrust, exhaustVel, burnTime, maxMassFlux, rocketImpulseClass, designPressureRatio, portThroatAreaRatio, propWeight, volumetricLoadingFraction, optimumAreaPerfExpansionMEOP, ratioInnerGrainAreaToThroatArea, exitMachNum, exitTemp] = astrothunderMATLAB(1.8, 3.239, 7.5, 5, 1.25)  
