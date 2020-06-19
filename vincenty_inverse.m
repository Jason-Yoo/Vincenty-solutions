function [distance,Alpha1,Alpha2]=vincenty_inverse(point1, point2)
% point[latitude  longitude ]
% distance between two points 
% WGS 84
a = 6378137;  % meters
f = 1 / 298.257223563;
b = 6356752.314245;  % meters; b = (1 - f)a

%MILES_PER_KILOMETER = 0.621371; 
MILES_PER_KILOMETER = 1;
MAX_ITERATIONS = 300;
CONVERGENCE_THRESHOLD = 1e-12;  % .000,000,000,001
convergenceFlag = 0;

% short-circuit coincident points
if (point1(1) == point2(1) && point1(2) == point2(2))
    distance = 0;
    Alpha1 = 0;
    Alpha2 = 0;
    return;
end
U1 = atan((1 - f) * tan(deg2rad(point1(1))));
U2 = atan((1 - f) * tan(deg2rad(point2(1))));
L = deg2rad(point2(2) - point1(2));
Lambda = L;

sinU1 = sin(U1);
cosU1 = cos(U1);
sinU2 = sin(U2);
cosU2 = cos(U2);

    while (--MAX_ITERATIONS>0)
        sinLambda = sin(Lambda);
        cosLambda = cos(Lambda);
        sinSigma = sqrt((cosU2 * sinLambda) ^ 2 + (cosU1 * sinU2 - sinU1 * cosU2 * cosLambda) ^ 2);
        if(sinSigma == 0)
             break;  % coincident points
        end
        cosSigma = sinU1 * sinU2 + cosU1 * cosU2 * cosLambda;
        sigma = atan2(sinSigma, cosSigma);
        sinAlpha = cosU1 * cosU2 * sinLambda / sinSigma;
        cosSqAlpha = 1 - sinAlpha ^ 2;

        cos2SigmaM = cosSigma - 2 * sinU1 * sinU2 / cosSqAlpha;
        if (isnan(cos2SigmaM))
            cos2SigmaM = 0;
        end
        C = f / 16 * cosSqAlpha * (4 + f * (4 - 3 * cosSqAlpha));
        LambdaPrev = Lambda;
        Lambda = L + (1 - C) * f * sinAlpha * (sigma + C * sinSigma *(cos2SigmaM + C * cosSigma *(-1 + 2 * cos2SigmaM ^ 2)));
        if (abs(Lambda - LambdaPrev) < CONVERGENCE_THRESHOLD)
            convergenceFlag = 1;
            break  % successful convergence
        end

    end
    if(convergenceFlag == 1)
        uSq = cosSqAlpha * (a ^ 2 - b ^ 2) / (b ^ 2);
        A = 1 + uSq / 16384 * (4096 + uSq * (-768 + uSq * (320 - 175 * uSq)));
        B = uSq / 1024 * (256 + uSq * (-128 + uSq * (74 - 47 * uSq)));
        deltaSigma = B * sinSigma * (cos2SigmaM + B / 4 * (cosSigma *(-1 + 2 * cos2SigmaM ^ 2) - B / 6 * cos2SigmaM *(-3 + 4 * sinSigma ^ 2) * (-3 + 4 * cos2SigmaM ^ 2)));
        s = b * A * (sigma - deltaSigma);
        Alpha1 = atan(cosU2 * sinLambda / cosU1 * sinU2 - sinU1 * cosU2 * cosLambda);
        Alpha2 = atan(cosU1 * sinLambda / (-sinU1) * cosU2 + cosU1 * sinU2 * cosLambda);
        distance = s * MILES_PER_KILOMETER;  % kilometers to miles
    else
        distance = 0;
        Alpha1 = 0;
        Alpha2 = 0;
    end

end