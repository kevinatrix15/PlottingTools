function [p] = fittingCoeffSwap(p0)

% This function takes the 1D coefficients and reassigns them for direct use
% in the 2D gaussian function

a1x = p0(1); b1x = p0(2); c1x = p0(3);
a2x = p0(4); b2x = p0(5); c2x = p0(6);
a1y = p0(7); b1y = p0(8); c1y = p0(9);
a2y = p0(10); b2y = p0(11); c2y = p0(12);
d = p0(13);

% Make sure a1's are always > 0 if a2's < 0
if a1x < 0 && a2x > 0
    a1x = a2x;
    a2x = p0(1);
end
if a1y < 0 && a2y > 0
    a1y = a2y;
    a2y = p0(7);
end

% Find how many peak values (a-coeffs) are negative
neg_count = 0;
if a1x < 0
    neg_count = neg_count+1;
end
if a2x < 0
    neg_count = neg_count+1;
end
if a1y < 0
    neg_count = neg_count+1;
end
if a2y < 0
    neg_count = neg_count+1;
end

% Determine which configuration to apply
case_flag = 0;                      % 1: single peak along 1 axis, double
                                    % in other; 2: either all peaks
                                    % concentric, or 2 separate peaks in
                                    % both x' and y'
                                    
tol = 3e-6;                         % TODO: Make more sophistocated
if neg_count == 1 && abs(b2y - b1y) < tol && abs(b2x - b1x) > tol
    % single peak in y', 2 in x'
        if a1x > a2x
            case_flag = '1a';
        else
            case_flag = '1b';
        end
elseif neg_count == 1 && abs(b2x-b1x) < tol && abs(b2y - b1y) > tol
        if a1y > a2y
            case_flag = '2a';
        else
            case_flag = '2b';
        end
%     end
elseif neg_count == 0
    case_flag = '3b';
else
    case_flag = '3a';
end
% if neg_count == 1
%     % single peak in y', 2 in x'
%     if abs(b2y - b1y) < tol && abs(b2x - b1x) > tol
%         if a1x > a2x
%             case_flag = '1a';
%         else
%             case_flag = '1b';
%         end
%     elseif abs(b2x-b1x) < tol && abs(b2y - b1y) > tol
%         if a1y > a2y
%             case_flag = '2a';
%         else
%             case_flag = '2b';
%         end
%     end
% elseif neg_count == 0
%     case_flag = '3b';
% else
%     case_flag = '3a';
% end

% Initialize replacement coefficients
p.a1 = 0.5*a1x;   p.a2 = 0.5*a2y;   p.a3 = 0.5*a1y;   p.a4 = 0.5*a2x;
p.b1 = b1x; p.c1 = c1x; p.b2 = b1y; p.c2 = c1y;
p.b3 = b1x; p.c3 = c1x; p.b4 = b2y; p.c4 = c2y;
p.b5 = b2x; p.c5 = c2x; p.b6 = b1y; p.c6 = c1y;
p.b7 = b2x; p.c7 = c2x; p.b8 = b2y; p.c8 = c2y;
p.d = d;


% Replace values based on configuration
if strcmp(case_flag, '1a')
    p.b4 = b1y;    p.c4 = c1y;
    p.b5 = b1x;    p.c5 = c1x;
    p.a4 = a2x;
elseif strcmp(case_flag, '1b')
    p.b3 = b2x;    p.c3 = c2x;
    p.b6 = b2y;    p.c6 = c2y;
    p.a1 = a1x;
elseif strcmp(case_flag, '2a')
    p.b3 = b2x;    p.c3 = c2x;
    p.b5 = b1x;    p.c5 = c1x;
    p.b7 = b1x;    p.c7 = c1x;
    p.b8 = b1y;    p.c8 = c1y;
    p.a2 = a2y;
elseif strcmp(case_flag, '2b')
    p.b1 = b2x;    p.c1 = c2x;
    p.b2 = b2y;    p.c2 = c2y;
    p.b3 = b2x;    p.c3 = c2x;
    p.b5 = b1x;    p.c5 = c1x;
    p.a3 = a1y;
elseif strcmp(case_flag, '3a')
    p.b3 = b2x;    p.c3 = c2x;
    p.b6 = b2y;    p.c6 = c2y;
    p.b7 = b1x;    p.c7 = c1x;
    p.b8 = b1y;    p.c8 = c1y;
elseif strcmp(case_flag, '3b')
    p.b3 = b2x;    p.c3 = c2x;
    p.b5 = b1x;    p.c5 = c1x;
else
    error('fittingCoeffSwap: case unaccounted for');
end

end

