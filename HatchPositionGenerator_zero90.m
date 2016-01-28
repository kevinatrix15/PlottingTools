%% Generates hatch files for layers alternating between 0 and 90 degree scan angles
clear all, clc, clf

hatch_spacing = 80e-6;             % [m]
stripe_width0 = 4e-3;             % Stripe width for 0 degree angles

part_x = 40e-3;
part_y = 2e-3;

num_stripes0 = ceil(part_x/stripe_width0);
% lines_per_stripe0 = round(part_y/hatch_spacing)+1;
lines_per_stripe0 = 24;

start_x = zeros(num_stripes0*lines_per_stripe0, 1);
start_y = zeros(num_stripes0*lines_per_stripe0, 1);
end_x = zeros(num_stripes0*lines_per_stripe0, 1);
end_y = zeros(num_stripes0*lines_per_stripe0, 1);

i = 1
for s = 1:num_stripes0
    for l = 1:lines_per_stripe0
        if mod(s,2)==1              % Odd stripes
            if l==1
                start_x(i) = (s-1)*stripe_width0;
                start_y(i) = 0;
                
                end_x(i) = start_x(i) + stripe_width0;
                end_y(i) = start_y(i);
            elseif mod(l,2) == 1        % Odd scan lines
                start_x(i) = end_x(i-1);
                start_y(i) = end_y(i-1) + hatch_spacing;
                
                end_x(i) = start_x(i) + stripe_width0;
                end_y(i) = start_y(i);
            else                        % Even scan lines
                start_x(i) = end_x(i-1);
                start_y(i) = end_y(i-1) + hatch_spacing;
                
                end_x(i) = start_x(i) - stripe_width0;
                end_y(i) = start_y(i);
            end
            
        else                            % Even stripes
            if l==1
                start_x(i) = (s-1)*stripe_width0 + stripe_width0;
                start_y(i) = 0;
                
                end_x(i) = start_x(i) + stripe_width0;
                end_y(i) = start_y(i);
            elseif mod(l,2) == 0        % Even scan lines
                start_x(i) = end_x(i-1);
                start_y(i) = end_y(i-1) + hatch_spacing;
                
                end_x(i) = start_x(i) + stripe_width0;
                end_y(i) = start_y(i);
            else                        % Odd scan lines
                start_x(i) = end_x(i-1);
                start_y(i) = end_y(i-1) + hatch_spacing;
                
                end_x(i) = start_x(i) - stripe_width0;
                end_y(i) = start_y(i);
            end
        end
        i = i+1;
    end
end

figure(1)
plot(start_x, start_y, 'bo', end_x, end_y, 'rx')
