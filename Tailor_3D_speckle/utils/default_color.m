function [a] = default_color(char)
if char == 'b'
    a = [0 0.4470 0.7410];
elseif char == 'r'
    a = [0.85 0.3250 0.0980];
elseif char == 'y'
    a = [0.9290 0.6940, 0.1250];
elseif char == 'p'
    a = [0.4940 0.1840 0.5560];
elseif char == 'g'
    a = [0.4660 0.6740 0.1880];
elseif char == 'c'
    a = [0.3010 0.7450 0.9330];
elseif char == 'k'
    a = [0 0 0];
elseif char == 'w' 
    a = [0.6350 0.0780 0.1840]; %more like wine color
elseif char == 'm'
    a = [1 0.07 0.65]; %more like wine color
else
    disp('Please type among b,r,y,p,g,c,w,m')
    a = [rand(1) rand(1) rand(1)];
end