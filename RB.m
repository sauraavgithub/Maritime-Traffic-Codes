function [RB] = RB(x0,y0,x,y,c0)

% Calculating Relative bearing.
SH = c0; %Ship heading/Course

if(x<=x0)
    if(y<=y0)
      T = pi + atan(abs((x-x0)/(y-y0)));
    else
      T = 2*pi - atan(abs((x-x0)/(y-y0)));
    end
    
else 
    if(y<=y0)
      T = pi/2 + atan(abs((y-y0)/(x-x0)));
    else 
      T = atan(abs((x-x0)/(y-y0)));
    end
end

 RB = abs(T -SH) ;       

end










