

function [ z_chk ] = zone_chk(x,y,utm_x,utm_y )

     z_chk=0;  %Initialising the variable to 0 value.
%      number = 0;
     %1
     %Checking for (x,y) to lie in zone demarcated by [utm_x,utm_y].
     %Following is for a zone with  6 vertices. For ease of checking the zone is divided in two with 4 vertices each.
     %Similarly analysis can be done for a zone with 4 vertices 
     
     if (y-utm_y(1))<=(utm_y(2)-utm_y(1))/(utm_x(2)-utm_x(1))*(x-utm_x(1))
               
         if (y-utm_y(2))<=(utm_y(3)-utm_y(2))/(utm_x(3)-utm_x(2))*(x-utm_x(2))
                 
             if (y-utm_y(3))>=(utm_y(4)-utm_y(3))/(utm_x(4)-utm_x(3))*(x-utm_x(3)) 
               
                 if (y-utm_y(4))>=(utm_y(1)-utm_y(4))/(utm_x(1)-utm_x(4))*(x-utm_x(4))            
                                  
                     z_chk =1;    
                    
                 end
             end
         end
                          
    end
     
    
%       if (y-utm_y(3))<=(utm_y(4)-utm_y(3))/(utm_x(4)-utm_x(3))*(x-utm_x(3))
%               
%          if (y-utm_y(4))<=(utm_y(5)-utm_y(4))/(utm_x(5)-utm_x(4))*(x-utm_x(4))
%                
%              if (y-utm_y(5))>=(utm_y(6)-utm_y(5))/(utm_x(6)-utm_x(5))*(x-utm_x(5)) 
%                 
%                  if (y-utm_y(6))<=(utm_y(3)-utm_y(6))/(utm_x(3)-utm_x(6))*(x-utm_x(6))            
%                                   
%                      z_chk =1;
%                        
%                  end
%              end
%          end
%                           
%        end
%      
% end

