function [ v_chk] = velocity_chk( sog)

if sog <= 0.1 % The sog limit can be set (here 0.1 knots) by the user. 
    v_chk = 0;
else
    v_chk = 1;
end

end

