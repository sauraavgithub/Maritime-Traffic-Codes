function [ dcpa,tcpa ] = dcpa1(dist,c_rt,b_rt,v_rt)

%Calculating DCPA and TCPA

if(b_rt~=pi)
    dcpa=abs(10*dist*sin(c_rt-b_rt));
    tcpa=abs(60*0.5144*dist*cos(c_rt-b_rt)/v_rt);
else
    dcpa=10*dist;
    tcpa=0;
end

end

